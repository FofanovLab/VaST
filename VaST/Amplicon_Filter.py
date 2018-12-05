import itertools as it
import logging
from collections import namedtuple

import numpy as np
import pandas as pd

from Pattern import Patterns
from utils import AMBIGUOUS_DICT

Amplicon = namedtuple('Amplicon', ['start', 'stop', 'genome', 'site_ids'])


class AmpliconFilter:
    def __init__(
            self, sites, var_matrix, flags, window, pz_size,
            pz_filter_length, pz_filter_percent, strict):
        self._patterns = Patterns()
        self._sites = pd.DataFrame(
            sites,
            columns=["Genome", "Start", "Stop"])
        self._sites['Start'] = self._sites['Start'].apply(pd.to_numeric)
        self._sites['Stop'] = self._sites['Stop'].apply(pd.to_numeric)
        self._sites = self._sites.groupby("Genome")
        self._var_matrix = var_matrix
        self._flags = flags
        self._window = window
        self._pz_size = pz_size
        self._pz_filter_length = pz_filter_length
        self._pz_filter_percent = pz_filter_percent
        self._logger = logging.getLogger(__name__)
        self._strict = strict

    def get_patterns(self):
        return self._patterns

    def filter_amplicons_get_patterns(self):
        self._logger.info("BEGIN Amplicon Filter")
        for genome, sites in self._sites:
            self._logger.info("Filtering sites in: %s", genome)
            sites = sites.sort_values("Start")
            pos = 0
            filtered_amp_count = 0
            accepted_amp_count = 0
            current_site = sites.iloc[pos].Start
            last_site = sites.iloc[-1].Start
            genome_size = len(self._flags[genome])
            go = True
            while go:
                amplicon = self._amplicon(genome, sites.loc[
                    (sites["Start"] >= current_site) &
                    (sites["Stop"] < current_site + self._window)],
                    genome_size)
                if amplicon is None:
                    pos += 1
                    filtered_amp_count += 1
                    if current_site == last_site:
                        go = False
                    else:
                        current_site = sites.iloc[pos].Start
                    continue
                self.update_pattern(amplicon, pos, genome_size)
                accepted_amp_count += 1
                if len(amplicon) == 1:
                    pos += 1
                else:
                    # check if shifting to the next site
                    # in the amplicon includes a new site
                    # if yes shift to next site in amplicon
                    # if no shift to next site outside of
                    # amplicon
                    next_site = amplicon.iloc[1].Start
                    next_amplicon = sites.loc[
                        (sites["Start"] >= next_site) &
                        (sites["Stop"] < next_site + self._window)]
                    if len(np.setdiff1d(next_amplicon.Start, amplicon.Start)):
                        pos += 1
                    else:
                        pos += len(amplicon)
                if pos >= len(sites):
                    go = False
                else:
                    current_site = sites.iloc[pos].Start

            self._logger.info("%s amplicon(s) filtered", filtered_amp_count)
            self._logger.info("%s amplicon(s) passed", accepted_amp_count )
        self._logger.info("FINISHED Amplicon Filter")
        return self._patterns

    def update_pattern(self, amplicon, position, genome_size):
        sites = self._var_matrix[
            position: position + len(amplicon)]
        amplicon = Amplicon(
            amplicon.iloc[0].Start,
            amplicon.iloc[-1].Stop,
            amplicon.iloc[0].Genome,
            list(amplicon.apply(
                lambda row: "::".join(
                    [str(r) for r in row]), axis=1)))

        if not self._strict:
            sites = _adjust_missing_if_vntr(sites)

            # TODO: Add check for integers in array and change
            # missing calls to any of the other
            sites = [[tuple(AMBIGUOUS_DICT[call]) if call in
                      AMBIGUOUS_DICT else tuple([call]) for call in site]
                     for site in sites]

            sites = map(list, zip(*sites))
            # Check if any sites got expanded
            if np.all(
                    np.equal(
                        [[len(site_i) for site_i in
                          site] for site in sites], 1)):
                self._patterns.add_unambiguous_amplicon(
                    sites, amplicon, genome_size)
            else:
                self._patterns.add_ambiguous_amplicon(
                    sites, amplicon, genome_size)

        else:
            self._patterns.add_unambiguous_amplicon(
                sites, amplicon, genome_size)
            # send directly to pattern

    def _amplicon(self, genome, amplicon, genome_size):
        if len(amplicon):
            start = amplicon.iloc[0].Start
            stop = amplicon.iloc[-1].Stop
            upstream, downstream = self._get_primer_zones(
                start, stop, genome, genome_size
            )
            pass_filter = self._filter(upstream, downstream)
            if pass_filter:
                return amplicon
            else:
                # do recursive call with last site removed
                return self._amplicon(
                    genome, amplicon.iloc[:-1], genome_size)
        else:
            return None

    def _get_primer_zones(
            self, amp_start, amp_stop, genome, genome_size):
        # subtract one from each index to get back to zeroindexing
        # TODO: Add option for circular chromosome to wrap around
        up_start = amp_start - self._pz_size - 1 if amp_start - self._pz_size > 1 else 0
        up_stop = amp_start - 1
        down_start = amp_stop
        down_stop = (
            amp_stop + self._pz_size if amp_stop +
            self._pz_size < genome_size else genome_size - 1)
        upstream_flags = np.array(
            self._flags[genome].iloc[up_start: up_stop].Flag, dtype=bool)
        downstream_flags = np.array(
            self._flags[genome].iloc[down_start: down_stop].Flag, dtype=bool)
        return upstream_flags, downstream_flags

    def _filter(self, upstream, downstream):
        upstream = np.array([sum(1 for _ in g[1])
                             for g in it.groupby(upstream) if np.all(g[0])], dtype=int)
        downstream = np.array([sum(1 for _ in g[1])
                               for g in it.groupby(downstream) if np.all(g[0])], dtype=int)
        # Dividing by pz_size even though targets near the beginning and
        # end may be shorter.
        upstream = np.divide(
            np.sum(
                upstream[upstream > self._pz_filter_length]),
            float(self._pz_size)) * 100
        downstream = np.divide(
            np.sum(
                downstream[downstream > self._pz_filter_length]),
            float(self._pz_size)) * 100
        if (upstream >= self._pz_filter_percent
                and downstream >= self._pz_filter_percent):
            return True
        else:
            return False


def _adjust_missing_if_vntr(sites):
    for site in sites:
        mask = [(s).isdigit() for s in site]
        if np.any(mask):
            site[site == "X"] = "*"
            AMBIGUOUS_DICT['*'] = list(np.unique(site[mask]))
    return sites
