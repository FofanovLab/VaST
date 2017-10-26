import logging
from collections import Counter
from functools import partial
from itertools import chain
from multiprocessing import Pool

import numpy as np
import pandas as pd
from ascii_graph import Pyasciigraph

from Pattern import Patterns, Resolution_Pattern
from utils import (cartesian_product, get_ambiguous_pattern, normalize,
                   parallel_apply)


class MinSet:
    def __init__(self, patterns, n_threads=1):

        self._logger = logging.getLogger(__name__)
        self._n_threads = n_threads
        self._patterns = patterns
        self._resolution_patterns = patterns.get_resolution_levels()
        self._max_resolution_score = self._get_max_resolution_score()
        self._current_is_ambiguous = None
        self._current_resolution_level = None
        self._update_patterns = None
        self._score_patterns = None
        self._selected_patterns = []
        self._selected_amplicons = []
        self._required_patterns = patterns.get_required_patterns()
        self._get_next_resolution_level()
        self._current_resolution_score = self._max_resolution_score
        self._current_resolution = self._init_resolution()
        self._resolution_index = self._get_resolution_index()
        self._pz_size = None

    def get_resolution_index(self):
        return self._resolution_index

    def get_resolution_groups(self):
        # Modify to work with ambiguous bases
        return np.unique(
            Counter(
                chain(
                    self._current_resolution)).values(),
            return_counts=True)

    def get_resolution_score(self):
        return self._current_resolution_score

    def get_selected_amplicons(self):
        return self._selected_amplicons

    def get_selected_patterns(self):
        return self._selected_patterns

    def get_minimum_spanning_set(
            self, max_loci, max_res, start, primer_zone_size):
        first = True
        go = True
        self._pz_size = primer_zone_size
        msg = "START Generating minimum spanning set {}".format(start)
        self._logger.info(msg)
        print msg
        while go:
            previous_resolution_score = self._current_resolution_score
            # score and pick minimum
            if self._required_patterns:
                next_pattern = self._required_patterns.pop(0)
                self._logger.info(
                    "Adding required pattern: %s", next_pattern)
            elif first:
                next_pattern = self._current_resolution_level.index[
                    start % len(self._current_resolution_level.index)]
                first = False
                self._logger.info(
                    "Adding first pattern starting with: %s",
                    next_pattern
                )
            else:
                next_pattern = self._get_next_pattern()
                self._logger.info(
                    "Adding pattern: %s", next_pattern
                )
            # if overlapping, remove and move to next pattern
            if self._overlapping_pattern(next_pattern):
                self._logger.info(
                    "Pattern removed due to overlapping amplicons"
                )
                self._remove_used_pattern(next_pattern)
                continue
            self._add_pattern_to_set([next_pattern])
            # check stopping points
            msg = "Starting Next Resolution Level"
            if previous_resolution_score == self._current_resolution_score:
                self._logger.info(
                    "No improvement with last pattern, removing pattern")
                # Check if there was no improvement with the last addition
                # Remove last addition
                self._selected_patterns = self._selected_patterns[:-1]
                self._selected_amplicons = self._selected_amplicons[:-1]
                if not self._get_next_resolution_level():
                    self._logger.info("Maximum resolution reached")
                    # stop if there is not another resolution level
                    go = False
                else:
                    print msg
                    self._logger.info(msg)
                continue
            if self._current_resolution_score == 0:
                self._logger.info(
                    "Maximum resolution achieved at current "
                    "resolution objective")
                # full resolution has been achieved
                if not self._get_next_resolution_level():
                    self._logger.info(
                        "Maximum Resolution achieved at final "
                        "resolution objective")
                    # stop if there is not another resolution level
                    go = False
                else:
                    print msg
                    self._logger.info(msg)
            if len(self._selected_patterns) == max_loci:
                self._logger.info("Maximum loci reached")
                # stop if max loci has been reached
                go = False
            if not len(self._resolution_patterns) and \
                    self._resolution_index >= max_res:
                self._logger.info(
                    "Maximum resolution percent reached")
                # check that we are on last resolution level
                # and resolution level has been reached
                go = False
            if not len(self._current_resolution_level):
                self._logger.info(
                    "No more patterns available"
                )
                go = False
            self._remove_used_pattern(next_pattern)
        return self._selected_patterns

    def _overlapping_pattern(self, pattern):
        # Current amplicons
        start, stop, genomes, ids = self._get_amplicon_positions(
            self._patterns.get_pattern_dic(
                self._selected_patterns))
        # New Amplicons
        new_start, new_stop, new_genomes, new_ids = self._get_amplicon_positions(
            self._patterns.get_pattern_dic([pattern]))
        if not start:  # First amplicon is being added
            self._selected_amplicons.append(new_ids[0])
            return False
        non_overlapping_amplicons = []
        for (s, e, g, i) in zip(new_start[0], new_stop[0], new_genomes[0], new_ids[0]):
            # TODO: make less stringent by removing
            # amplicon's from previously added patterns
            # if there are alternative amplicons
            if not np.any([np.any(a) for a in _get_overlap_array(
                    start, stop, genomes, s, e, g)]):
                non_overlapping_amplicons.append(i)
        if non_overlapping_amplicons:
            self._selected_amplicons.append(non_overlapping_amplicons)
            return False
        else:
            return True

    def _get_amplicon_positions(self, pattern_dic):
        start = []
        stop = []
        genomes = []
        ids = []
        for k, v in pattern_dic.iteritems():
            start.append([int(key) - self._pz_size for key in v.keys()])
            stop.append(
                [int(value['s']) + self._pz_size for value in v.values()])
            genomes.append([int(value['g']['id']) for value in v.values()])
            ids.append(v.keys())
        return start, stop, genomes, ids

    def _add_pattern_to_set(self, pattern_ids, get_progress=True):
        ''' update current resolution with new pattern '''
        for pattern_id in pattern_ids:
            self._current_resolution = \
                self._current_resolution_level.loc[pattern_id]
            self._current_resolution_score = \
                _ambiguous_score_helper(self._current_resolution)
            self._get_resolution_index()
            self._selected_patterns.append(pattern_id)
            if get_progress:
                self._get_progress()
            self._update_patterns()

    def _get_next_resolution_level(self):
        if len(self._resolution_patterns):
            next_level = self._resolution_patterns.pop(0)
            self._current_is_ambiguous = next_level.is_ambiguous()
            self._current_resolution_level = next_level.get_resolution_pattern()
            if self._current_is_ambiguous:
                self._logger.info(
                    "Current resolution objective is ambiguous")
                self._update_patterns = self._ambiguous_update_patterns
                self._score_patterns = self._ambiguous_score
            else:
                self._logger.info(
                    "Current resolution objective is unambiguous")
                self._update_patterns = self._unambiguous_update_patterns
                self._score_patterns = self._unambiguous_score
            self._sort_resolution_by_score()
            self._init_resolution()
            previous_selected_patterns = self._selected_patterns
            self._selected_patterns = []
            self._add_pattern_to_set(
                previous_selected_patterns, get_progress=False)
            return True
        else:
            return False

    def _sort_resolution_by_score(self):
        self._current_resolution_level['Scores'] = self._score_patterns()
        self._current_resolution_level = \
            self._current_resolution_level.sort_values('Scores').drop(
                'Scores', axis=1)

    def _init_resolution(self):
        return pd.DataFrame(
            {0: [(0,) for i in xrange(
                len(self._current_resolution_level.columns))]},
            index=self._current_resolution_level.columns)

    def _update(self, update_helper_func):
        self._current_resolution_level = parallel_apply(
            self._current_resolution_level,
            update_helper_func,
            self._n_threads, 1,
            **{"current_vector":
                self._current_resolution.values}
        )

    def _ambiguous_update_patterns(self):
        self._update(_ambiguous_update_helper)

    def _unambiguous_update_patterns(self):
        self._update(_unambiguous_update_helper)

    def _get_next_pattern(self):
        return self._score_patterns().idxmin()
        # return scores[scores == scores.min()].sample().index[0]

    def _ambiguous_score(self):
        return self._score(_ambiguous_score_helper)

    def _unambiguous_score(self):
        return self._score(_unambiguous_score_helper)

    def _score(self, score_helper_func):
        score = parallel_apply(
            self._current_resolution_level,
            score_helper_func,
            self._n_threads, 1)
        return score

    def _remove_used_pattern(self, pattern_id):
        self._current_resolution_level.drop(
            pattern_id, axis=0, inplace=True)

    def _get_max_resolution_score(self):
        n = self._resolution_patterns[-1].get_group_number()
        return n * (n - 1)

    def _get_resolution_index(self):
        self._resolution_index = (1 - normalize(
            self._current_resolution_score,
            self._max_resolution_score)) * 100

    def _get_progress(self):
        graph = Pyasciigraph()
        data, counts = self.get_resolution_groups()
        labels = ["Group size {}".format(i)
                  if i > 1 else "Fully Resolved" for i in data]
        progress_string = ""
        for line in graph.graph(
                "Pattern # {}".format(len(self._selected_patterns)),
                zip(labels, counts)):
            line = line.replace("#", "=")
            progress_string += line + "\n"
        percent_complete = self._resolution_index / 100.0
        ascii_percent = int(55 * percent_complete)
        progress_string += "=" * 79 + \
            "\nPercent Complete {0:0.1f}% ".format(percent_complete*100)
        progress_string += "|" * ascii_percent
        progress_string += "-" * (55 - ascii_percent)
        progress_string += "|\n" + "=" * 79 + "\n\n"
        print progress_string.encode('utf-8')
        self._logger.info(progress_string)


def _cartesian_product(pattern_vector, current_vector):
    cart_prod = partial(cartesian_product, dtype=int)
    p_update = map(
        cart_prod,
        [[p, c] for p, c in zip(
            pattern_vector, current_vector)])
    return [[tuple(r) for r in rr] for rr in p_update]


def _ambiguous_update_helper(pattern_vector, current_vector):
    p_update = _cartesian_product(
        pattern_vector, current_vector)

    a = set([p[0] for p in p_update if len(p) == 1])
    a = [[i] for i in a]
    amb_pattern = partial(
        get_ambiguous_pattern,
        feature_categories=a)
    return map(amb_pattern, p_update)


def _unambiguous_update_helper(pattern_vector, current_vector):
    p_update = _cartesian_product(
        pattern_vector, current_vector)
    _, pattern = np.unique(p_update, return_inverse=True, axis=0)
    return [(p,) for p in list(pattern)]


def _ambiguous_score_helper(pattern):
    n_strains_in_feature_category = Counter(chain(*pattern))
    return sum([s * s - s for s in
                n_strains_in_feature_category.values()])


def _unambiguous_score_helper(pattern):
    _, n_strains_in_feature_category = np.unique(
        pattern, return_counts=True)
    return sum(
        [s * s - s for s in
         n_strains_in_feature_category])


def _get_overlap_array(
        current_start, current_stop, current_genome,
        new_start, new_stop, new_genome):
    overlap_array = []
    for start, stop, genome in zip(
            current_start, current_stop, current_genome):
        overlap_array.append(
            _check_for_overlap(
                start, stop, genome, new_start, new_stop, new_genome))
    return overlap_array


def _check_for_overlap(
        start_1, stop_1, genome_1, start_2, stop_2, genome_2):
    return np.logical_and.reduce((
        np.equal(genome_1, genome_2),
        np.less_equal(start_1, stop_2),
        np.less_equal(start_2, stop_1)))
