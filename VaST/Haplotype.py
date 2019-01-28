import json
import logging
import numpy as np
import itertools as it
from utils import parse_flag_file
from analyze import (
    remove_extra_loci, get_resolution,
    get_summary_data, get_haplotype_matrix,
    get_resolution_matrix)


class Haplotype:
    def __init__(self, patterns, minimum_spanning_set,
                 flag_file_path, primer_zone_size, variant_matrix, sep):
        self._logger = logging.getLogger(__name__)
        self._variant_matrix = variant_matrix
        self._sep = sep
        self._minimum_spanning_set = minimum_spanning_set
        self._selected_patterns = \
            self._minimum_spanning_set.get_selected_patterns()
        self._selected_amplicons = \
            self._minimum_spanning_set.get_selected_amplicons()
        self._patterns = patterns
        self._pattern_dic = patterns.get_pattern_dic(
            self._selected_patterns)
        self._pattern_df = patterns.get_pattern_df(
            self._selected_patterns)
        # Selected amps is the group of amplicons for a pattern that
        # were not removed due to overlap with other amplicons
        # in the minimum spanning set.
        self._pattern_dic = self._get_selected_amplicons()
        self._get_flags(
            flag_file_path, int(primer_zone_size))

    def _get_flags(self, flag_file_path, primer_zone_size):
        flag_df = parse_flag_file(flag_file_path)
        for pattern, amplicons in self._pattern_dic.iteritems():
            for amplicon, chars in amplicons.iteritems():
                genome = chars['g']['name']
                genome_size = int(chars['g']['length'])
                start = int(amplicon)
                stop = int(chars['s'])
                up_start = start - primer_zone_size - 1 if start - primer_zone_size > 1 else 0
                up_stop = start - 1
                down_start = stop
                down_stop = (stop + primer_zone_size if 
                    stop + primer_zone_size < genome_size
                    else genome_size - 1)
                upstream_flags = np.array(
                    flag_df[flag_df.Genome == genome].iloc[up_start: up_stop].Flag, dtype=int)
                downstream_flags = np.array(
                    flag_df[flag_df.Genome == genome].iloc[down_start: down_stop].Flag, dtype=int)
                upstream_count = np.array([sum([not value for value in run[1]]) for
                    run in it.groupby(np.array(upstream_flags, dtype=bool))],
                    dtype=int)
                downstream_count = np.array(
                    [sum([not value for value in run[1]]) for
                    run in it.groupby(np.array(downstream_flags, dtype=bool))],
                    dtype=int)
                percent_ok = (
                    (np.sum(upstream_count) + np.sum(downstream_count))/float(
                        len(upstream_flags) + len(downstream_flags)) * 100)
                med_size = np.median(np.append(upstream_count, downstream_count))
                self._pattern_dic[pattern][amplicon]['primer_zone'] = {
                    'upstream': ",".join(np.array(upstream_flags, dtype=str)),
                    'downstream': ",".join(np.array(downstream_flags, dtype=str)),
                    'percent_ok': percent_ok,
                    'med_size': med_size
                }


    def _get_selected_amplicons(self):
        new_dic = {}
        for pattern, sel_amplicons in zip(
                self._selected_patterns, self._selected_amplicons):
            all_amplicons = self._pattern_dic[pattern]
            new_dic[pattern] = {
                k: v for k, v in all_amplicons.iteritems() if k in sel_amplicons}
        return new_dic

    def write_haplotype(self, file_name):
        self._logger.info("Writing haplotype to %s", file_name)
        self._pattern_df.to_csv(file_name)

    def write_json(self, file_name):
        self._logger.info(
            "Writing minimum spanning set amplicons to %s",
            file_name)
        with open(file_name, 'w') as out:
            out.write(json.dumps(self._pattern_dic))
    
    def write_haplotype_matrix(self, file_name):
        self._logger.info("Writing haplotype matrix to %s", file_name)


    def write_suggested_amplicons(self, file_name):
        self._logger.info(
            "Writing suggested amplicons to %s", file_name
        )
        with open(file_name, 'w') as out:
            out.write()

    def write_output(self, haplotype_output, pattern_output, amplicon_output):
        best_loci = remove_extra_loci(self._pattern_dic)
        pattern_order = self._pattern_df.columns
        haplotype_matrix = get_haplotype_matrix(
            pattern_order, best_loci, self._variant_matrix, self._sep)
        self._logger.info("Writing haplotype matrix to %s", haplotype_output)
        haplotype_matrix.to_csv(haplotype_output)
        self._logger.info("Writing pattern matrix to %s", pattern_output)
        scores, patterns = get_resolution(self._pattern_df)
        pattern_matrix = get_resolution_matrix(
            self._pattern_df.index, pattern_order, patterns)
        pattern_matrix.to_csv(pattern_output)
        self._logger.info("Writing amplicon matrix to %s", amplicon_output)
        amplicon_matrix = get_summary_data(best_loci, scores, pattern_order)
        amplicon_matrix.to_csv(amplicon_output, index=None)






    def write_summary(self, file_name):
        self._logger.info(
            "Writing summary to %s", file_name
        )
        with open(file_name, 'w') as out:
            out.write(
                "Minimum set size: {}\n".format(
                    len(self._selected_patterns)))
            out.write(
                "Resolution Index: {:0.2f}%\n".format(
                    self._minimum_spanning_set.get_resolution_index())
            )
            group_size, counts = \
                self._minimum_spanning_set.get_resolution_groups()
            out.write("Group Size Breakdown:\n")
            labels = [
                "Group(s) of size {}".format(i)
                if i > 1 else "Strain(s) Fully Resolved"
                for i in group_size]
            for label, count in zip(labels, counts):
                out.write("{0} {1}\n".format(count, label))
