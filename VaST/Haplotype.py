import logging
import json
import pandas as pd
import numpy as np


class Haplotype:
    def __init__(self, patterns, minimum_spanning_set):
        self._minimum_spanning_set = minimum_spanning_set
        self._selected_patterns = \
            self._minimum_spanning_set.get_selected_patterns()
        self._patterns = patterns
        self._pattern_dic = patterns.get_pattern_dic(
            self._selected_patterns)
        self._pattern_df = patterns.get_pattern_df(
            self._selected_patterns)
        self._logger = logging.getLogger(__name__)


# TODO: Add summary
# TODO: Add upstream and downstream flags to amplicon dic


    def write_haplotype(self, file_name):
        self._logger.info("Writing haplotype to %s", file_name)
        self._pattern_df.to_csv(file_name)

    def write_json(self, file_name):
        self._logger.info(
            "Writing minimum spanning set amplicons to %s",
            file_name)
        with open(file_name, 'w') as out:
            out.write(json.dumps(self._pattern_dic))

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





