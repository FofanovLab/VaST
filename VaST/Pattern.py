import copy
import json
import logging
from functools import partial
from itertools import combinations

import numpy as np
import pandas as pd

from utils import cartesian_product, get_ambiguous_pattern


class Patterns:
    def __init__(self):
        self._patterns = {}
        self._genomes = {}
        self._amp_2_pattern = {}
        self._logger = logging.getLogger(__name__)
        self._strains = []
        self._resolution_levels = []
        self._pattern_df = None
        self._required_patterns = []

    def get_required_patterns(self):
        return self._required_patterns

    def get_pattern_dic(self, patterns):
        return {
            pattern: self._patterns[pattern]
            for pattern in patterns}

    def get_pattern_df(self, patterns):
        return self._pattern_df[patterns]

    def get_resolution_levels(self):
        return copy.deepcopy(self._resolution_levels)

    def to_json(self, file_name, strains):
        print self._genomes
        d = {"strains": strains,
             "genomes": {v['id']: v for v in self._genomes.values()},
             "patterns": self._patterns}
        self._logger.info("Writing patterns to: {}".format(file_name))
        with open(file_name, 'w') as out:
            out.write(json.dumps(d))

    def load_patterns(self, file_name):
        self._logger.info("Reading data from JSON: %s", file_name)
        try:
            with open(file_name, 'rU') as json_file:
                pattern_data = json.load(json_file)
        except ValueError:
            self._logger.error("Could not parse JSON: %s", file_name)
            raise
        self._strains = pattern_data['strains']
        self._genomes = pattern_data['genomes']
        self._patterns = pattern_data['patterns']
        self._get_pattern_dataframe()
        self._rekey_pattern_dict()
        self._get_amp_2_pattern()

    def add_required_sites(self, required_sites):
        self._logger.info(
            "Including the following sites in the solution: %s",
            ", ".join(required_sites))
        for site in required_sites:
            print "Site", required_sites
            print self._amp_2_pattern.keys()
            try:
                pattern_id = self._amp_2_pattern[site]['pattern_id']
                print "pattern id", pattern_id
                self._required_patterns.append(pattern_id)
            except KeyError:
                self._logger.warning("Required site not found: %s",
                                     site)
        self._required_patterns = list(set(self._required_patterns))

    def remove_sites(self, exclude_sites):
        self._logger.info(
            "Dropping the following sites: %s", ", ".join(exclude_sites))
        for site in exclude_sites:
            try:
                pattern_id = self._amp_2_pattern[site]['pattern_id']
                if len(self._patterns[pattern_id]) > 1:
                    # Only delete the amplicon if there are more
                    # than one amplicon for the pattern
                    del self._patterns[pattern_id][
                        self._amp_2_pattern[site]['site_id']]
                else:
                    # if this is the only amplicon for the pattern
                    # delete the pattern from dictionary and dataframe
                    del self._patterns[pattern_id]
                    self._pattern_df.drop(pattern_id, inplace=True)
                del self._amp_2_pattern[site]
            except KeyError:
                self._logger.warning(
                    "Excluded site not found: %s", site)

    def remove_strains(self, exclude_strains):
        """
        Removes strains from dataframe and strain list
        """
        self._logger.info("Removing strains")
        # check that all of excluded strains are in strains
        in_strains = np.intersect1d(
            self._strains, exclude_strains)
        if len(in_strains) != len(exclude_strains):
            self._logger.warning(
                "The following excluded strains were not "
                "found: %s",
                ", ".join(
                    np.setdiff1d(
                        exclude_strains, in_strains)))
        self._strains = np.setdiff1d(self._strains, exclude_strains)
        self._pattern_df = self._pattern_df[self._strains]
        self._logger.info(
            "%s strains are being used in analysis.",
            len(self._strains))
        self._logger.info(
            "The following strains were removed:\n%s",
            "\n".join(in_strains))

    def add_new_pattern(self, pattern, amplicon, genome_size):
        pattern_str = self._pattern_to_string(pattern)
        if amplicon.genome not in self._genomes:
            self._genomes[amplicon.genome] = {
                "id": len(self._genomes), "length": genome_size,
                "name": amplicon.genome}

        if pattern_str in self._patterns:
            self._patterns[pattern_str][amplicon.start] = {'s': amplicon.stop,
                                                           'g': self._genomes[amplicon.genome],
                                                           'sites': amplicon.site_ids}
        else:
            self._patterns[pattern_str] = {amplicon.start: {'s': amplicon.stop,
                                                            'g': self._genomes[amplicon.genome],
                                                            'sites': amplicon.site_ids}}

    def add_ambiguous_amplicon(self, features, amplicon, genome_size):
        features = [[tuple(r) for r in cartesian_product(row)]
                    for row in features]
        feature_categories = set([a[0] for a in features if len(a) == 1])
        feature_categories = [[a] for a in feature_categories]
        amb_pattern = partial(
            get_ambiguous_pattern,
            feature_categories=feature_categories)
        pattern = map(amb_pattern, features)
        self.add_new_pattern(pattern, amplicon, genome_size)

    def add_unambiguous_amplicon(self, features, amplicon, genome_size):
        # Requires numpy >= 1.13
        _, pattern = np.unique(features, return_inverse=True, axis=0)
        pattern = [(p,) for p in list(pattern)]
        self.add_new_pattern(pattern, amplicon, genome_size)

    def set_resolution(self, alt_resolution_file, stop_at_res):
        self._logger.info("Setting Resolution")
        full_resolution = self._get_full_resolution()
        if alt_resolution_file is not None:
            try:
                alt_resolution = pd.read_csv(
                    alt_resolution_file, index_col=0, header=None)
            except Exception:
                self._logger.error(
                    "Unable to parse resolution file: %s",
                    alt_resolution_file
                )
                raise
            # Make sure group ids are uniform
            for column in alt_resolution:
                alt_resolution[column] = np.unique(
                    alt_resolution[column],
                    return_inverse=True)[1]
            # Join full resolution and resolution dataframe
            # using a right join. Any NaNs indicate an invalid input.

            resolution = alt_resolution.join(
                full_resolution, how="right")
            if resolution.isnull().values.any():
                e_message = "Resolution file does not include all strains, "\
                    "NaNs produced."
                self._logger.error(e_message)
                raise ValueError(e_message)
            resolution = self._is_valid_resolution(resolution)

            # If stop at res, remove last column
            if stop_at_res:
                resolution = resolution.drop(
                    'Full_res',
                    axis=1, errors="ignore"
                )
            resolution.columns = [
                "Level_{}".format(i + 1) for i in range(
                    len(resolution.columns))]

            for col in resolution:
                self._resolution_levels.append(
                    Resolution_Pattern(
                        resolution[col],
                        self._get_copy_of_patterns()))
        else:
            self._resolution_levels.append(
                Resolution_Pattern(
                    full_resolution,
                    self._get_copy_of_patterns()
                )
            )

    def _get_pattern_dataframe(self):
        ''' Converts pattern dictionary into dataframe '''
        self._pattern_df = pd.DataFrame(index=self._strains)
        for i, pattern in enumerate(self._patterns.keys()):
            pattern = self._string_to_pattern(pattern)
            self._pattern_df[i] = pattern

    def _get_amp_2_pattern(self):
        ''' Reverse look up for amplicon to pattern '''
        for pattern_id, amps in self._patterns.iteritems():
            for site_id, amp in amps.iteritems():
                for site in amp['sites']:
                    self._amp_2_pattern[site] = {
                        'pattern_id': pattern_id,
                        'site_id': site_id}

    def _rekey_pattern_dict(self):
        rekey_dic = {}
        for i, amplicons in enumerate(self._patterns.values()):
            rekey_dic[i] = amplicons
        self._patterns = rekey_dic

    def _pattern_to_string(self, pattern):
        return ",".join(
            [" ".join(
                [str(p) for p in pp]) for pp in pattern])

    def _string_to_pattern(self, pat_str):
        return [tuple(
            [int(pp) for pp in p.split(" ")]) for p in pat_str.split(",")]

    def _get_copy_of_patterns(self):
        return self._pattern_df.copy()

    def _get_full_resolution(self):
        ''' return full resolution dataframe '''
        return pd.DataFrame(
            {"Full_res": range(len(self._strains))},
            index=self._strains)

    def _is_valid_resolution(self, resolution):
        # drop resolution columns that are the same
        drop_columns = []
        for c1, c2 in combinations(resolution.columns, 2):
            if same_group_pattern(
                    resolution[c1], resolution[c2], resolution.index):
                drop_columns.append(c2)
        resolution = resolution.drop(drop_columns, axis=1)
        if drop_columns:
            self._logger.warning(
                "Dropping duplicate resolution column in resolution file: %s",
                ", ".join(drop_columns))

        # Verify that columns are in increasing resolution order
        # the max group id should increase in each column
        if not resolution.max().equals(resolution.max().sort_values()):
            e_message = "Resolution is not increasing across columns, "\
                        "check resolution file"
            self._logger.error(e_message)
            raise ValueError(e_message)

        # Verify that the relationship is hierarchical meaning each
        # group is a subset of the same parent group
        for i, res_level in enumerate(resolution.columns[1:]):
            grouped = resolution.groupby([res_level])
            for name, group in grouped:
                if len(group.ix[:, i].unique()) > 1:
                    e_message = "Resolution is not hierarchical. "\
                                "Each group should be a subset of "\
                                "the same parent group. Issue in "\
                                "resolution level: %s", res_level
                    self._logger.error(e_message)
                    raise ValueError(e_message)
        return resolution


class Resolution_Pattern:
    def __init__(self, resolution, pattern_df):
        self._resolution = pd.DataFrame(resolution)
        self._pattern_df = pattern_df
        self._ambiguous = False
        self._set_resolution_pattern()

    def get_resolution_pattern(self):
        return self._pattern_df.copy().T

    def get_group_number(self):
        return int(self._resolution.max() + 1)

    def is_ambiguous(self):
        return self._ambiguous

    def _set_resolution_pattern(self):
        new_df = []
        temp_df = self._resolution.join(
            self._pattern_df)
        for name, group in temp_df.groupby(
                self._resolution.columns[0]):
            new_column = []
            for column in group[group.columns[1:]]:
                combined_pattern = ()
                for value in group[column]:
                    combined_pattern += value
                new_value = tuple(set(combined_pattern))
                # set ambiguous flag if multiple values in tuple
                if not self._ambiguous and len(new_value) > 1:
                    self._ambiguous = True
                new_column.append(new_value)
            new_df.append(new_column)
        self._pattern_df = pd.DataFrame(
            new_df, columns=self._pattern_df.columns, dtype=int)


def same_group_pattern(v1, v2, index):
    """Check if two arrays have the same pattern of similarity

    v1 and v2 and index are arrays of equal length. Each value in index is
    considered a label for the values in v1 and v2 in the corresponding
    position. Labels in index are grouped together if they share similar values
    in v1 and in v2. The sorted list of groups is compared to see if they are
    the equivalent.

    Arguments:
        v1 (array_like): First array (len(v1) = len(v2) = len(index))
        v2 (array_like): Second array
        index (array_like): Array of labels

    Returns:
        bool

    Raises:
        ValueError if each of the arguments are not the same length
    """
    if len(v1) == len(v2) == len(index):
        index = np.array(index)
        v1_groups = [
            ",".join(index[np.where(v1 == i)[0]]) for i in np.unique(v1)]
        v2_groups = [
            ",".join(index[np.where(v2 == i)[0]]) for i in np.unique(v2)]
        v1_groups.sort()
        v2_groups.sort()
        return v1_groups == v2_groups
    else:
        raise ValueError("Each array must be the same length")
