
from __future__ import absolute_import, print_function, division
import argparse
import sys
import json
import logging
import os
import ast
import pandas as pd
import numpy as np
from glob import glob
from collections import Counter, defaultdict
from itertools import chain, product, starmap
from utils import file_type, path_type, config_logging


LOG = logging.getLogger(__name__)

def parse_args(args=None):
    parser = argparse.ArgumentParser(
        prog='VaST Analyze Results',
        description="Analyze and format pattern selection ouput",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "project_dir", metavar='PROJECT_DIR', type=path_type,
        help="Path to pattern selection directory"
    )

    parser.add_argument(
        "var_matrix_path", metavar="VAR_MATRIX_PATH", type=file_type,
        help="Path to variant site matrices."
    )

    parser.add_argument(
        "--sep", type=str, default="tab",
        choices=['tab', 'comma', 'space'],
        help="Specify the file delimiter of variant site matrix"
    )

    parser.add_argument(
        "--log", type=str, default="./VaST_analysis.log",
        help="Path to log file."
    )
    return vars(parser.parse_args(args))


def expand_files(project_dir, path_to_file):
    """
    Returns a list of all files in a given project dir with a given wildcard
    path.
    
    Input: 
        project_dir: Path to pattern selection project directory
        path_to_file: Subpath with wildcard file name
    Output:
        List of matching file paths
    """
    return glob(os.path.join(project_dir, path_to_file))

def get_file_from_project(project_dir, path_to_file):
    """
    Return file from pattern selection project.

    Input:
        project_dir: Path to pattern selection project directory
        path_to_file: Subpath with wildcard file name

    Output:
        Valid file path to project file

    Raises:
        ValueError: if multiple files match wildcard
        FileNotFoundError: if path is not a valid file
        IndexError: if no files match wildcard
    """
    expanded_files = expand_files(project_dir, path_to_file)
    if len(expanded_files) > 1:
        LOG.error(
            "Multiple input files found, directory "\
            "should only contain one: {}".format(
                expanded_files))
        raise ValueError("Multiple input files")
        
    try:
        if os.path.isfile(expanded_files[0]):
            return expanded_files[0]
        else:
            raise FileNotFoundError
    except IndexError:
        LOG.error("No input file found")
        raise

# TODO: replace this with pick_best_loci to pick better primer zones
def remove_extra_loci(patterns):
    """
        If a pattern has more than one possible loci, return only one
        Input:
            patterns: A dictionary of patterns
        Output: A dictionary of patterns with only a single loci per pattern.
    """
    return {pattern: {list(loci.keys())[0]: loci[list(loci.keys())[0]]}
            for pattern, loci in patterns.items()}

def pick_best_loci(patterns):
    """
    For patterns with multiple loci, pick the loci with the most
    conserved primer zones
    
    Input: 
        patterns: A dictionary of patterns
    Output: A dictionary with a single loci from each pattern
    """
    
    return {pattern: sorted(loci,
        key=lambda k: (loci[k]['primer_zone']['percent_ok'],
        loci[k]['primer_zone']['med_size']),
        reverse=True)[0] for pattern, loci in patterns.items()}

def get_ordered_patterns(order, patterns):
    """
    Place the sites in the pattern dictionary into the order that they were
    added. 
    Input:
        order: a list of pattern ids providing the order
        patterns: a dictionary of patterns
    Output: Nested array with columns for sites, their order, and amplicon id
    """
    order_dict = {pattern: order + 1 for order, pattern in enumerate(order)}
    result = []
    for pattern, loci in patterns.items():
        for locus, data in loci.items():
            for site in data['sites']:
                result.append([site.strip(), order_dict[pattern], pattern])
    return np.array(result, dtype=str)


def get_haplotype_matrix(pattern_order, pattern_json, var_matrix_path, sep):
    """
    Builds a dataframe with columns for site id, amplicon order, amplicon id,
    and for each genome in the SNP matrix which includes the base call for
    each genome for each site in the minimum spanning set.

    Input:
        pattern_order: a list containing the pattern ids in the order in which
                       they were added to the minimum spanning set.
        pattern_json: a dictionary of patterns
        var_matrix_path: path to the variant site matrix which contains the
                         base calls for all variant sites.
    Output: A dataframe
    """
    index_cols = get_ordered_patterns(pattern_order, pattern_json)
    var_matrix = pd.read_csv(
        var_matrix_path,
        index_col=0,
        sep=sep).loc[index_cols[:,0]]
    index_df = pd.DataFrame(
        index_cols[:,1:], index=index_cols[:,0], dtype=int,
        columns = ["Amplicon_Order", "Amplicon_ID"])
    return pd.merge(
        index_df, var_matrix, how="left",
        left_index=True, right_index=True).sort_values(by="Amplicon_Order")

def write_dataframe_to_file(df, outpath, **kwargs):
    """
    Write a pandas database as a csv to provided file path
    Input:
        df: pandas dataframe
        outpath: path to write csv
    Ouput: None

    Raises:
        FileNotFoundError: if not a valid path
    """
    try:
        df.to_csv(outpath, **kwargs)
    except FileNotFoundError:
        LOG.error("Cannot write to: {}".format(outpath))
        raise


def get_resolution(patterns_df):
    current = format_column(patterns_df.iloc[:,0])
    scores = []
    patterns = []
    for col in patterns_df:
        current = get_feature_categories(
            combine_patterns(current, format_column(patterns_df[col])))
        patterns.append([c for c in current])
        scores.append(get_resolution_score(current))
    n_strains = patterns_df.shape[0]
    scores = [
        (1 - (s / (n_strains**2 - n_strains))) * 100 for s in scores]
    return scores, patterns
        

def format_column(col):
    try:
        fmt = [ast.literal_eval(c) for c in col]
    except ValueError:
        return col



def get_resolution_score(pattern_vec):
    r = Counter(list(chain(*pattern_vec))).values()
    return sum([s**2 - s for s in r])

def combine_patterns(previous, current):
    return list(starmap(
        lambda x, y: list(set(product(x, y))),
        zip(previous, current)))

def get_feature_categories(combined_patterns):
    id_generator = counter()
    feature_dict = {}
    feature_list = []
    # Set up features for non-ambiguous strains
    for features in combined_patterns:
        if len(features) == 1:
            if features[0] not in feature_dict:
                feature_dict[features[0]] = next(id_generator)
    # Set up features that overlap in ambiguous strains
    ambiguous_features = defaultdict(list)
    for features in combined_patterns:   
        if len(features) > 1:
            for f, feature in enumerate(features):
                if feature not in feature_dict:
                    ambiguous_features[feature].append(f)
    for k, v in ambiguous_features.items():
        if len(v) > 1:
            feature_dict[k] = next(id_generator)
    for features in combined_patterns:
        current_features = []
        for feature in features:
            try:
                current_features.append(feature_dict[feature])
            except KeyError:
                pass
        # In case of completely unique ambiguous pattern
        if not current_features:
            current_features.append(next(id_generator))
        feature_list.append(tuple(set(current_features)))

    return feature_list

def get_summary_data(pattern_dict, resolution_scores, order):
    """
    Creates dataframe summarizing data in pattern_dict
    Input:
        pattern_dict: a dictionary of patterns
        resolution_scores: a list of resolution scores
                           in the same order as order
        order: the order of patterns
    Ouput: a dataframe
    """
    result = []
    for pattern, score in zip(order, resolution_scores):
        for locus, data in pattern_dict[pattern].items():
            result.append([
                pattern, # pattern id
                data['g']['name'], # genome name
                locus, # start position
                data['s'], # stop position
                int(data['s']) - int(locus) + 1, # size of target
                score, # Resolution score
                len(data['sites']), # number of sites
                " ".join([str(d) for d in data['sites']]), # list of sites
                data['primer_zone']['upstream'], # upstream primer zone
                data['primer_zone']['downstream']])
    return pd.DataFrame(result,
        columns=[
            "Amplicon_ID",
            "Sequence_ID",
            "Start_Position",
            "End_Position",
            "Target_Size",
            "Resolution_Score",
            "Num_of_Sites",
            "Sites",
            "Upstream_Primerzone",
            "Downstream_Primerzone"])

def counter(start=0):
    while True:
        yield start
        start += 1


def get_resolution_matrix(strains, pattern_ids, patterns):
    return pd.DataFrame(patterns, columns=strains, index=pattern_ids).T


def main(project_dir, var_matrix_path, log, sep):
    config_logging(log, "INFO")
    LOG.info("Starting Analysis")
    sep = {'comma': ',',
           'space': ' ',
           'tab': '\t'}[sep]
    pattern_json = json.loads(open(get_file_from_project(
            project_dir, "./minimum_spanning_set/amplicons_*.json"),
            'r').read())
    best_loci = remove_extra_loci(pattern_json)
    haplotype_df = pd.read_csv(get_file_from_project(
        project_dir, "./minimum_spanning_sehaplotype_*.csv"), index_col=0)
    pattern_order = list(haplotype_df.columns)
    write_dataframe_to_file(get_haplotype_matrix(
        pattern_order, best_loci, var_matrix_path, sep),
        os.path.join(project_dir,
            "minimum_spanning_set/haplotypes_matrix.csv"))
    scores, patterns = get_resolution(haplotype_df)
    write_dataframe_to_file(
        get_summary_data(best_loci, scores, pattern_order),
        os.path.join(project_dir,
        "minimum_spanning_set/amplicon_matrix.csv"), index=False)
    write_dataframe_to_file(
        get_resolution_matrix(
            haplotype_df.index, pattern_order, patterns),
        os.path.join(project_dir,
        "minimum_spanning_set/pattern_matrix.csv"))
    

        
if __name__ == "__main__":
    main(**parse_args(sys.argv[1:]))

