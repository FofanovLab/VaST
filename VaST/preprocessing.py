import argparse
import logging
import os
import time
from multiprocessing import Pool

import numpy as np
import pandas as pd

from Amplicon_Filter import AmpliconFilter
from utils import (AMBIGUOUS, BASES, History, Project_Directory,
                   config_logging, file_type, path, percent, positive_int,
                   read_list_file)


def _get_strains_from_file(matrix, sep="\t"):
    """Get strain names from file header

    Arguments:
    matrix (str): full path and file name of matrix

    Keyword Arguments:
    sep (str): the delimiter character used in the file
                (default: "\t")

    Returns:
    Numpy str array of strain names

    """
    logger = logging.getLogger(__name__)
    logger.info(
        "Reading strains from file: %s",
        matrix)
    try:
        with open(matrix, 'rU') as infile:
            header = infile.readline().strip().split(sep)
    except IOError:
        logger.error(
            "Cannot open matrix file: %s", matrix)
        raise
    strains = np.array(header[1:], dtype=str)
    logger.info("Read %s strains from file", len(strains))
    logger.debug(
        "Strains:\n%s", "\n".join(strains))
    return strains


def _get_sites_from_file(matrix, sep="\t"):
    """Get sites from first column in matrix

    Arguments:
    matrix (str): full path and file name of matrix

    Keyword Arguments:
    sep (str): the delimiter character used in the file
                (default: "\t")

    Returns:
    Numpy tuple array of sites split by "::"

    Raises:
    IOError if file cannot be opened

    """
    logger = logging.getLogger(__name__)
    logger.info(
        "Reading sites from matrix file: %s",
        matrix)
    try:
        with open(matrix, 'rU') as matrix_handle:
            sites = np.genfromtxt(
                matrix_handle, usecols=0, delimiter=sep,
                dtype=str, skip_header=True)
    except IOError:
        logger.error(
            "Cannot open matrix file: %s", matrix)
        raise
    sites = np.array(sites, dtype=str)
    logger.info("Read %s sites from file", len(sites))
    sites = np.array([s.split("::") for s in sites], dtype=str)
    return sites


def _check_inputs(pz_size, pz_filter_length,
                  strain_cutoff, n_strains):
    logger = logging.getLogger(__name__)
    logger.info("Checking input parameters")
    # Raise error if the size of the primer fragment
    # cutoff is larger than size of primer zone,
    # this would cause all of the variable regions
    # to be filtered.
    try:
        if pz_filter_length > pz_size:
            raise ValueError(
                "The pz_filter_length is larger "
                "than the primer zone")
    except ValueError:
        logger.error("", exc_info=True)
        raise

    if strain_cutoff > n_strains:
        logger.warning("The core cutoff if larger than the "
                       "number of strains being analyzed")
    logger.info("Done checking input parameters")


def _remove_strains(exclude_strains, strains):
    """
    Returns array of strains that should be included
    in analysis
    """
    logger = logging.getLogger(__name__)
    logger.info("Removing excluded strains from analysis")
    try:
        exclude = read_file_list(exclude_strains)
    except IOError:
        logger.error(
            "Cannot open excluded strain file: %s",
            exclude_strains)
        raise
    # check that all of excluded strains are in strains
    in_strains = np.intersect1d(strains, exclude)
    if len(in_strains) != len(exclude):
        logger.warning(
            "The following strains were not "
            "found in the variant matrix: %s",
            ", ".join(
                np.setdiff1d(
                    exclude, in_strains)))
    strains = np.setdiff1d(strains, exclude)
    logger.info(
        "%s strains are being used in analysis.",
        len(strains))
    logger.debug(
        "The following strains were removed:\n%s",
        "\n".join(in_strains))
    logger.debug(
        "The following strains remain:\n%s",
        "\n".join(strains))
    return strains


def _set_parameters(**kwargs):
    for key, value in [("strict", False),
                       ("skip_filter", False),
                       ("sep", "\t"),
                       ("window", 50),
                       ("strain_cutoff", 1),
                       ("pz_size", 200),
                       ("pz_filter_percent", 25),
                       ("pz_filter_length", 25),
                       ("exclude_strains", None),
                       ("n_threads", 1)]:
        if key not in kwargs:
            kwargs[key] = value
        # Important for when this parameter is used
        # for checking overlap in pattern_selection
        if kwargs['skip_filter']:
            kwargs['pz_size'] = 0
    return kwargs


def _parse_var_matrix(var_matrix_path, strains, sep):
    return pd.read_csv(
        var_matrix_path, sep=sep,
        usecols=strains, dtype=str).apply(
            lambda x: x.str.upper(), axis=0).as_matrix()


def _get_flags_helper(chunk):
    logger = logging.getLogger(__name__)
    logger.info(
        "On chunk: %s - %s", chunk.index.min(),
        chunk.index.max())
    chunk = chunk.apply(
        lambda x: x.str.upper(), axis=1)
    return chunk.apply(
        _non_conserved_strain_count, axis=1)


def _get_flags_from_full_matrix(
        full_matrix_path, strains, strain_cutoff,
        sep, n_threads, outfile):
    ''' 
    Get dictionary of flagged primer zone regions 
    Returns a dictionary key: Genome ID 
    value: dataframe of flags for each site
    '''
    logger = logging.getLogger(__name__)
    logger.info("BEGIN parsing full genome matrix")
    full_matrix_strains = _get_strains_from_file(
        full_matrix_path, sep)
    full_matrix_sites = _get_sites_from_file(
        full_matrix_path, sep)
    # Check that all of the strains are in full strain matrix
    missing_strains = np.setdiff1d(
        strains, full_matrix_strains)
    if len(missing_strains):
        err = "Missing strains in full matrix: {}".format(
            ", ".join(missing_strains))
        logger.error(err)
        raise IOError(err)

    chunk_size = len(
        full_matrix_sites) / n_threads if len(
            full_matrix_sites) / n_threads < 50000 else 50000
    full_matrix = pd.read_csv(
        full_matrix_path, delimiter=sep, usecols=strains,
        iterator=True, chunksize=chunk_size, dtype=str)

    # For now assume that all sites are there
    # TODO: fill in missing sites with flags and warn

    pool = Pool(n_threads)
    flags = pd.concat(pool.map(_get_flags_helper, full_matrix))

    # flags = np.array([])
    # for c, chunk in enumerate(full_matrix):
    #     logger.debug("On chunk number: %s", c)
    #     chunk = chunk.apply(
    #     lambda x: x.str.upper(), axis=1)
    #     flags = np.append(
    #         flags, np.array(
    #             chunk.apply(
    #                 _non_conserved_strain_count, axis=1),
    #                 dtype=int))

    # Adjust for strain cutoff
    flags = np.less(flags, strain_cutoff)
    # Split into dictionary by genome
    # TODO: Add check that no sites are missing
    flags = pd.DataFrame({"Genome": full_matrix_sites[:, 0],
                          "Site": full_matrix_sites[:, 1],
                          "Flag": flags})
    logger.info("Writing flags to %s", outfile)
    flags[['Site', 'Genome', 'Flag']].to_csv(outfile, index=False)

    flag_dic = {}
    for name, group in flags.groupby('Genome'):
        flag_dic[name] = group
    return flag_dic


def _non_conserved_strain_count(call_data):
    unique, counts = np.unique(
        call_data.str.upper(), return_counts=True)
    if len(unique) is 1 and unique[0] in BASES:
        return 0
    flag_count = 0
    snp_list = []
    for u, c in zip(unique, counts):
        if u in AMBIGUOUS:
            flag_count += c
        else:
            snp_list.append(c)
    if len(snp_list) > 1:
        # Count minority snps
        flag_count += sum(snp_list) - max(snp_list)
    return flag_count


def _get_unambiguous_sites(calls):
    return not np.any(
        map(
            lambda call: [a in call for a in AMBIGUOUS],
            np.unique(calls)))


def _remove_ambiguous_sites(var_matrix, sites):
    unambiguous_sites = np.array(
        map(_get_unambiguous_sites, var_matrix),
        dtype=bool)
    return (var_matrix[unambiguous_sites],
            sites[unambiguous_sites])


def preprocessing(project_directory, var_matrix_path,
                  full_matrix_path, **kwargs):
    """
    Groups variant sites into amplicon windows and filters
    out any amplicons that do not have well conserved
    upstream and downstream primer regions
    """
    logger = logging.getLogger(__name__)
    logger.info("BEGIN Preprocessing")
    args = _set_parameters(**kwargs)
    start_time = time.time()

    history = History(
        project_directory.make_new_file(
            "history", "preprocessing_history"),
        "Preprocessing",
        project_directory.timestamp,
        param_dict=args)

    # Get strains and sites from matrix
    strains = _get_strains_from_file(var_matrix_path)
    sites = _get_sites_from_file(var_matrix_path)

    # Remove excluded strains
    if args["exclude_strains"] is not None:
        strains = _remove_strains(args["exclude_strains"],
                                  strains)

    var_matrix = _parse_var_matrix(
        var_matrix_path, strains, args["sep"])

    if args["strict"]:
        n_sites_before = len(sites)
        var_matrix, sites = _remove_ambiguous_sites(
            var_matrix, sites)
        logger.info("Strict Mode: {} sites with ambiguous "
                    "or missing data were removed".format(
                        n_sites_before - len(sites)))

    history.add_path("Variant Site Matrix File", var_matrix_path)
    history.add_path("Full Genome Matrix File", full_matrix_path)
    history.add_parameter("Number of Sites", len(sites))
    history.add_parameter("Number of Strains", len(strains))

    _check_inputs(args["pz_size"],
                  args["pz_filter_length"],
                  args["strain_cutoff"],
                  len(strains))
    # TODO: save flag dict as JSON that can be provided as an option
    # instead of the full_matrix
    flag_file = project_directory.make_new_file(
        "flags", "primer_zone_flags", "csv")
    history.add_path("Primer Zone Flags", flag_file)

    flag_dic = {}
    if not args['skip_filter']:
        flag_dic = _get_flags_from_full_matrix(
            full_matrix_path, strains,
            args["strain_cutoff"], args["sep"],
            args["n_threads"],
            flag_file)

    amplicon_filter = AmpliconFilter(
        sites, var_matrix, flag_dic,
        args['window'], args['pz_size'],
        args['pz_filter_length'],
        args['pz_filter_percent'],
        args['strict'])

    if not args['skip_filter']:
        patterns = amplicon_filter.filter_amplicons_get_patterns()

    else:
        sites = pd.DataFrame(
            sites,
            columns=["Genome", "Start", "Stop"])
        for site in xrange(len(sites)):
            amplicon_filter.update_pattern(
                sites, site)
        patterns = amplicon_filter.get_patterns()
    # Write patterns to a json file
    pattern_json_file = project_directory.make_new_file(
        "patterns", "patterns", "json")
    patterns.to_json(pattern_json_file, list(strains))
    history.add_path("PATTERN JSON", pattern_json_file)

    # Write history
    logger.info("FINISHED Preprocessing")
    run_time = time.time() - start_time
    history.add_other("Run Time", run_time)
    history.write()

    # TODO: Write summary file
    # total number of amplicons, how many genomes
    # Average number of sites per amplicon
    # Number of amplicons with same pattern etc.


if __name__ == "__main__":

    PARSER = argparse.ArgumentParser(
        prog='AMPLICON_FILTER',
        description="Parses a variant site matrix to find "
                    "candidate amplicons",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    PARSER.add_argument(
        "project_name", metavar='PROJECT_NAME', type=str,
        help="Project name and output file prefix"
    )

    PARSER.add_argument(
        "project_dir", metavar='PROJECT_DIR', type=path,
        help="Project directory"
    )

    PARSER.add_argument(
        "var_matrix_path", metavar="VAR_MATRIX_PATH", type=file_type,
        help="Path to variant site matrices"
    )

    PARSER.add_argument(
        "full_matrix_path", metavar="FULL_MATRIX_PATH", type=file_type,
        help="Path to full genome matrix (includes calls at every site "
             "in the genome)"
    )

    PARSER.add_argument(
        '-d', "--sep", type=str, default="tab",
        choices=['tab', 'comma', 'space'],
        help="Specify the file delimiter of variant site matrix"
    )

    PARSER.add_argument(
        '-w', "--window", type=positive_int, default=50,
        help="Size of window for variable region. "
             "(HINT: Should be small enough to fit in desired "
             "amplicon with room for upstream and downstream "
             "primer zone)"
    )

    PARSER.add_argument(
        '-s', "--strict", action='store_true', default=False,
        help="Strict mode ignores sites with missing or ambiguous data"
    )

    PARSER.add_argument(
        "--skip_filter", action='store_true', default=False,
        help="Skip amplicon filter and just run pattern discovery"
    )

    PARSER.add_argument(
        '-c', "--strain_cutoff", type=positive_int, default=1,
        help="The number of strains at a primer zone site that can have "
             "a non-conserved call before the site is flagged"
    )

    PARSER.add_argument(
        '-z', "--pz_size", type=positive_int, default=200,
        help="Size of the primer zone"
    )

    PARSER.add_argument(
        '-l', "--pz_filter_length", type=positive_int, default=25,
        help="The length of un-flagged primer zone segments that count "
             "toward the primer zone filter percent"
    )

    PARSER.add_argument(
        '-p', "--pz_filter_percent", type=percent, default=25,
        help="The percent of primer zone positions that must be present "
             "in un-flagged segments of the primer zone that are longer "
             "than the primer zone filter length"
    )

    PARSER.add_argument(
        '--log', type=str, default="INFO",
        choices=['DEBUG', 'INFO', 'WARNING'],
        help="Set logging level")

    PARSER.add_argument(
        "-x", "--exclude_strains", type=file, default=None,
        help="Path to file containing a list of strains to exclude "
             "from analysis in a single column"
    )

    PARSER.add_argument(
        "-t", "--threads", type=positive_int, default=1,
        help="Number of threads for multiprocessing"
    )

    ARGS = PARSER.parse_args()

    # Make project directory instance
    # Set up project tree
    # Project_Name_Y_M_D_H_M_S
    # ---- patterns
    # ---- history
    # ---- logs
    PROJECT_DIR = Project_Directory(
        ARGS.project_dir, ARGS.project_name,
        ["patterns", "logs", "history", "flags"])

    config_logging(
        os.path.join(
            PROJECT_DIR.get_sub_directory("logs"),
            "{0}_preprocessing.log".format(
                PROJECT_DIR.project_name)), ARGS.log)

    SEP = {'comma': ',',
           'space': ' ',
           'tab': '\t'}[ARGS.sep]

    PARAMS = {"strict": ARGS.strict,
              "skip_filter": ARGS.skip_filter,
              "window": ARGS.window,
              "strain_cutoff": ARGS.strain_cutoff,
              "pz_size": ARGS.pz_size,
              "pz_filter_percent": ARGS.pz_filter_percent,
              "pz_filter_length": ARGS.pz_filter_length,
              "exclude_strains": ARGS.exclude_strains,
              "sep": SEP,
              "n_threads": ARGS.threads}

    preprocessing(
        PROJECT_DIR, ARGS.var_matrix_path,
        ARGS.full_matrix_path, **PARAMS)
