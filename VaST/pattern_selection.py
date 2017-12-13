import argparse
import logging
import os.path
import time

import numpy as np

from Haplotype import Haplotype
from Minimum_Spanning_Set import MinSet
from Pattern import Patterns
from utils import (History, Project_Directory, config_logging, file_type, path,
                   percent, positive_int, read_list_file)


def _get_minimum_spanning_set(
        patterns, reps, max_loci, max_res, n_threads, primer_zone_size):
    minimum_spanning_set_objs = [
        MinSet(
            patterns, n_threads)
        for _ in xrange(reps)]
    mss_size = [len(mss.get_minimum_spanning_set(
        max_loci, max_res, rep, primer_zone_size))
        for rep, mss in enumerate(minimum_spanning_set_objs)]
    res_scores = [
        mss.get_resolution_score() for mss in minimum_spanning_set_objs]
    return sorted(
        zip(minimum_spanning_set_objs, res_scores, mss_size),
        key=lambda x: (x[1], x[2]))[0][0]


def _check_inputs(max_loci, required_loci, exclude_loci):
    logger = logging.getLogger(__name__)
    # If there are required loci, check that max_loci is larger
    try:
        if max_loci < len(required_loci):
            raise ValueError("Maximum number of loci is less than or equal "
                             "to the number of required loci. "
                             "No more loci can be added")
    except ValueError:
        logger.error("", exc_info=True)
        raise

    # check that include loci and exclude loci are mutually exclusive
    try:
        intersect = np.intersect1d(exclude_loci, required_loci)
        if len(intersect):
            raise ValueError("Excluded loci and required loci are not "
                             "mutually exclusive: {}".format(
                                 ", ".join(intersect)))
    except ValueError:
        logger.error("", exc_info=True)
        raise


def _set_parameters(**kwargs):
    for key, value in [("res", None),
                       ("stop_at_res", False),
                       ("n_threads", 1),
                       ("max_loci", np.inf),
                       ("max_res", 100),
                       ("reps", 1),
                       ("exclude_loci", []),
                       ("required_loci", []),
                       ("exclude_strains", [])]:
        if key not in kwargs:
            kwargs[key] = value
    return kwargs


def pattern_selection(project_directory, **kwargs):
    logger = logging.getLogger(__name__)
    logger.info("BEGIN Pattern Selection")
    args = _set_parameters(**kwargs)
    start_time = time.time()
    _check_inputs(args['max_loci'],
                  args['required_loci'],
                  args['exclude_loci'])
    history = History(
        project_directory.make_new_file(
            "history", "pattern_selection_history"),
        "Pattern_Selection",
        project_directory.timestamp,
        param_dict=args)

    preprocessing_history = History(
        project_directory.get_parent_subdirectory_file(
            "history",
            "preprocessing_history_{}.txt".format(
                project_directory.get_parent_directory_timestamp())),
        "Preprocessing",
        exists=True)

    # Get JSON file path from preprocessing step
    json_file = preprocessing_history.get_path("PATTERN JSON")

    # Get flag file path from preprocessing step
    flag_file = preprocessing_history.get_path("PRIMER ZONE FLAGS")
    primer_zone_size = preprocessing_history.get_parameter("PZ_SIZE")


    history.add_path("PATTERN JSON", json_file)
    logger.info("Reading from pattern JSON: %s", json_file)
    # Read in pattern JSON
    patterns = Patterns()
    patterns.load_patterns(json_file)
    if len(args['exclude_loci']):
        patterns.remove_sites(args['exclude_loci'])
    if len(args['required_loci']):
        patterns.add_required_sites(args['required_loci'])
    if len(args['exclude_strains']):
        patterns.remove_strains(args['exclude_strains'])
    patterns.set_resolution(args['res'], args['stop_at_res'])
    best_set = _get_minimum_spanning_set(
        patterns, args['reps'], args['max_loci'],
        args['max_res'], args['n_threads'],
        int(preprocessing_history.get_parameter("PZ_SIZE")))

    haplotype_file = project_directory.make_new_file(
        "minimum_spanning_set", "haplotype", "csv")
    amplicon_json = project_directory.make_new_file(
        "minimum_spanning_set", "amplicons", "json")
    summary_file = project_directory.make_new_file(
        "summary", "summary")

    haplotype = Haplotype(patterns, best_set, flag_file, primer_zone_size)

    haplotype.write_haplotype(haplotype_file)
    history.add_path("Haplotype File", haplotype_file)

    haplotype.write_json(amplicon_json)
    history.add_path("Amplicon JSON", amplicon_json)

    haplotype.write_summary(summary_file)
    history.add_path("Summary", summary_file)

    logger.info("FINISHED Pattern Selection")
    run_time = time.time() - start_time
    history.add_other("Run Time", run_time)
    history.write()


if __name__ == "__main__":

    PARSER = argparse.ArgumentParser(
        prog='PATTERN_SELECTION',
        description="Find the minimum set of patterns that "
                    "uniquely identify all strains",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    PARSER.add_argument(
        "project_name", metavar='PROJECT_NAME', type=str,
        help="Project name and output file prefix"
    )

    PARSER.add_argument(
        "project_dir", metavar='PROJECT_DIR', type=path,
        help="Path to preprocessing directory"
    )

    PARSER.add_argument(
        "--res", type=file_type, default=None,
        help="Path to file which defines the desired level of resolution for "
             "each of the strains (the default if full strain resolution) "
             "The resolution file should be a csv table "
             "with each row starting with a strain name (exactly as it "
             "appears in the variant site matrix) followed by a group ID. Strains "
             "with the same group ID will not be resolved. Multiple columns "
             "of group IDs may be present but must start with the lowest "
             "resolution grouping and increase to higher resolution groupings "
             "and must be consistent."
    )

    PARSER.add_argument(
        "--stop_at_res", action="store_true",
        help="Force solution to stop once the maximum level specified in "
             "resolution file has been reached. Otherwise, the solution will "
             "continue towards full strain resolution once it has achieved "
             "the level of resolution defined in the file. It will only "
             "terminate early when MAX_LOCI is reached or no more resolution "
             "is achievable."
    )

    PARSER.add_argument(
        '--log', type=str, default="INFO",
        choices=['DEBUG', 'INFO', 'WARNING'],
        help="Set logging level")

    PARSER.add_argument(
        '-t', '--threads', type=int, default=1,
        help="Set number of threads"
    )

    PARSER.add_argument(
        '--reps', type=positive_int, default=1,
        help="Run N sets at a time and pick the best"
    )

    EARLY_TERM = PARSER.add_argument_group(
        title="Early Termination",
        description="Options for creating early stopping points"
    )

    EARLY_GROUP = EARLY_TERM.add_mutually_exclusive_group(required=False)
    EARLY_GROUP.add_argument(
        "--max_loci", type=positive_int, default=np.inf,
        help="Stop after MAX_LOCI is reached. If required loci are "
             "added, MAX_LOCI must be greater than or equal to the "
             "number of required loci."
    )

    EARLY_GROUP.add_argument(
        "--max_res", type=percent, default=100,
        help="Stop when the resolution has reached MAX_RES percent "
             "of total resolution."
    )

    REQ_LOCI_GROUP = PARSER.add_argument_group(
        title="Add required loci to solution",
        description="Add loci that must be included in the solution "
                    "either from command line or from file. Names should "
                    "match the locus ID given in the variant site matrix."
    )

    REQ_GROUP = REQ_LOCI_GROUP.add_mutually_exclusive_group(required=False)
    REQ_GROUP.add_argument(
        "--required_loci", nargs='*', type=str, default=[],
        help="List of locus IDs that must be included in solution."
    )

    REQ_GROUP.add_argument(
        "--req_loci_file", type=file_type, default=None,
        help="Path to file containing list of locus IDs "
             "(in a single column) that must be included "
             "in solution."
    )

    EXC_LOCI_GROUP = PARSER.add_argument_group(
        title="Exclude loci from consideration",
        description="Add loci either from command line or from file. "
                    "These loci will removed from consideration when "
                    "building the solution. Names should match "
                    "the locus ID given in the variant site matrix."
    )

    EXC_GROUP = EXC_LOCI_GROUP.add_mutually_exclusive_group()
    EXC_GROUP.add_argument(
        "--exclude_loci", nargs='*', type=str, default=[],
        help="List of loci that must be excluded from solution"
    )

    EXC_GROUP.add_argument(
        "--excl_loci_file", type=file_type, default=None,
        help="Path to file containing list of locus IDs "
             "(in a single column) that must be excluded "
             "from solution."
    )

    EXC_STRAIN_GROUP = PARSER.add_argument_group(
        title="Exclude strains from consideration",
        description="Add strains either from command line or from file. "
                    "Names should match the variant site matrix."
    )

    EXC_STRAINS = EXC_STRAIN_GROUP.add_mutually_exclusive_group()
    EXC_STRAINS.add_argument(
        "--exclude_strains", nargs='*', type=str, default=[],
        help="List of strains to be excluded"
    )

    EXC_STRAINS.add_argument(
        "--excl_strains_file", type=file_type, default=None,
        help="Path to file containing list of strains "
             "(in a single column) that should be excluded"
    )

    ARGS = PARSER.parse_args()

    PROJECT_DIR = Project_Directory(
        ARGS.project_dir,
        ARGS.project_name,
        ["summary", "logs", "history", "minimum_spanning_set"],
        ["patterns", "flags", "history"])

    config_logging(
        os.path.join(
            PROJECT_DIR.get_sub_directory("logs"),
            "{0}_pattern_selection.log".format(
                PROJECT_DIR.project_name)), ARGS.log)

    LOGGER = logging.getLogger(__name__)

    try:
        if ARGS.req_loci_file is not None:
            ARGS.required_loci = read_list_file(ARGS.req_loci_file)
    except IOError:
        LOGGER.error("Cannot open required loci file: %s",
                     ARGS.req_loci_file)
        raise

    try:
        if ARGS.excl_loci_file is not None:
            ARGS.exclude_loci = read_list_file(ARGS.excl_loci_file)
    except IOError:
        LOGGER.error("Cannot open excluded loci file: %s",
                     ARGS.excl_loci_file)
        raise

    try:
        if ARGS.excl_strains_file is not None:
            ARGS.exclude_strains = read_list_file(ARGS.excl_strains_file)
    except IOError:
        LOGGER.error("Cannot open exclude strains file: %s",
                     ARGS.excl_strains_file)
        raise

    PARAMS = {"res": ARGS.res,
              "stop_at_res": ARGS.stop_at_res,
              "n_treads": ARGS.threads,
              "max_loci": ARGS.max_loci,
              "max_res": ARGS.max_res,
              "reps": ARGS.reps,
              "exclude_loci": ARGS.exclude_loci,
              "required_loci": ARGS.required_loci,
              "exclude_strains": ARGS.exclude_strains
              }

    pattern_selection(PROJECT_DIR, **PARAMS)
