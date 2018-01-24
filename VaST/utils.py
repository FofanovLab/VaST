import argparse
import datetime
import glob
import logging
import os.path
from functools import partial
from itertools import combinations
from multiprocessing import Pool

import numpy as np
import pandas as pd

BASES = ['A', 'C', 'T', 'G']
AMBIGUOUS_DICT = {
    'A': ['A'],
    'B': ['C', 'G', 'T'],
    'C': ['C'],
    'D': ['A', 'G', 'T'],
    'G': ['G'],
    'H': ['A', 'C', 'T'],
    'K': ['G', 'T'],
    'M': ['A', 'C'],
    'N': ['G', 'A', 'T', 'C'],
    'R': ['A', 'G'],
    'S': ['C', 'G'],
    'T': ['T'],
    'V': ['A', 'C', 'G'],
    'W': ['A', 'T'],
    'X': ['G', 'A', 'T', 'C', '-'],
    'Y': ['C', 'T']}
AMBIGUOUS = np.setdiff1d(AMBIGUOUS_DICT.keys(), BASES)


class Project_Directory:
    def __init__(self, path, project_name, sub_dirs, parent_dirs=[]):
        self.path = path
        self.timestamp = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
        self.project_name = "{0}_{1}".format(project_name, self.timestamp)
        self.project_path = os.path.join(
            os.path.abspath(path), "{0}/".format(self.project_name))
        self.sub_dirs = {
            sub_dir: self.prepend_path(sub_dir) for sub_dir in sub_dirs}
        self._make_directory_tree()
        self.parent_dirs = {
            parent_dir: os.path.abspath(os.path.join(self.project_path, '..', parent_dir))
            for parent_dir in parent_dirs}

    def get_sub_directory(self, sub_directory_name):
        return self.sub_dirs[sub_directory_name]

    def get_parent_subdirectory(self, parent_directory_name):
        return self.parent_dirs[parent_directory_name]

    def get_parent_subdirectory_file(
            self, parent_directory_name, file_name):
        return glob.glob(os.path.join(
            self.get_parent_subdirectory(
                parent_directory_name),
            file_name))[0]

    def get_parent_directory_timestamp(self):
        return "_".join(self.get_parent_directory().split('/')[-1].split("_")[-2:])

    def get_parent_directory(self):
        return os.path.abspath(os.path.join(self.project_path, '../'))

    def prepend_path(self, sub_directory_name):
        return os.path.join(self.project_path, "{}/".format(sub_directory_name))

    def make_subdirectory(self, sub_directory):
        os.mkdir(sub_directory)

    def new_sub_directory(self, sub_directory_name):
        self.sub_dirs[sub_directory_name] = self.prepend_path(
            sub_directory_name)

    def make_new_file(self, sub_directory_name, file_name, extension="txt"):
        return os.path.join(
            self.get_sub_directory(sub_directory_name),
            "{0}_{1}.{2}".format(file_name, self.timestamp, extension))

    def _make_directory_tree(self):
        """ Creates a directory tree. The root of the tree is the Project_Directory
        instance's project name with a unique timestamp
        """
        try:
            os.mkdir(self.project_path)
        except OSError:
            raise OSError("Path does not exist: {}".format(self.path))
        for sub_dir in self.sub_dirs.values():
            self.make_subdirectory(sub_dir)

class History:
    def __init__(self, history_path, program_name,
                 timestamp=None, path_dict={}, param_dict={},
                 other_dict={}, exists=False):
        self._file = history_path
        self._program_name = program_name
        self._path_dict = path_dict
        self._param_dict = param_dict
        self._other_dict = other_dict
        self._timestamp = timestamp
        if exists:
            self._read_history()

    def add_path(self, path_name, path):
        self._path_dict[path_name] = path

    def add_parameter(self, parameter_name, parameter):
        self._param_dict[parameter_name] = parameter

    def add_other(self, other_name, other):
        self._other_dict[other_name] = other

    def get_parameter(self, parameter_name):
        return self._param_dict[parameter_name]

    def get_path(self, path_name):
        return self._path_dict[path_name]

    def get_other(self, other_name):
        return self._other_dict[other_name]

    def write(self):
        with open(self._file, 'w') as history_out:
            history_out.write(
                "PROGRAM: {0}\nTIMESTAMP: {1}\n".format(
                    self._program_name, self._timestamp))
            form = "{0}: {1}\n"
            if self._path_dict:
                self._write_dict(history_out, "PATHS", self._path_dict)
            if self._param_dict:
                self._write_dict(history_out, "PARAMETERS", self._param_dict)
            if self._other_dict:
                self._write_dict(history_out, "OTHER", self._other_dict)

    def _write_dict(self, handle, name, dic):
        handle.write("{}:\n".format(name))
        for name, value in dic.iteritems():
            handle.write(
                "{0}: {1}\n".format(name.upper(), value))

    def _read_history(self):
        with open(self._file, 'rU') as infile:
            dic_type = {
                "PATHS": self._path_dict,
                "PARAMETERS": self._param_dict,
                "OTHER": self._other_dict}
            current_dic = None
            infile.readline()
            infile.readline()
            for line in infile:
                key, value = line.strip().split(":")
                if not value and key in dic_type:
                    current_dic = dic_type[key]
                else:
                    current_dic[key] = value.strip()


# ARGPARSE types
def positive_int(input_val):
    ''' Make a positive int type for argparse'''
    try:
        input_val = int(input_val)
        if input_val <= 0:
            raise ValueError
    except ValueError:
        raise argparse.ArgumentTypeError("Not a positive integer")
    return input_val


def percent(input_val):
    ''' Make percent type for argparse'''
    try:
        input_val = float(input_val)
        if input_val < 0 or input_val > 100:
            raise ValueError
    except ValueError:
        raise argparse.ArgumentTypeError("Not a percent")
    return input_val


def path_type(input_path):
    if not os.path.isdir(input_path):
        raise argparse.ArgumentTypeError("Not a valid path")
    return os.path.abspath(input_path)


def file_type(input_file):
    if not os.path.isfile(input_file):
        raise argparse.ArgumentTypeError("Not a valid file path")
    return os.path.abspath(input_file)

#########


def parse_NASP(nasp_file_name, sep="\t", split="#SNPcall", **kwargs):
    """
    kwargs: Arguments passed to pandas.read_csv()

    Returns:
    Pandas dataframe with metadata after split removed. Each row is a locus
    and each column is a strain

    Raises:
    ValueError if split not found in header indicating an inappropriate file
    format
    """

    with open(nasp_file_name, 'r') as infile:
        header = infile.readline().split(sep)
    try:
        index = header.index(split)
    except ValueError:
        raise ValueError("Not a valid NASP file format")

    return pd.read_csv(
        nasp_file_name, sep=sep, usecols=range(0, index),
        index_col=0, **kwargs)


def config_logging(log_file_name, level):
    """Configure logging

    Arguments:
        log_file_name (str): Path to log file
        level (str): Logging level

    Returns:
        Logging object

    Raises:
        IOError if log file path is not valid
    """

    if os.path.isdir(os.path.dirname(log_file_name)):
        logging.basicConfig(
            filename=log_file_name,
            datefmt='%m/%d/%Y %I:%M:%S %p',
            level=getattr(logging, level),
            filemode='w',
            format='%(asctime)s %(levelname)s: [%(name)s] %(message)s')
    else:
        raise IOError("Invalid log file path")


def read_list_file(file_name, sep="\n", dtype=str):
    with open(file_name, 'rU') as handle:
        return list(np.atleast_1d(np.genfromtxt(
            handle, delimiter=sep, dtype=dtype)))


def parallel_apply(df, func, threads, axis, **kwargs):
    num_partitions = df.shape[axis] // threads
    df_split = np.array_split(df, num_partitions)
    pool = Pool(threads)
    func_partial = partial(func, **kwargs)
    apply_function = partial(
        parallel_apply_helper, func=func_partial, axis=axis)
    df = pd.concat(pool.map(apply_function, df_split))
    pool.close()
    pool.join()
    return df


def normalize(value, _max, _min=0):
    return (value - _min) / float(_max - _min)


def parallel_apply_helper(df, func, axis):
    return df.apply(func, axis=axis)


def get_locus_groups(row, alphabet):
    size = len(row)
    pattern = [[] for i in xrange(size)]
    if np.all(np.in1d(row, ["A", "T", "C", "G"])):
        uni = np.unique(row, return_inverse=True)[1]
        for i in xrange(max(uni) + 1):
            pos = np.where(uni == i)[0]
            for p in pos:
                pattern[p] = pos[pos != p]
    else:
        row = [set(alphabet[base]) for base in row]
        for i, j in combinations(range(size), 2):
            if row[i].intersection(row[j]):
                pattern[i].append(j)
                pattern[j].append(i)
    return [tuple(p) for p in pattern]


def cartesian_product(arr, out=None, dtype="S5"):    
    """
    Generate a cartesian product of input arrays.
        Parameters
        From https://gist.github.com/hernamesbarbara/68d073f551565de02ac5
        ----------
        arrays : list of array-like
            1-D arrays to form the cartesian product of.
        Returns
        -------
        out : ndarray
            2-D array of shape (M, len(arrays)) containing cartesian products
            formed of input arrays.
        Examples
        --------
        >>> cartesian(([1, 2, 3], [4, 5], [6, 7]))
        array([[1, 4, 6],
                [1, 4, 7],
                [1, 5, 6],
                [1, 5, 7],
                [2, 4, 6],
                [2, 4, 7],
                [2, 5, 6],
                [2, 5, 7],
                [3, 4, 6],
                [3, 4, 7],
                [3, 5, 6],
                [3, 5, 7]])
        """
    arrays = [np.asarray(x) for x in arr]

    n = np.prod([x.size for x in arrays])
    if out is None:
        out = np.zeros([n, len(arrays)], dtype=dtype)

    m = n / arrays[0].size
    out[:, 0] = np.repeat(arrays[0], m)
    if arrays[1:]:
        cartesian_product(
            arrays[1:], out=out[0:m, 1:], dtype=dtype)
        for j in xrange(1, arrays[0].size):
            out[j * m:(j + 1) * m, 1:] = out[0:m, 1:]
    return out


def get_ambiguous_pattern(feature, feature_categories):
    pattern = []
    for i, fc in enumerate(feature_categories):
        for j in feature:
            if j in fc:
                pattern.append(i)
                break
    if not len(pattern):
        feature_categories.append(feature)
        return (len(feature_categories) - 1, )
    return tuple(set(pattern))

def parse_flag_file(flag_file_path, sep=","):
    ''' Parses csv flag file with headers: 
    Site, Genome, Flag. '''
    flag_df = pd.read_csv(flag_file_path, sep=sep)
    # Check headers are correct
    headers = [h.lower() for h in flag_df.columns.values]
    if not np.all(np.in1d(["site", "genome", "flag"], headers)):
        raise IOError("Flag file is missing expected headers")
    return flag_df

