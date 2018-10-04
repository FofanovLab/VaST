import argparse
import os
import pandas as pd
from utils import path_type, file_type

def get_outpath(path, infile, ext="csv"):
    return os.path.join(path,
        "{0}.{1}".format(
        os.path.splitext(
            os.path.basename(infile))[0],
            ext))
    
def parse_index(index):
    return "{0}::{1}".format(
        index,
        index.split("::")[-1])

class Matrix:
    def __init__(self, path):
        self.path = path
        self.matrix = self.get_matrix()

    def write_matrix(self, path, sep=","):
        self.matrix.to_csv(path, sep=sep)
        

class Master(Matrix):
    def __init__(self, path):
        Matrix.__init__(self, path)

    def get_matrix(self):
        return parse_NASP(self.path)

class SNP(Matrix):
    def __init__(self, path):
        Matrix.__init__(self, path)

    def get_matrix(self):
        matrix = parse_NASP(self.path)
        matrix.index = matrix.index.map(parse_index)
        return matrix

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


if __name__ == "__main__":

    PARSER = argparse.ArgumentParser(
        prog='NASP Converter',
        description="Converts NASP matrix into VaST input files",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    PARSER.add_argument(
        "outpath", metavar="OUTPATH", type=path_type,
        help="Output directory path"
    )

    PARSER.add_argument(
        "--master", type=file_type, default=None,
        help="Path to master SNP matrix to convert"
    )

    PARSER.add_argument(
        "--snp", type=file_type, default=None,
        help="Path to SNP matrix to convert"
    )

    ARGS = PARSER.parse_args()

    if ARGS.master != None:
        master = Master(ARGS.master)
        master.write_matrix(
            get_outpath(
                ARGS.outpath, ARGS.master))

    if ARGS.snp != None:
        snp = SNP(ARGS.snp)
        snp.write_matrix(get_outpath(
            ARGS.outpath, ARGS.snp
        ))
