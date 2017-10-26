# VaST: **Va**riant Site **S**train **T**yper
## Description
VaST finds the minimum number of variant site targets that are required to differentiate a collection of strains. VaST is specifically designed to assist with generating panels for amplicon sequencing (AmpSeq).

## Clone VaST
```
git clone https://TaraFurstenau@bitbucket.org/fofanovlab/vast.git
cd vast
```

## Create Conda Environment
```
conda env create -f vast_env.yml
source activate vast_env
```

## Preprocessing Module
The preprocessing module includes the Amplicon Filter and the Pattern Discovery steps and outputs a JSON file which can be passed to the Pattern Selection Module. The Amplicon Filter treates each variant site as a potential amplicon, combining adjacent sites as necessary and filters out any amplicons that may be difficult to amplify in all strains. During Pattern Discovery, amplicons are divided into groups based on their resolution patterns which describes how the strains vary at the amplicon.

### Inputs
VaST must be provided a variant site matrix (`VAR_MATRIX_PATH`) where each row represents a genomic site that varies across the columns of strains; the values in the marix characterize the state of each strain at the variable sites. Many different types of genomic variation can be included in this matrix (SNPs, indels, VNTRs) provided that the variable region is short enough to be captured in an AmpSeq reaction. The first column of the variant site matrix should contain a genome identifier, a start position, and an end position separated by two colons, (e.g genome123::115::115). The start and end position should be the same for SNPs and for VNTRs, the stopping position should be based on the longest repeat.

To run the Amplicon Filter Module, VaST requires information about the regions upstream and downstream of each of the variable sites. Therefore, either a full genome matrix (`FULL_MATRIX_PATH`) must be provided which should include a call for each position in the genome for all of the strains or a previously generated flag file (`FLAG_FILE_PATH`) must be provided. The full genome matrix can be generated through the alignment of genome assemblies to a reference genome (same reference that was used to identify variant sites) or from VCF files that contain calls for each position in the genome. The first column of the full genome matrix should have a genome ID and the position separated by two colons, (e.g genome123::115). See `./example` for examples.


### Usage
```
$ python ./VaST/preprocessing.py --help
usage: AMPLICON_FILTER [-h]
                       (--full_matrix_path FULL_MATRIX_PATH | --flag_file_path FLAG_FILE_PATH)
                       [-d {tab,comma,space}] [-w WINDOW] [-s]
                       [-c STRAIN_CUTOFF] [-z PZ_SIZE] [-l PZ_FILTER_LENGTH]
                       [-p PZ_FILTER_PERCENT] [--log {DEBUG,INFO,WARNING}]
                       [-x EXCLUDE_STRAINS] [-t THREADS]
                       PROJECT_NAME PROJECT_DIR VAR_MATRIX_PATH

Parses a variant site matrix to find candidate amplicons

positional arguments:
  PROJECT_NAME          Project name and output file prefix
  PROJECT_DIR           Project directory
  VAR_MATRIX_PATH       Path to variant site matrices

optional arguments:
  -h, --help            show this help message and exit
  -d {tab,comma,space}, --sep {tab,comma,space}
                        Specify the file delimiter of variant site matrix
                        (default: tab)
  -w WINDOW, --window WINDOW
                        Size of window for variable region. (HINT: Should be
                        small enough to fit in desired amplicon with room for
                        upstream and downstream primer zone) (default: 50)
  -s, --strict          Strict mode ignores sites with missing or ambiguous
                        data (default: False)
  -c STRAIN_CUTOFF, --strain_cutoff STRAIN_CUTOFF
                        The number of strains at a primer zone site that can
                        have a non-conserved call before the site is flagged
                        (default: 1)
  -z PZ_SIZE, --pz_size PZ_SIZE
                        Size of the primer zone (default: 200)
  -l PZ_FILTER_LENGTH, --pz_filter_length PZ_FILTER_LENGTH
                        The length of un-flagged primer zone segments that
                        count toward the primer zone filter percent (default:
                        25)
  -p PZ_FILTER_PERCENT, --pz_filter_percent PZ_FILTER_PERCENT
                        The percent of primer zone positions that must be
                        present in un-flagged segments of the primer zone that
                        are longer than the primer zone filter length
                        (default: 25)
  --log {DEBUG,INFO,WARNING}
                        Set logging level (default: INFO)
  -x EXCLUDE_STRAINS, --exclude_strains EXCLUDE_STRAINS
                        Path to file containing a list of strains to exclude
                        from analysis in a single column (default: None)
  -t THREADS, --threads THREADS
                        Number of threads for multiprocessing (default: 1)

Filter file type:
  Provide full genome matrix or an already calculated flag file

  --full_matrix_path FULL_MATRIX_PATH
                        Path to full genome matrix (includes calls at every
                        site in the genome) (default: None)
  --flag_file_path FLAG_FILE_PATH
                        Path to precomputed flag file (default: None)
```

### Parameters
| Parameter  | Description  | Notes  |
|---|---|---|
| Strict mode  | VaST ignores missing or ambiguous data in input matrix  |  Speeds up preprocessing ut some sites are lost |
|  Window Size | Maximum distance between adjacent sites that can be combined into a single amplicon  |The desired amplicon length should be considered when setting the window size. A larger window may increase the number of variant sites that are included in the amplicons making them more efficient   |
| Primer Zone Size  | Size of the region upstream and downstream of the target to evaluate in the amplicon filter  | The primer zones begin immediately before the first and immediately after the last target site in the window, so the maximum amplicon size is 2 x primer zone size + window size. A smaller primer zone may limit the number of primer options.  |
|Strain Cutoff |The number of strains at a primer zone site that can have a non-conserved call before the site is flagged | A strain cutoff greater than one will not guarantee that the primer zone sequences are conserved across all of the strains but it may be appropriate in cases where one or a few strains have low sequence coverage|
|Primer Zone Filter Percent|The percent of primer zone positions that must be present in un-flagged segments of the primer zone that are longer than the primer zone filter length| A higher primer zone filter percent will increase the total number of primer options in amplicons that pass the filter|
|Primer Zone Filter Length | The length of un-flagged primer zone segments that count toward the primer zone filter percent |The primer zone filter length should be at least as long as the minimum acceptable primer length to ensure that conserved primers can be found within the primer zone.|

## Pattern Selection Module
The Pattern Selection Module uses the patterns from Pattern Discovery and finds the minimum number of patterns that maximize strain resolution.

### Inputs
The Pattern Selection Module requires a project name and the path to the root of the Preprocessing project directory that contains the desired pattern JSON file. The output of the Pattern Selection Module will be placed in a directory within the Preprocessing project directory.

### Usage

```
$ python ./VaST/pattern_selection.py --help
usage: PATTERN_SELECTION [-h] [--res RES] [--stop_at_res]
                         [--log {DEBUG,INFO,WARNING}] [-t THREADS]
                         [--reps REPS]
                         [--max_loci MAX_LOCI | --max_res MAX_RES]
                         [--required_loci [REQUIRED_LOCI [REQUIRED_LOCI ...]]
                         | --req_loci_file REQ_LOCI_FILE]
                         [--exclude_loci [EXCLUDE_LOCI [EXCLUDE_LOCI ...]] |
                         --excl_loci_file EXCL_LOCI_FILE]
                         [--exclude_strains [EXCLUDE_STRAINS [EXCLUDE_STRAINS ...]]
                         | --excl_strains_file EXCL_STRAINS_FILE]
                         PROJECT_NAME PROJECT_DIR

Find the minimum set of patterns that uniquely identify all strains

positional arguments:
  PROJECT_NAME          Project name and output file prefix
  PROJECT_DIR           Path to preprocessing directory

optional arguments:
  -h, --help            show this help message and exit
  --res RES             Path to file which defines the desired level of
                        resolution for each of the strains (the default if
                        full strain resolution) The resolution file should be
                        a csv table with each row starting with a strain name
                        (exactly as it appears in the variant site matrix)
                        followed by a group ID. Strains with the same group ID
                        will not be resolved. Multiple columns of group IDs
                        may be present but must start with the lowest
                        resolution grouping and increase to higher resolution
                        groupings and must be consistent. (default: None)
  --stop_at_res         Force solution to stop once the maximum level
                        specified in resolution file has been reached.
                        Otherwise, the solution will continue towards full
                        strain resolution once it has achieved the level of
                        resolution defined in the file. It will only terminate
                        early when MAX_LOCI is reached or no more resolution
                        is achievable. (default: False)
  --log {DEBUG,INFO,WARNING}
                        Set logging level (default: INFO)
  -t THREADS, --threads THREADS
                        Set number of threads (default: 1)
  --reps REPS           Run N sets at a time and pick the best (default: 1)

Early Termination:
  Options for creating early stopping points

  --max_loci MAX_LOCI   Stop after MAX_LOCI is reached. If required loci are
                        added, MAX_LOCI must be greater than or equal to the
                        number of required loci. (default: inf)
  --max_res MAX_RES     Stop when the resolution has reached MAX_RES percent
                        of total resolution. (default: 100)

Add required loci to solution:
  Add loci that must be included in the solution either from command line or
  from file. Names should match the locus ID given in the variant site
  matrix.

  --required_loci [REQUIRED_LOCI [REQUIRED_LOCI ...]]
                        List of locus IDs that must be included in solution.
                        (default: [])
  --req_loci_file REQ_LOCI_FILE
                        Path to file containing list of locus IDs (in a single
                        column) that must be included in solution. (default:
                        None)

Exclude loci from consideration:
  Add loci either from command line or from file. These loci will removed
  from consideration when building the solution. Names should match the
  locus ID given in the variant site matrix.

  --exclude_loci [EXCLUDE_LOCI [EXCLUDE_LOCI ...]]
                        List of loci that must be excluded from solution
                        (default: [])
  --excl_loci_file EXCL_LOCI_FILE
                        Path to file containing list of locus IDs (in a single
                        column) that must be excluded from solution. (default:
                        None)

Exclude strains from consideration:
  Add strains either from command line or from file. Names should match the
  variant site matrix.

  --exclude_strains [EXCLUDE_STRAINS [EXCLUDE_STRAINS ...]]
                        List of strains to be excluded (default: [])
  --excl_strains_file EXCL_STRAINS_FILE
                        Path to file containing list of strains (in a single
                        column) that should be excluded (default: None)
```
