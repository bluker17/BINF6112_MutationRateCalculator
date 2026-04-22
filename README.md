# BINF6112_MutationRateCalculator
**UNCC BINF6112 - Programming II April 2nd Challenges**

## Authors:
**Denis Jacob Machado**
dmachado@charlotte.edu

**Bobby Luker**
rluker@charlotte.edu
UNCC ID: 801484356

## Program Description
Given a FASTA file containing multiple query sequences, a reference FASTA file, and parameters for match, mismatch, and indel costs, the mutation rate calculator will determine alignments between each query sequence and the reference sequence. It will then return all results in a CSV file.

## License: 
**GNU General Public License Version 3**

The GNU GPL is a license that ensures code is open-source. GNU GPL allows others to utilize, modify, or distribute code. If other users modify the code, then these users are expected to share their changes to the code under a GNU GPL to maintain the open-source integrity of the code.

[Project URL](https://github.com/bluker17/BINF6112_MutationRateCalculator)


## Project File Structure:

```
└── 📁src
    └── 📁consense
        ├── __init__.py
        ├── consense.py
    └── 📁mut_calc
        ├── __init__.py
        ├── mut_rate.py
    └── 📁newu2
        ├── __init__.py
        ├── num.py
        ├── read.py
        ├── write.py
    └── 📁score
        ├── __init__.py
        ├── align_score.py
└── 📁testing_materials
    └── 📁example_data
        ├── reference.fasta
        ├── sequences.fasta
    └── 📁example_outputs
    ├── run_test.sh
├── environment-alternative.yml
├── environment.yml
├── LICENSE
├── main.py
└── README.md
```
## Overview:
`src` contains the `consense`, `mut_calc`, `newu2`, and `score` packages.

* `consense`: Contains `consense.py` which determines the consensus sequence between the query and reference se quences.

* `mut_calc`: Contains `mut_rate.py` which calculates the mutation rate between the query and referecne sequences. It also returns all relevant information for the two sequences in a pandas data frame. 

* `newu2`: Contains multiple modules to handle input and output files.
    1. `num.py` completes a Needleman-Wunsch alignment between query and reference sequences. 
    2. `read.py` reads input FASTA files.
    3. `write.py` writes output CSV files using the generated pandas data frame from `mut_rate.py`

* `score`: Contains `align_score.py` which calculates the alignment score following Needleman-Wunsch alignment.

`main.py` executes all modules to generate results.

`testing_materials` contains `example_data` and `example_outputs` subdirectories, as well as `run_test.sh`.

* `run_test.sh`: Bash script that executes a test run for the user.

* `example_outputs`: Contains all generated CSV files from `run_test.sh`.

* `example_data` contains default files for the program.
    1. `sequences.fasta`: Default query sequences FASTA file.
    2. `reference.fasta`: Default reference FASTA file.


## Testing Instructions:
04/06/2026

### Installation

For MacOS:
```bash
conda env create -f environment.yml
```
For Windows/Linux:
```bash
conda env create -f environment-alternative.yml
```
Conda will automatically create an environment named mut_calculator with all the specified packages and versions.

### Usage

1. Activate the environment:
```bash
conda activate mut_calculator
```
2. Run following command to test:
```bash
bash testing/run_test.sh
```

#### Command-Line Arguments:
| Argument                | Description                                  | Default         |
| ----------------------- | -------------------------------------------- | --------------- |
| `--indel`| Cost for InDels | -1         |
| `--mismastch` | Cost for mismatches                           | -1|
| `--match`       | Cost for matches             | 1|
| `-r`, `--reference_file` | Reference FASTA filepath                    | reference.fasta|
| `-i`, `--input_file` | Query sequences FASTA filepath                    | sequences.fasta|
| `-0`, `--output_file` | Output CSV filepath                    | alignments.csv|

Expected Output:

- Prints alignment progress. 
- Upon alignment completion, the generated output CSV filepath is printed.

## Contributions
Dr. Jacob Machado provided extensive groundwork for modules found in `newu2` and execution of these modules in `main.py`.

Bobby created modules in `consense`, `mut_calc`, and `score` to generate results from Dr. Jacob Machado's alignments from `newu2`. Additionally, `main.py` was updated to handle these modules, as well as updated `num.py` to successfully align multiple query sequences to a reference sequence.

## References
OpenAI's ChatGPT model GPT-5.3 was used to guide coding decisions on creating a consensus sequence and how to handle numpy arrays for calculating mutation rate.