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

[Project URL](https://github.com/bluker17/BINF6112_MutationRateCalculator)

## License: 
**GNU General Public License Version 3**

Review `LICENSE` for more information.

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

## References:

### Python Standard Library

**`argparse`**  
Python Software Foundation. (2024). *argparse — Parser for command-line options, arguments and sub-commands*. Python 3 Documentation.  
https://docs.python.org/3/library/argparse.html  

Used for parsing command-line arguments and handling CLI input configuration.

---

**`csv`**
csv Development Team. (2024). *csv — CSV File Reading and Writing*. Python 3 Documentation.
https://docs.python.org/3/library/csv.html

Used for reading and writing comma-separated value (CSV) files, providing tools to efficiently parse and generate tabular data in text format.

---

**`pandas`**
pandas Development Team. (May 11, 2026). *pandas documentation*. Python Data Analysis Library.
https://pandas.pydata.org/docs/

Used for data manipulation and analysis in Python, providing high-performance DataFrame and Series structures for working with tabular and time-series data.

---

**`pathlib`**  
Python Software Foundation. (2024). *pathlib — Object-oriented filesystem paths*. Python 3 Documentation.  
https://docs.python.org/3/library/pathlib.html  

Used for platform-independent file and directory path handling.

---

**`re`**
re Development Team. (2024). *re — Regular expression operations*. Python 3 Documentation.
https://docs.python.org/3/library/re.html

Used for working with regular expressions in Python, enabling pattern matching, searching, splitting, and text manipulation on strings.

---

**`sys`**  
Python Software Foundation. (2024). *sys — System-specific parameters and functions*. Python 3 Documentation.  
https://docs.python.org/3/library/sys.html  

Used for interacting with interpreter-level functionality such as command-line arguments and program exit handling.

### AI assistance:

This project was developed with the help of [ChatGPT-5.5](https://chatgpt.com) by [OpenAI](https://openai.com).

ChatGPT assisted with:
- Code architecture and implementation
- Debugging and code review

All generated code was reviewed and tested by the author.

