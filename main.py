#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
main.py

Driver script for running the Needleman-Wunsch algorithm demonstration.
Uses the newu2 module for alignment functions.
Uses score module to collect alignment scores.
Uses consense module to generate consensus sequences.
Uses mut_calc module to calculate mutation rates.
"""

import sys, argparse
from pathlib import Path

from src.newu2.read import read_fa
from src.newu2.write import append_csv, write_csv
from src.newu2.num import parse_fa, fill_matrix, trace_matrix
from src.consense.consense import get_consensus
from src.consense.consense import create_consensus
from src.score.align_score import get_bottom_right_value
from src.mut_calc.mut_rate import mutation_rate


class AlignmentRunner:
	"""
	AlignmentRunner is responsible for managing sequence alignment tasks.
	It loads input sequences, applies the Needleman-Wunsch algorithm, 
	and outputs results.
	"""

	def __init__(self, fasta_file: str, aln_file: str, reference_file: str, match:int=1, mismatch:int=-1, indel:int=-1) -> None:
		"""
		Initialize the runner with a fasta file.

		Parameters
		----------
		fasta_file : str
			Path to the input fasta file containing sequences to align.
		"""
		self.fasta_file = fasta_file
		sequences = read_fa(fasta_file)
		self.sequences = sequences

		reference  = read_fa(reference_file)
		self.ref_id = list(reference.keys())[0]
		self.ref_seq = list(reference.values())[0]


		self.match = match
		self.mismatch = mismatch
		self.indel = indel
		self.aln_file = aln_file

		self.mat = None
		self.seq1 = None
		self.id1 = None

	def fill(self) -> None:
		"""
		Filling in the table.
		"""
		self.mat = fill_matrix(
			matrix = self.mat,
			seq1_array = self.seq1,
			seq2_array = self.ref_seq,
			match = self.match,
			mismatch = self.mismatch,
			indel = self.indel
		)

	def trace(self) -> None:
		"""
		Tracing back to produce the alignment.
		"""
		self.aln = trace_matrix(
			matrix = self.mat,
			seq1_array = self.seq1,
			seq2_array = self.ref_seq,
		)

	def consense_pair(self) -> None:
		seq1, ref_seq = ["".join(seq) for seq in self.aln]
		self.consensus = get_consensus(seq1, ref_seq)

	def report_mat_cost(self) -> None:
		self.alncost = get_bottom_right_value(self.mat)

def parse_args() -> argparse.Namespace:
	"""
	Parse command-line arguments for the program.

	Returns
	-------
	argparse.Namespace
		Object containing parsed arguments.
	"""
	parser = argparse.ArgumentParser(
		description="Needleman-Wunsch global sequence alignment demo"
	)
	parser.add_argument(
		"-i", "--input_file",
		required=False,
		default="testing_materials/example_data/sequences.fasta",
		type=str,
		help="Path to input fasta file containing sequences"
	)
	parser.add_argument(
		"-o", "--output_file",
		required=False,
		default="output/alignments.csv",
		type=str,
		help="Path to output fasta file containing the sequence alignment"
	)
	parser.add_argument(
		"--match",
		required=False,
		default=1,
		type=int,
		help="Cost for matches"
	)
	parser.add_argument(
		"--mismatch",
		required=False,
		default=-1,
		type=int,
		help="Cost for mismatches"
	)
	parser.add_argument(
		"--indel",
		required=False,
		default=-1,
		type=int,
		help="Cost for InDels"
	)
	parser.add_argument(
		"-r", "--reference_file",
		required=False,
		default="testing_materials/example_data/reference.fasta",
		type=str,
		help="Path to reference sequence file"
	)

	return parser.parse_args()

def validate_args(args: argparse.Namespace) -> None:
	"""
	Validate the parsed command-line arguments.

	Parameters
	----------
	args : argparse.Namespace
		Object containing parsed arguments.

	Raises
	------
	ValueError
		If any of the arguments are invalid.
	"""
	files= (".fasta", ".fa", ".fna", ".fas")
	if not args.input_file.endswith(tuple(files)):
		raise ValueError("Input file must be a FASTA file with .fasta, .fa, .fna, or .fas extension.")
	
	if not Path(args.input_file).exists():
		raise ValueError(f"Input file '{args.input_file}' does not exist.")
	if not Path(args.reference_file).exists():
		raise ValueError(f"Reference file '{args.reference_file}' does not exist.")
	
	if not args.output_file.endswith(".csv"):
		raise ValueError("Output file must be a CSV file with .csv extension.")

	dir_path = Path(args.output_file).parent
	if not dir_path.is_dir():
		raise FileNotFoundError(f"Directory for output file does not exist: {dir_path}")
	
	if not args.reference_file.endswith(tuple(files)):
		raise ValueError("Reference file must be a FASTA file with .fasta, .fa, .fna, or .fas extension.")
	if read_fa(args.reference_file) is None or len(read_fa(args.reference_file)) > 1:
		raise ValueError("Reference file must have exactly one sequence")
	

def main() -> None:
	"""
	Main function to execute the alignment workflow.

	Returns
	-------
	None
	"""
	args = parse_args()
	validate_args(args)

	runner = AlignmentRunner(
		fasta_file=args.input_file,
		aln_file=args.output_file,
		match=args.match,
		mismatch=args.mismatch,
		indel=args.indel,
		reference_file = args.reference_file
	)

	write_csv(args.output_file)
	for id1, seq in runner.sequences.items():
		runner.mat, runner.seq1, runner.ref_seq = parse_fa(seq, runner.ref_seq)

		check_input_message = f"""
Aligning sequences:
	>{id1}
	{''.join(seq)}
	>{runner.ref_id}
	{''.join(runner.ref_seq)}
	"""
		sys.stdout.write(check_input_message)
		runner.fill()
		runner.trace()

		runner.consense_pair()
		consensus_sequence = create_consensus(runner.consensus)

		runner.report_mat_cost()

		mut_rate = mutation_rate(runner.aln, id1, seq, runner.ref_id, ''.join(runner.ref_seq), runner.alncost, consensus_sequence)
		append_csv(mut_rate, args.output_file)


	sys.stdout.write(f"\nAll alignments written into {args.output_file}\n\n")


if __name__ == "__main__":
	main()