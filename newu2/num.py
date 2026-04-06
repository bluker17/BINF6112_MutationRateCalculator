#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
num.py

Matrix and sequence utilities for the newu2 module.
This script defines the function parse_fa to transform
a formatted FASTA string into numpy arrays and an alignment matrix.
"""

import numpy as np

def parse_fa(seq1: str, seq2: str) -> np.ndarray:
    """
    Initialize an alignment matrix for two sequences.

    Parameters
    ----------
    seq1 : str
        The first sequence (query).
    seq2 : str
        The second sequence (reference).

    Returns
    -------
    tuple[np.ndarray, np.ndarray, np.ndarray]
		- Empty matrix of shape initialized to zeros, with dimensions (len(seq2)+1, len(seq1)+1)
			- Columns correspond to sequence 1; rows correspond to sequence 2.
		- First sequence as a NumPy array of characters
		- Second sequence as a NumPy array of characters

        A zero-initialized matrix of shape (len(seq2)+1, len(seq1)+1)
    """
    seq1_array = np.array(list(seq1))
    seq2_array = np.array(list(seq2))
    rows = len(seq2_array) + 1
    cols = len(seq1_array) + 1
    matrix = np.zeros((rows, cols), dtype=int)
    return matrix, seq1_array, seq2_array

def fill_matrix(
	matrix: np.ndarray,
	seq1_array: np.ndarray,
	seq2_array: np.ndarray,
	match: int = 1,
	mismatch: int = -1,
	indel: int = -1
) -> np.ndarray:
	"""
	Fill the scoring matrix using the Needleman–Wunsch down-pass algorithm.

	Parameters
	----------
	matrix : np.ndarray
		A zero-initialized matrix of shape (len(seq2)+1, len(seq1)+1).
		Columns correspond to sequence 1; rows correspond to sequence 2.
	seq1_array : np.ndarray
		NumPy array of characters representing sequence 1.
	seq2_array : np.ndarray
		NumPy array of characters representing sequence 2.
	match : int, optional
		Score for a character match (default = +1).
	mismatch : int, optional
		Score for a character mismatch (default = -1).
	indel : int, optional
		Penalty for an insertion/deletion (default = -1).

	Returns
	-------
	np.ndarray
		The filled scoring matrix.

	Notes
	-----
	- The first cell (0,0) is initialized to 0.
	- The first row and column are filled using only indel penalties.
	- Each subsequent cell (i,j) is filled using:
		- diagonal = matrix[i-1, j-1] + (match or mismatch)
		- up       = matrix[i-1, j] + indel
		- left     = matrix[i, j-1] + indel
	  The final cell value is max(diagonal, up, left).
	"""
	
	rows, cols = matrix.shape
	
	# Initialize the first row and first column with indel penalties
	for i in range(1, rows):
		matrix[i, 0] = matrix[i - 1, 0] + indel
	for j in range(1, cols):
		matrix[0, j] = matrix[0, j - 1] + indel
		
	# Fill the matrix row by row
	for i in range(1, rows):
		for j in range(1, cols):
			# Compare corresponding characters
			if seq2_array[i - 1] == seq1_array[j - 1]:
				score_diagonal = matrix[i - 1, j - 1] + match
			else:
				score_diagonal = matrix[i - 1, j - 1] + mismatch
				
			# Indel scores
			score_up = matrix[i - 1, j] + indel
			score_left = matrix[i, j - 1] + indel
			
			# Take maximum
			max_value = max(score_diagonal, score_up, score_left)
			matrix[i, j] = max_value
			
	return matrix

def trace_matrix(
	matrix: np.ndarray,
	seq1_array: np.ndarray,
	seq2_array: np.ndarray
) -> np.ndarray:
	"""
	Trace back through the filled Needleman–Wunsch matrix to reconstruct
	the optimal alignment path.

	Parameters
	----------
	matrix : np.ndarray
		The filled Needleman–Wunsch scoring matrix.
	seq1_array : np.ndarray
		NumPy array of characters representing sequence 1.
	seq2_array : np.ndarray
		NumPy array of characters representing sequence 2.

	Returns
	-------
	np.ndarray
		A NumPy array with two rows:
			- Row 0: aligned sequence 1
			- Row 1: aligned sequence 2
		Both aligned sequences have the same length.

	Notes
	-----
	- Tracing starts from the bottom-right corner.
	- At each step, we move in the direction that matches the cell's value:
		* Diagonal (match/mismatch)
		* Left  (gap in sequence 2)
		* Up    (gap in sequence 1)
	- When multiple paths are equally optimal:
		1. Prefer the diagonal move.
		2. If diagonal is not optimal, prefer the left move.
		3. Otherwise move up.
	"""
	i, j = matrix.shape[0] - 1, matrix.shape[1] - 1
	aligned_seq1 = []
	aligned_seq2 = []
	
	while i > 0 or j > 0:
		# Handle boundaries first
		if i == 0:
			aligned_seq1.append(seq1_array[j - 1])
			aligned_seq2.append('-')
			j -= 1
			continue
		elif j == 0:
			aligned_seq1.append('-')
			aligned_seq2.append(seq2_array[i - 1])
			i -= 1
			continue
		
		current = matrix[i, j]
		diag = matrix[i - 1, j - 1]
		up = matrix[i - 1, j]
		left = matrix[i, j - 1]
		
		# Determine which move leads to current
		if current == diag + 1 or current == diag - 1:
			# Diagonal move (priority 1)
			aligned_seq1.append(seq1_array[j - 1])
			aligned_seq2.append(seq2_array[i - 1])
			i -= 1
			j -= 1
		elif current == left - 1:
			# Move left (gap in seq2, priority 2)
			aligned_seq1.append(seq1_array[j - 1])
			aligned_seq2.append('-')
			j -= 1
		else:
			# Move up (gap in seq1, priority 3)
			aligned_seq1.append('-')
			aligned_seq2.append(seq2_array[i - 1])
			i -= 1
			
	aligned_seq1.reverse()
	aligned_seq2.reverse()
	
	return np.array([aligned_seq1, aligned_seq2])