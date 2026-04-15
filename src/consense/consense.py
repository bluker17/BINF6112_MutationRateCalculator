#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

IUPAC = {
    'A': {'A'},
    'C': {'C'},
    'G': {'G'},
    'T': {'T'},
    'R': {'A', 'G'},
    'Y': {'C', 'T'},
    'S': {'G', 'C'},
    'W': {'A', 'T'},
    'K': {'G', 'T'},
    'M': {'A', 'C'},
    'B': {'C', 'G', 'T'},
    'D': {'A', 'G', 'T'},
    'H': {'A', 'C', 'T'},
    'V': {'A', 'C', 'G'},
    'N': {'A', 'C', 'G', 'T'},
}

frozen_IUPAC = {frozenset(v): k for k, v in IUPAC.items()}

def get_consensus(seq1: np.ndarray, seq2: np.ndarray) -> np.ndarray:
	"""
	Generate a consensus sequence from two aligned sequences.

	Parameters
	----------
	seq1 : np.ndarray
		First aligned sequence (array of characters).
	seq2 : np.ndarray
		Second aligned sequence (array of characters).

	Returns
	-------
	np.ndarray
		An array where each position is a set containing:
		  - One element if both positions match.
		  - Two elements if the positions differ.
	"""
	if len(seq1) != len(seq2):
		raise ValueError("Aligned sequences must have the same length.")

	consensus = []
	for a, b in zip(seq1, seq2):
		if a == b:
			consensus.append({a})
		else:
			consensus.append({a, b})

	return np.array(consensus, dtype=object)

def create_consensus(consensus_pair: np.ndarray) -> str:
	"""
	Creates a consensus sequence from two aligned sequences.

	Parameters
	---------
	seq1 : np.ndarray
		First aligned sequence (array of characters).
	seq2 : np.ndarray
		Second aligned sequence (array of characters).

	Returns
	-------
	str
		The consensus sequence as a string.

	"""
	consensus_sequence = ''
	for base_pair in consensus_pair:
		if '-' in base_pair:
			cleaned = base_pair - {'-'}
		else:
			cleaned = base_pair

		if len(cleaned) == 0:
			consensus_sequence += '-'

		consensus_sequence += frozen_IUPAC.get(frozenset(cleaned))

	return consensus_sequence