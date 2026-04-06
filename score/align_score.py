#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

def get_bottom_right_value(matrix: np.ndarray):
	"""
	Return the value of the bottom-right cell in a NumPy matrix.

	Parameters
	----------
	matrix : np.ndarray
		A 2D NumPy array.

	Returns
	-------
	The value stored in the bottom-right corner of the matrix.
	"""
	if matrix.ndim != 2:
		raise ValueError("Input must be a 2D NumPy array.")
	return matrix[-1, -1]