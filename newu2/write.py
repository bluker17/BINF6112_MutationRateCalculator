#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Utility to write sequences and their alignment information into a CSV-formatted file.
"""
import csv
import pandas as pd


def write_csv(filepath: str) -> None:
    """
    Create and write the headers into a file in CSV format.

    Parameters
    ----------
    filepath : str
        Path to the output CSV file.

    Returns
    -------
    None
    """
    
    headers = ["Sequence_ID", "Sequence", "Reference_ID", "Reference_Sequence", "Consensuse_Sequence", "Alignment_Cost", "Mutation_Rate", "Mutation_Counts", "Alignment_Length"]

    with open(filepath, "w", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(headers)

def append_csv(data_frame: pd.DataFrame, filepath: str) -> None:
    """
    Append data from a pandas DataFrame containing alignment information to a CSV file.

    Parameters
    ----------
    pd.DataFrame
        DataFrame containing data in columns: Sequence_ID, Sequence, Reference_ID, Reference_Sequence, Alignment_Cost, Consensus_Sequence, Mutation_Rate, Mutations, Alignment_Length.
    filepath : str
        Path to the output CSV file.

    Returns
    -------
    None
    """
    
    data_frame.to_csv(filepath, mode="a", index=False, header=False)