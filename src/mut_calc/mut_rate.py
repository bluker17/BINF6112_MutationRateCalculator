#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd

def mutation_rate(aln: np.ndarray, seq_id: str, seq: str, ref_id: str, ref_seq: str, cost: int, con: str) -> pd.DataFrame:
    """
    Calculate the mutation rate between a query sequence and a reference sequence based on their alignment.

    Parameters
    ----------
    aln: np.ndarray
        A 2D NumPy array containing the aligned sequences (query and reference).
    seq_id: str
        Identifier for the query sequence.
    seq: str
        The query sequence.
    ref_id: str
        Identifier for the reference sequence.
    ref_seq: str
        The reference sequence.
    cost: int
        The alignment cost (score) for the aligned sequences.
    con: str
        The consensus sequence derived from the alignment.

    Returns
    -------
    pd.DataFrame
        A DataFrame containing the mutation rate and all other related information for the aligned sequences.
    """
    mutations = 0
    
    for i in range(len(aln[0])):
        if aln[0][i] != aln[1][i]:
            mutations += 1
    
    rate = mutations / len(aln[0])
    
    df = pd.DataFrame([{
        "sequence_id": seq_id,
        "sequence": seq,
        "reference_id": ref_id,
        "reference_sequence": ref_seq,
        "consensus_sequence": con,
        "alignment_cost": cost,
        "mutation_rate": rate,
        "mutations": mutations,
        "alignment_length": len(aln[0])
    }])
    
    return df