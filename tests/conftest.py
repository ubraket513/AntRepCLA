"""Shared fixtures and helpers."""

from pathlib import Path

import networkx as nx
import pandas as pd
import pytest

REAL_DATASET = Path(__file__).resolve().parent.parent / "data" / "igblast_results.tsv"

requires_real_dataset = pytest.mark.skipif(
    not REAL_DATASET.is_file(),
    reason="real IgBLAST dataset not present",
)


def make_df(rows):
    """Build an igblast-shaped frame from (cdr3, v_call, j_call) triples."""
    return pd.DataFrame([
        {
            'cdr3': cdr3,
            'cdr3_aa': 'X' * (len(cdr3) // 3),
            'cdr_length': len(cdr3),
            'v_call': set(v),
            'j_call': set(j),
        }
        for cdr3, v, j in rows
    ])


def genuine_edges(graph):
    """Edge set with self-loops removed."""
    return {frozenset(e) for e in graph.edges() if e[0] != e[1]}


def components(graph):
    return sorted(sorted(c) for c in nx.connected_components(graph))
