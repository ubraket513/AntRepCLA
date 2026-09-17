"""Construction of the Hamming graph over CDR3 sequences.

Vertices are distinct CDR3 nucleotide sequences. Two vertices are joined when

1. the CDR3s have the same length,
2. their Hamming distance is within ``hamming_threshold_ratio`` of that length,
3. they share at least one V gene and at least one J gene (alleles stripped).

Connected components of the resulting graph are the clonal lineages.
"""

from __future__ import annotations

from collections import defaultdict
from itertools import combinations

import networkx as nx
import pandas as pd
from loguru import logger


def hamming_distance(str1: str, str2: str) -> int:
    """Number of differing characters between two equal-length strings."""
    return sum(symb1 != symb2 for symb1, symb2 in zip(str1, str2))


def collapse_to_unique_cdr3(records: dict) -> dict:
    """Collapse reads to distinct CDR3s, unioning gene calls across reads.

    A repertoire carries heavy read-level redundancy -- roughly fourfold in the
    reference dataset -- and the graph is keyed on CDR3 sequence, so comparing
    the same pair of distinct sequences once per read pair is wasted work.

    The union matters: the same CDR3 is sometimes reported with different gene
    calls on different reads. Keeping only the first read's calls would hide
    routes through which the intersection test could have succeeded, silently
    dropping edges.
    """
    collapsed: dict = {}
    for rec in records.values():
        cdr3 = rec['cdr3']
        entry = collapsed.get(cdr3)
        if entry is None:
            collapsed[cdr3] = {
                'cdr3': cdr3,
                'cdr_length': rec['cdr_length'],
                'v_call': set(rec['v_call']),
                'j_call': set(rec['j_call']),
            }
        else:
            entry['v_call'].update(rec['v_call'])
            entry['j_call'].update(rec['j_call'])
    return collapsed


def build_blocks(records: dict) -> dict:
    """Index records by every ``(cdr_length, v_gene, j_gene)`` they carry.

    Conditions 1 and 3 of the edge predicate are equality constraints, so they
    double as blocking keys: a pair can only be joined if it shares such a
    triple. A record with several V or J calls is indexed under each
    combination, which is what keeps the scan exact rather than approximate.
    """
    blocks: dict = defaultdict(list)
    for key, rec in records.items():
        for v_gene in rec['v_call']:
            for j_gene in rec['j_call']:
                blocks[(rec['cdr_length'], v_gene, j_gene)].append(key)
    return blocks


def build_hamming_graph(
    igblast_result: pd.DataFrame,
    threshold_ratio: float = 0.1,
) -> nx.Graph:
    """Build the Hamming graph from normalised IgBLAST results.

    Only pairs sharing a blocking key are compared, which reduces the ~50
    million pairs of a 10,000-read repertoire to a few tens of thousands
    without changing the resulting graph.

    Args:
        igblast_result: Frame with 'cdr3', 'cdr_length', 'v_call', 'j_call'.
        threshold_ratio: Maximum Hamming distance as a fraction of CDR3 length.

    Returns:
        The Hamming graph, with distinct CDR3 sequences as nodes.
    """
    graph = nx.Graph()
    unique = collapse_to_unique_cdr3(igblast_result.to_dict('index'))

    # Isolated CDR3s are still lineages of size one, so every sequence is a node.
    graph.add_nodes_from(unique)

    blocks = build_blocks(unique)

    # A pair sharing several blocks would otherwise be compared repeatedly.
    candidates = set()
    for members in blocks.values():
        if len(members) < 2:
            continue
        for key_a, key_b in combinations(members, 2):
            candidates.add((key_a, key_b) if key_a < key_b else (key_b, key_a))

    for key_a, key_b in candidates:
        # The blocking key already guarantees equal length and gene overlap.
        if hamming_distance(key_a, key_b) <= threshold_ratio * unique[key_a]['cdr_length']:
            graph.add_edge(key_a, key_b)

    logger.info(
        f"Hamming graph: {graph.number_of_nodes():,} nodes, "
        f"{graph.number_of_edges():,} edges, "
        f"{len(candidates):,} candidate pairs compared."
    )
    return graph


def build_hamming_graph_reference(
    igblast_result: pd.DataFrame,
    threshold_ratio: float = 0.1,
) -> nx.Graph:
    """Build the graph by comparing every pair of reads.

    This is the original quadratic formulation, kept as the correctness oracle
    for :func:`build_hamming_graph` and used only by the test suite.

    It emits a self-loop for each pair of distinct reads carrying an identical
    CDR3, since it iterates over reads while keying nodes by sequence. Those
    self-loops are artifacts and do not affect connected components; the
    blocking implementation does not produce them.
    """
    graph = nx.Graph()
    records = igblast_result.to_dict('index')

    for key_a, key_b in combinations(records.keys(), 2):
        seq_a = records[key_a]
        seq_b = records[key_b]

        graph.add_node(seq_a['cdr3'])
        graph.add_node(seq_b['cdr3'])

        if seq_a['cdr_length'] != seq_b['cdr_length']:
            continue
        if hamming_distance(seq_a['cdr3'], seq_b['cdr3']) > threshold_ratio * seq_a['cdr_length']:
            continue
        if not seq_a['v_call'].intersection(seq_b['v_call']):
            continue
        if not seq_a['j_call'].intersection(seq_b['j_call']):
            continue

        graph.add_edge(seq_a['cdr3'], seq_b['cdr3'])

    return graph
