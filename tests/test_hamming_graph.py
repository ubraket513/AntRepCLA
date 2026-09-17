"""Correctness of the blocking-key Hamming graph.

The reference is the quadratic implementation, which compares every pair of
reads. The two must agree on nodes, on genuine edges, and on connected
components -- the last being what the lineage analysis actually consumes.

They differ in exactly one respect: the reference emits a self-loop whenever
two distinct reads carry an identical CDR3, because it iterates over reads
while keying nodes by CDR3 sequence. Those are artifacts and are excluded when
comparing.
"""

import networkx as nx
import pytest

from antibody_repertoire.graph.hamming import (
    build_hamming_graph,
    build_hamming_graph_reference,
    collapse_to_unique_cdr3,
    hamming_distance,
)

from conftest import components, genuine_edges, make_df, requires_real_dataset


def assert_equivalent(df):
    """Both implementations agree on nodes, real edges, and components."""
    expected = build_hamming_graph_reference(df)
    actual = build_hamming_graph(df)

    assert set(actual.nodes()) == set(expected.nodes())
    assert genuine_edges(actual) == genuine_edges(expected)
    assert components(actual) == components(expected)
    assert not list(nx.selfloop_edges(actual)), "blocked version must not emit self-loops"


def test_hamming_distance_counts_mismatches():
    assert hamming_distance("AAAA", "AAAA") == 0
    assert hamming_distance("AAAA", "AAAT") == 1
    assert hamming_distance("AAAA", "TTTT") == 4


def test_identical_sequences_form_one_node():
    """Repeated reads of one CDR3 collapse to a single node, no self-loop."""
    df = make_df([
        ('AAAACCCCGGGGTTTTAAAACCCCGGGGTTTT', ['V1'], ['J1']),
        ('AAAACCCCGGGGTTTTAAAACCCCGGGGTTTT', ['V1'], ['J1']),
    ])
    graph = build_hamming_graph(df)
    assert graph.number_of_nodes() == 1
    assert graph.number_of_edges() == 0


def test_close_sequences_are_joined():
    """One mismatch in 32bp is within the 10% threshold."""
    df = make_df([
        ('AAAACCCCGGGGTTTTAAAACCCCGGGGTTTT', ['V1'], ['J1']),
        ('AAAACCCCGGGGTTTTAAAACCCCGGGGTTTA', ['V1'], ['J1']),
    ])
    assert build_hamming_graph(df).number_of_edges() == 1
    assert_equivalent(df)


def test_distant_sequences_are_not_joined():
    df = make_df([
        ('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA', ['V1'], ['J1']),
        ('TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT', ['V1'], ['J1']),
    ])
    assert build_hamming_graph(df).number_of_edges() == 0
    assert_equivalent(df)


def test_threshold_ratio_is_honoured():
    """A wider threshold admits an edge the default rejects."""
    df = make_df([
        ('AAAACCCCGGGGTTTTAAAACCCCGGGGTTTT', ['V1'], ['J1']),
        ('AAAACCCCGGGGTTTTAAAACCCCGGGGTTAA', ['V1'], ['J1']),
    ])
    assert build_hamming_graph(df, threshold_ratio=0.0).number_of_edges() == 0
    assert build_hamming_graph(df, threshold_ratio=0.5).number_of_edges() == 1


def test_different_lengths_never_joined():
    df = make_df([
        ('AAAACCCCGGGGTTTTAAAACCCCGGGGTTTT', ['V1'], ['J1']),
        ('AAAACCCCGGGGTTTTAAAACCCCGGGGTTT', ['V1'], ['J1']),
    ])
    assert build_hamming_graph(df).number_of_edges() == 0
    assert_equivalent(df)


def test_gene_mismatch_blocks_edge():
    """Near-identical CDR3s stay apart when V and J calls do not overlap."""
    df = make_df([
        ('AAAACCCCGGGGTTTTAAAACCCCGGGGTTTT', ['V1'], ['J1']),
        ('AAAACCCCGGGGTTTTAAAACCCCGGGGTTTA', ['V2'], ['J2']),
    ])
    assert build_hamming_graph(df).number_of_edges() == 0
    assert_equivalent(df)


def test_partial_gene_overlap_allows_edge():
    """A single shared V and J call satisfies the intersection predicate."""
    df = make_df([
        ('AAAACCCCGGGGTTTTAAAACCCCGGGGTTTT', ['V1', 'V9'], ['J1', 'J9']),
        ('AAAACCCCGGGGTTTTAAAACCCCGGGGTTTA', ['V9'], ['J9']),
    ])
    assert build_hamming_graph(df).number_of_edges() == 1
    assert_equivalent(df)


def test_gene_calls_union_across_duplicate_reads():
    """Collapsing must union gene calls, not take the first read's.

    The shared CDR3 is reported with V1 on one read and V9 on another. Keeping
    only the first would hide the V9 route and lose the edge to the third
    sequence.
    """
    df = make_df([
        ('AAAACCCCGGGGTTTTAAAACCCCGGGGTTTT', ['V1'], ['J1']),
        ('AAAACCCCGGGGTTTTAAAACCCCGGGGTTTT', ['V9'], ['J9']),
        ('AAAACCCCGGGGTTTTAAAACCCCGGGGTTTA', ['V9'], ['J9']),
    ])
    graph = build_hamming_graph(df)
    assert graph.number_of_nodes() == 2
    assert graph.number_of_edges() == 1
    assert_equivalent(df)


def test_collapse_unions_gene_calls():
    records = make_df([
        ('AAAA', ['V1'], ['J1']),
        ('AAAA', ['V9'], ['J9']),
    ]).to_dict('index')
    collapsed = collapse_to_unique_cdr3(records)
    assert len(collapsed) == 1
    assert collapsed['AAAA']['v_call'] == {'V1', 'V9'}
    assert collapsed['AAAA']['j_call'] == {'J1', 'J9'}


def test_transitive_lineage_forms_one_component():
    """Chained near-neighbours collapse into a single lineage."""
    df = make_df([
        ('AAAACCCCGGGGTTTTAAAACCCCGGGGTTTT', ['V1'], ['J1']),
        ('AAAACCCCGGGGTTTTAAAACCCCGGGGTTTA', ['V1'], ['J1']),
        ('AAAACCCCGGGGTTTTAAAACCCCGGGGTTAA', ['V1'], ['J1']),
    ])
    assert nx.number_connected_components(build_hamming_graph(df)) == 1
    assert_equivalent(df)


def test_isolated_sequences_still_appear_as_nodes():
    df = make_df([
        ('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA', ['V1'], ['J1']),
        ('TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT', ['V2'], ['J2']),
        ('GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG', ['V3'], ['J3']),
    ])
    graph = build_hamming_graph(df)
    assert graph.number_of_nodes() == 3
    assert graph.number_of_edges() == 0
    assert_equivalent(df)


@pytest.mark.slow
@requires_real_dataset
def test_matches_reference_on_real_dataset():
    """End-to-end equivalence on the full 10k-read repertoire."""
    from antibody_repertoire.io.igblast import load_igblast_results
    from conftest import REAL_DATASET

    assert_equivalent(load_igblast_results(REAL_DATASET))
