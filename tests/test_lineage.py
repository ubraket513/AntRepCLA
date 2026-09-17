"""Lineage extraction, statistics and gene usage."""

import networkx as nx
import pandas as pd

from antibody_repertoire.analysis.lineage import (
    assign_lineage_ids,
    get_cdr3_aa_from_largest_lineage,
    get_clonal_lineages,
    get_lineage_stats,
    get_v_gene_usage,
)

from conftest import make_df, requires_real_dataset


def test_lineages_are_ordered_largest_first():
    graph = nx.Graph()
    graph.add_edges_from([('a', 'b'), ('a', 'c')])   # size 3
    graph.add_edges_from([('d', 'e')])               # size 2
    graph.add_node('f')                              # size 1

    sizes = [len(c) for c in get_clonal_lineages(graph)]
    assert sizes == [3, 2, 1]


def test_assign_lineage_ids_indexes_by_position():
    df = make_df([
        ('AAAA', ['V1'], ['J1']),
        ('TTTT', ['V1'], ['J1']),
        ('GGGG', ['V1'], ['J1']),
    ])
    lineages = [{'AAAA', 'TTTT'}, {'GGGG'}]
    ids = assign_lineage_ids(df, lineages)
    assert list(ids) == [0, 0, 1]


def test_assign_lineage_ids_marks_unassigned_reads_na():
    df = make_df([('AAAA', ['V1'], ['J1']), ('CCCC', ['V1'], ['J1'])])
    ids = assign_lineage_ids(df, [{'AAAA'}])
    assert ids[0] == 0
    assert pd.isna(ids[1])


def test_lineage_ids_are_stable_for_equal_sets():
    """Membership must not depend on set iteration order.

    Building keys with ``tuple(set)`` was order-dependent: two equal sets can
    produce different tuples. Integer positions cannot drift this way.
    """
    df = make_df([('AAAA', ['V1'], ['J1'])])
    forward = {'AAAA', 'TTTT', 'GGGG'}
    backward = set()
    for cdr3 in ['GGGG', 'TTTT', 'AAAA']:
        backward.add(cdr3)

    assert list(assign_lineage_ids(df, [forward])) == list(assign_lineage_ids(df, [backward]))


def test_stats_separate_read_counts_from_unique_cdr3s():
    """The largest lineage has 2 unique CDR3s but 3 supporting reads."""
    df = make_df([
        ('AAAA', ['V1'], ['J1']),
        ('AAAA', ['V1'], ['J1']),
        ('TTTT', ['V1'], ['J1']),
        ('GGGG', ['V2'], ['J2']),
    ])
    lineages = [{'AAAA', 'TTTT'}, {'GGGG'}]
    stats = get_lineage_stats(df, lineages, min_lineage_size=2)

    assert stats['number of clonal lineages'] == 2
    assert stats['number of unique CDR3s in the largest clonal lineage'] == 2
    assert stats['number of sequences in the largest clonal lineage'] == 3
    assert stats['number of lineages represented by at least 2 unique CDR3s'] == 1
    assert stats['number of lineages represented by at least 2 sequences'] == 1


def test_stats_on_empty_lineages():
    df = make_df([('AAAA', ['V1'], ['J1'])])
    stats = get_lineage_stats(df, [])
    assert stats['number of clonal lineages'] == 0
    assert stats['number of sequences in the largest clonal lineage'] == 0


def test_usage_counts_every_assigned_gene():
    """A read with ambiguous calls contributes to each gene it names."""
    df = make_df([
        ('AAAA', ['V1', 'V2'], ['J1']),
        ('TTTT', ['V1'], ['J1']),
    ])
    graph = nx.Graph()
    graph.add_nodes_from(['AAAA', 'TTTT'])

    usage = get_v_gene_usage(df, graph)
    assert usage == {'V1': 2, 'V2': 1}


def test_usage_ignores_reads_absent_from_graph():
    df = make_df([
        ('AAAA', ['V1'], ['J1']),
        ('TTTT', ['V2'], ['J1']),
    ])
    graph = nx.Graph()
    graph.add_node('AAAA')

    assert get_v_gene_usage(df, graph) == {'V1': 1, 'V2': 0}


def test_weblogo_query_pads_to_equal_length():
    df = pd.DataFrame([
        {'cdr3': 'AAAA', 'cdr3_aa': 'CAR', 'cdr_length': 4,
         'v_call': {'V1'}, 'j_call': {'J1'}},
        {'cdr3': 'TTTT', 'cdr3_aa': 'CARWG', 'cdr_length': 4,
         'v_call': {'V1'}, 'j_call': {'J1'}},
    ])
    lines = get_cdr3_aa_from_largest_lineage(df, [{'AAAA', 'TTTT'}]).split('\n')

    assert len(lines) == 2
    assert len({len(line) for line in lines}) == 1, "WebLogo needs equal lengths"


def test_weblogo_query_empty_without_lineages():
    df = make_df([('AAAA', ['V1'], ['J1'])])
    assert get_cdr3_aa_from_largest_lineage(df, []) == ""


@requires_real_dataset
def test_real_dataset_statistics_are_unchanged():
    """Pin the published figures for the reference repertoire.

    These are the numbers the original implementation produced and that the
    report is written against.
    """
    from antibody_repertoire.graph.hamming import build_hamming_graph
    from antibody_repertoire.io.igblast import load_igblast_results
    from conftest import REAL_DATASET

    df = load_igblast_results(REAL_DATASET)
    lineages = get_clonal_lineages(build_hamming_graph(df))
    stats = get_lineage_stats(df, lineages)

    assert stats == {
        'number of clonal lineages': 623,
        'number of unique CDR3s in the largest clonal lineage': 103,
        'number of sequences in the largest clonal lineage': 501,
        'number of lineages represented by at least 10 unique CDR3s': 37,
        'number of lineages represented by at least 10 sequences': 272,
    }
