"""Clonal lineage extraction and statistics."""

from __future__ import annotations


import networkx as nx
import pandas as pd
from loguru import logger


def get_clonal_lineages(graph: nx.Graph) -> list[set[str]]:
    """Return connected components as lineages, largest first."""
    lineages = sorted(nx.connected_components(graph), key=len, reverse=True)
    logger.info(f"Identified {len(lineages):,} clonal lineages.")
    return lineages


def assign_lineage_ids(
    igblast_result: pd.DataFrame,
    lineages: list[set[str]],
) -> pd.Series:
    """Map each read to the index of the lineage its CDR3 belongs to.

    Lineages are identified by their position in ``lineages`` (0 being the
    largest). An integer label is used rather than the lineage's member set
    because ``tuple(set)`` ordering is not determined by content -- two equal
    sets can produce different tuples -- which makes set-derived keys unsafe
    to compare or group on.

    Returns:
        An integer Series aligned to ``igblast_result``; reads whose CDR3 is in
        no lineage are ``NA``.
    """
    cdr3_to_lineage = {
        cdr3: index
        for index, lineage in enumerate(lineages)
        for cdr3 in lineage
    }
    return igblast_result['cdr3'].map(cdr3_to_lineage).astype('Int64')


def get_lineage_stats(
    igblast_result: pd.DataFrame,
    lineages: list[set[str]],
    min_lineage_size: int = 10,
) -> dict[str, int]:
    """Summarise the lineage structure of a repertoire.

    Counts are reported at two granularities, which differ because a single
    CDR3 is typically observed on many reads: unique CDR3s measure the
    diversity within a lineage, read counts measure its abundance.

    Args:
        igblast_result: Normalised IgBLAST frame.
        lineages: Lineages, largest first, from :func:`get_clonal_lineages`.
        min_lineage_size: Cutoff for the "well represented" tallies.

    Returns:
        A dict of human-readable metric names to counts.
    """
    lineage_ids = assign_lineage_ids(igblast_result, lineages)

    largest = lineages[0] if lineages else set()
    reads_per_lineage = lineage_ids.value_counts()

    stats = {
        'number of clonal lineages':
            len(lineages),
        'number of unique CDR3s in the largest clonal lineage':
            len(largest),
        'number of sequences in the largest clonal lineage':
            int((lineage_ids == 0).sum()) if lineages else 0,
        f'number of lineages represented by at least {min_lineage_size} unique CDR3s':
            sum(1 for lineage in lineages if len(lineage) >= min_lineage_size),
        f'number of lineages represented by at least {min_lineage_size} sequences':
            int((reads_per_lineage >= min_lineage_size).sum()),
    }

    logger.info("Computed clonal lineage statistics.")
    return stats


def get_v_gene_usage(
    igblast_result: pd.DataFrame,
    graph: nx.Graph,
) -> dict[str, int]:
    """Count reads per V gene, restricted to CDR3s present in the graph.

    A read contributes to every V gene it was assigned, so the counts sum to
    more than the number of reads when calls are ambiguous.
    """
    in_graph = set(graph.nodes())

    usage: dict[str, int] = {
        gene: 0
        for calls in igblast_result['v_call']
        for gene in calls
    }

    for v_call, cdr3 in zip(igblast_result['v_call'], igblast_result['cdr3']):
        if cdr3 in in_graph:
            for gene in v_call:
                usage[gene] += 1

    logger.info(f"Computed usage across {len(usage):,} V genes.")
    return usage


def get_cdr3_aa_from_largest_lineage(
    igblast_result: pd.DataFrame,
    lineages: list[set[str]],
) -> str:
    """Build a WebLogo query from the largest lineage's CDR3 amino acids.

    WebLogo expects equal-length sequences, so sequences are padded to the
    mean length with '*'. Internal gaps use the same character.

    Returns:
        Newline-separated amino acid sequences, or an empty string if there
        are no lineages.
    """
    if not lineages:
        logger.warning("No lineages available; returning an empty query.")
        return ""

    largest = pd.DataFrame(sorted(lineages[0]), columns=['cdr3'])
    merged = pd.merge(igblast_result, largest, how='inner', on='cdr3')

    aa_seqs = merged['cdr3_aa'].str.replace(' ', '*').tolist()
    if not aa_seqs:
        logger.warning("Largest lineage has no amino acid sequences.")
        return ""

    # Pad to the longest sequence, not the mean: padding to the mean leaves
    # every above-average sequence unpadded and the output ragged, which
    # WebLogo rejects.
    target_length = max(map(len, aa_seqs))
    aa_seqs = [seq.ljust(target_length, '*') for seq in aa_seqs]

    logger.info(f"Extracted {len(aa_seqs):,} CDR3 amino acid sequences.")
    return '\n'.join(aa_seqs)
