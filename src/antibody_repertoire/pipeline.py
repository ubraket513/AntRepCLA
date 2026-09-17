"""End-to-end clonal lineage analysis."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import networkx as nx
import pandas as pd
from loguru import logger

from antibody_repertoire.analysis.lineage import (
    get_cdr3_aa_from_largest_lineage,
    get_clonal_lineages,
    get_lineage_stats,
    get_v_gene_usage,
)
from antibody_repertoire.config import Config
from antibody_repertoire.graph.hamming import build_hamming_graph
from antibody_repertoire.io.cache import load_cache, save_cache
from antibody_repertoire.io.igblast import load_igblast_results


@dataclass
class AnalysisResult:
    """Everything the pipeline produces, for plotting or further analysis."""

    igblast_result: pd.DataFrame
    graph: nx.Graph
    lineages: list[set[str]]
    stats: dict[str, int]
    v_gene_usage: dict[str, int]
    weblogo_query: str


def build_graph(igblast_result: pd.DataFrame, config: Config) -> nx.Graph:
    """Build the Hamming graph, reusing a cached one when allowed."""
    cache_path = config.paths.hamming_graph_cache

    if config.use_cache:
        cached = load_cache(cache_path)
        if cached is not None:
            return cached

    graph = build_hamming_graph(igblast_result, config.hamming_threshold_ratio)
    save_cache(graph, cache_path)
    return graph


def run_analysis(igblast_path: Path | None = None, config: Config | None = None) -> AnalysisResult:
    """Run the full analysis and return its results.

    No files are written apart from caches; rendering and reporting are the
    caller's concern, which keeps this callable from a notebook.
    """
    config = config or Config()
    igblast_path = igblast_path or config.paths.igblast_results

    igblast_result = load_igblast_results(igblast_path)
    graph = build_graph(igblast_result, config)
    lineages = get_clonal_lineages(graph)

    stats = get_lineage_stats(igblast_result, lineages, config.min_lineage_size)
    usage = get_v_gene_usage(igblast_result, graph)
    weblogo_query = get_cdr3_aa_from_largest_lineage(igblast_result, lineages)

    save_cache(usage, config.paths.usage_counts_cache)
    logger.success("Analysis complete.")

    return AnalysisResult(
        igblast_result=igblast_result,
        graph=graph,
        lineages=lineages,
        stats=stats,
        v_gene_usage=usage,
        weblogo_query=weblogo_query,
    )
