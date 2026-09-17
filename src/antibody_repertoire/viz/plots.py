"""Plots summarising lineage structure and gene usage."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import seaborn as sns
from loguru import logger

from antibody_repertoire.viz.display import has_display  # selects the backend

import matplotlib.pyplot as plt  # noqa: E402  (must follow backend selection)


def _finish(fig_path: Path | None, show: bool) -> None:
    """Save, optionally show, and always close the current figure."""
    plt.tight_layout()

    if fig_path is not None:
        fig_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(fig_path, dpi=150)
        logger.info(f"Plot saved to {fig_path}")

    if show:
        if has_display():
            plt.show()
        else:
            logger.info("No display available; skipping interactive show.")

    plt.close()


def plot_v_gene_usage(
    usage: dict[str, int],
    output_path: Path | None = None,
    show: bool = False,
) -> None:
    """Bar plot of read counts per V gene, most used first.

    Raises:
        ValueError: If ``usage`` is empty.
    """
    if not usage:
        raise ValueError("usage must be a non-empty mapping of gene to count.")

    ordered = dict(sorted(usage.items(), key=lambda item: item[1], reverse=True))

    plt.figure(figsize=(14, 8))
    sns.barplot(
        x=list(ordered.keys()),
        y=list(ordered.values()),
        hue=list(ordered.keys()),
        palette="viridis",
        legend=False,
    )

    plt.title('V Gene Usage Statistics', fontsize=16, fontweight='bold')
    plt.xlabel('V Genes', fontsize=14, fontweight='bold')
    plt.ylabel('Read Count', fontsize=14, fontweight='bold')
    plt.xticks(rotation=90, fontsize=10)
    plt.yticks(fontsize=10)
    plt.grid(axis='y', linestyle='--', alpha=0.7)

    _finish(output_path, show)


def plot_lineage_stats(
    stats: dict[str, int],
    output_path: Path | None = None,
    show: bool = False,
) -> None:
    """Horizontal bar chart of the lineage summary metrics."""
    if not stats:
        raise ValueError("stats must be a non-empty mapping.")

    df = pd.DataFrame(list(stats.items()), columns=['Metric', 'Value'])

    plt.figure(figsize=(12, 6))
    sns.barplot(
        data=df,
        x='Value',
        y='Metric',
        hue='Metric',
        palette='viridis',
        legend=False,
        orient='h',
    )

    plt.title('Clonal Lineage Statistics', fontsize=16, fontweight='bold')
    plt.xlabel('Value', fontsize=14, fontweight='bold')
    plt.ylabel('')

    for index, value in enumerate(df['Value']):
        plt.text(value, index, f' {value}', va='center', fontsize=11)

    _finish(output_path, show)
