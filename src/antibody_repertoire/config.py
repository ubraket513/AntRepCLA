"""Project-wide paths and analysis parameters.

Every filesystem location and tunable threshold lives here, so that a caller
can relocate outputs (a scratch directory on a cluster node, say) by building
one :class:`Config` instead of hunting for string literals.

Paths default to a layout rooted at the current working directory, matching how
the pipeline has always been run from the repository root.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

# Columns read from the IgBLAST AIRR-format (outfmt 19) output.
IGBLAST_COLUMNS = [
    'v_call', 'd_call', 'j_call',
    'cdr3', 'cdr3_aa',
    'cdr3_start', 'cdr3_end',
]

# Gene-call columns whose allele suffixes are stripped and which are split
# into sets of gene names.
GENE_CALL_COLUMNS = ['v_call', 'd_call', 'j_call']


@dataclass(frozen=True)
class Paths:
    """Filesystem layout for inputs, caches and outputs."""

    root: Path = field(default_factory=Path.cwd)

    @property
    def data(self) -> Path:
        return self.root / "data"

    @property
    def cache(self) -> Path:
        return self.root / "cache"

    @property
    def logs(self) -> Path:
        return self.root / "logs"

    @property
    def plots(self) -> Path:
        return self.root / "plots"

    @property
    def igblast_results(self) -> Path:
        return self.data / "igblast_results.tsv"

    @property
    def hamming_graph_cache(self) -> Path:
        return self.cache / "HammingGraph.pkl"

    @property
    def usage_counts_cache(self) -> Path:
        return self.cache / "UsageCounts.pkl"

    def ensure_dirs(self) -> None:
        """Create the writable output directories if they are missing."""
        for directory in (self.cache, self.logs, self.plots, self.data):
            directory.mkdir(parents=True, exist_ok=True)


@dataclass(frozen=True)
class Config:
    """Parameters controlling lineage inference.

    Attributes:
        paths: Filesystem layout.
        hamming_threshold_ratio: Two CDR3s may be joined when their Hamming
            distance is at most this fraction of the CDR3 length.
        min_lineage_size: Cutoff used when reporting how many lineages are
            "well represented", applied to both unique CDR3s and read counts.
        use_cache: Whether expensive intermediates may be loaded from disk.
    """

    paths: Paths = field(default_factory=Paths)
    hamming_threshold_ratio: float = 0.1
    min_lineage_size: int = 10
    use_cache: bool = False
