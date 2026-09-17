"""Reading and normalising NCBI IgBLAST output."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from loguru import logger

from antibody_repertoire.config import GENE_CALL_COLUMNS, IGBLAST_COLUMNS


def load_igblast_results(filepath: str | Path) -> pd.DataFrame:
    """Load IgBLAST AIRR output and normalise it for lineage analysis.

    Three transformations are applied:

    * allele suffixes are stripped from gene calls, so that ``IGHV3-23*01``
      and ``IGHV3-23*02`` are treated as the same gene;
    * the CDR3 length is derived from the reported start/end coordinates;
    * comma-separated gene calls become sets, since a read may be assigned
      several equally plausible genes and lineage membership only requires
      that two reads share one.

    Args:
        filepath: Path to the IgBLAST results TSV (``-outfmt 19``).

    Returns:
        A frame with the CDR3 sequence, its length, and set-valued V/D/J calls.

    Raises:
        FileNotFoundError: If ``filepath`` does not exist.
    """
    filepath = Path(filepath)
    logger.info(f"Loading IgBLAST results from {filepath}")

    if not filepath.is_file():
        raise FileNotFoundError(f"Input file '{filepath}' does not exist.")

    df = pd.read_csv(filepath, sep='\t', usecols=IGBLAST_COLUMNS)

    for col in GENE_CALL_COLUMNS:
        df[col] = df[col].str.replace(r"\*\d+", '', regex=True)

    df['cdr_length'] = df['cdr3_end'] - df['cdr3_start'] + 1
    df = df.drop(columns=['cdr3_start', 'cdr3_end'])

    for col in GENE_CALL_COLUMNS:
        df[col] = df[col].apply(lambda x: set(str(x).split(",")) if pd.notnull(x) else set())

    logger.info(f"Loaded {len(df):,} sequences ({df['cdr3'].nunique():,} unique CDR3s).")
    return df
