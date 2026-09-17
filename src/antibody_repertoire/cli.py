"""Command-line entry point."""

from __future__ import annotations

import argparse
from pathlib import Path

from antibody_repertoire.config import Config, Paths
from antibody_repertoire.logging import configure_logging
from antibody_repertoire.pipeline import run_analysis


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="antibody-repertoire",
        description=(
            "Clonal lineage analysis of antibody repertoire sequencing data, "
            "from NCBI IgBLAST output."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--igblast",
        type=Path,
        required=True,
        help="Path to the NCBI IgBLAST results TSV (-outfmt 19).",
    )
    parser.add_argument(
        "--cache",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Reuse cached intermediates instead of recomputing.",
    )
    parser.add_argument(
        "--usage-plot",
        type=str,
        default=None,
        help="Basename for the V gene usage plot, written under plots/.",
    )
    parser.add_argument(
        "--weblogo-query",
        type=str,
        default=None,
        help="Basename for the WebLogo query file, written under data/.",
    )
    parser.add_argument(
        "--stats-plot",
        type=str,
        default="clonal_lineage_stats",
        help="Basename for the lineage statistics plot, written under plots/.",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Display plots interactively; ignored where there is no display.",
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=None,
        help="Project root for data/, cache/, logs/ and plots/.",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Console log verbosity.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    paths = Paths(root=args.root) if args.root else Paths()
    config = Config(paths=paths, use_cache=args.cache)
    config.paths.ensure_dirs()

    configure_logging(config.paths.logs, level=args.log_level)

    result = run_analysis(igblast_path=args.igblast, config=config)

    # Plots are imported here so that a headless run never needs matplotlib
    # unless it is actually going to draw something.
    from antibody_repertoire.viz.plots import plot_lineage_stats, plot_v_gene_usage

    for metric, value in result.stats.items():
        print(f"{value:>8,}  {metric}")

    if args.stats_plot:
        plot_lineage_stats(
            result.stats,
            config.paths.plots / f"{args.stats_plot}.png",
            show=args.show,
        )

    if args.usage_plot:
        plot_v_gene_usage(
            result.v_gene_usage,
            config.paths.plots / f"{args.usage_plot}.png",
            show=args.show,
        )

    if args.weblogo_query:
        out = config.paths.data / f"{args.weblogo_query}.txt"
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(result.weblogo_query)
        print(f"\nWebLogo query written to {out}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
