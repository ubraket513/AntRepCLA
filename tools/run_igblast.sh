#!/usr/bin/env bash
# Run NCBI IgBLAST over a FASTA and produce the AIRR TSV the analysis reads.
#
# This is the surviving half of the old SLURM script: IgBLAST is a separate
# program and is needed whatever the analysis itself is written in.
#
#   tools/run_igblast.sh <query.fasta> <output.tsv> [threads]
set -euo pipefail

query="${1:?usage: run_igblast.sh <query.fasta> <output.tsv> [threads]}"
output="${2:?usage: run_igblast.sh <query.fasta> <output.tsv> [threads]}"
threads="${3:-$(nproc)}"

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
igblast="$root/ncbi-igblast-1.22.0"
[[ -d "$igblast" ]] || { echo "IgBLAST not found at $igblast" >&2; exit 1; }

query="$(readlink -f "$query")"
mkdir -p "$(dirname "$output")"
output="$(readlink -f "$output")"

# igblastn resolves internal_data/ and database/ relative to its own directory,
# so the working directory has to be the IgBLAST tree; paths pointing back out
# are made absolute above rather than relying on a fixed number of "..".
cd "$igblast"
bin/igblastn \
    -germline_db_V database/my_seq_V \
    -germline_db_D database/my_seq_D \
    -germline_db_J database/my_seq_J \
    -organism human \
    -domain_system imgt \
    -query "$query" \
    -outfmt 19 \
    -out "$output" \
    -auxiliary_data optional_file/human_gl.aux \
    -num_threads "$threads"

echo "wrote $output"
