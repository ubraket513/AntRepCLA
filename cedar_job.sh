#!/bin/bash
#SBATCH --mem-per-cpu=2.0G
#SBATCH --time=00:08:00
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48

# Abort on the first failure rather than pressing on with a half-built env.
set -euo pipefail

module load python/3.12

# The virtualenv must be created and activated in *this* shell. Doing it inside
# an `srun ... << EOF` heredoc would confine the activation to that subshell,
# leaving the analysis step below running against the bare system interpreter.
virtualenv --no-download "$SLURM_TMPDIR/env"
source "$SLURM_TMPDIR/env/bin/activate"

pip install --no-index --upgrade pip
pip install --no-index -r requirements.txt
pip install --no-deps -e .

# Define input and output paths
QUERY="data/PRJNA324093_Dnr4_10k.fasta"
OUTPUT="data/igblast_results.tsv"
DB_DIR="database"
USAGE_PLOT="usage_plot"
WEBLOGO_QUERY="weblogo_query"
AUXILIARY="optional_file/human_gl.aux"
CORES=${SLURM_CPUS_PER_TASK:-48}

cd ncbi-igblast-1.22.0

# run ncbi-igblast-1.22.0
bin/igblastn \
    -germline_db_V "$DB_DIR/my_seq_V" \
    -germline_db_D "$DB_DIR/my_seq_D" \
    -germline_db_J "$DB_DIR/my_seq_J" \
    -organism human \
    -domain_system imgt \
    -query "../$QUERY" \
    -outfmt 19 \
    -out "../$OUTPUT" \
    -auxiliary_data "$AUXILIARY" \
    -num_threads "$CORES"

cd ..

# run analysis script
antibody-repertoire --igblast "$OUTPUT" --no-cache --usage-plot "$USAGE_PLOT" --weblogo-query "$WEBLOGO_QUERY"
