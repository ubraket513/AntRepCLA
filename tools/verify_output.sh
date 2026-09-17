#!/usr/bin/env bash
# Prove the analysis still reproduces its frozen reference output exactly.
#
# The golden files in tests/golden/ were produced by the original Python
# implementation, before the port, and are committed. This script runs the
# binary over the same inputs and requires byte-identical output. Any
# difference is a failure -- there is no tolerance to tune, because the
# pipeline's output is discrete: an edge is present or it is not, and a missing
# edge silently splits a lineage.
#
# The Python code is gone; these files are what remains of it, and they are the
# reason it could be deleted safely.
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
binary="$root/build/antibody-repertoire"
golden="$root/tests/golden"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

[[ -x "$binary" ]] || { echo "building..."; make -C "$root" >/dev/null; }

status=0
check() {
    local label="$1" expected="$2" actual="$3"
    if diff -q "$expected" "$actual" >/dev/null 2>&1; then
        printf '  \033[32mIDENTICAL\033[0m  %s\n' "$label"
    else
        printf '  \033[31mDIFFERS\033[0m    %s\n' "$label"
        diff "$expected" "$actual" | head -10
        status=1
    fi
}

echo "Real dataset (10,000 reads):"
"$binary" --igblast "$root/data/igblast_results.tsv" \
    --dump-lineages "$work/lineages.tsv" \
    --dump-usage "$work/v_gene_usage.tsv" \
    --weblogo-query "$work/weblogo_query.txt" \
    >"$work/stats.txt" 2>/dev/null
for name in stats.txt lineages.tsv v_gene_usage.tsv weblogo_query.txt; do
    check "$name" "$golden/$name" "$work/$name"
done

# Forty lineages all of size three: every ordering decision in this fixture is
# a tie, which the real dataset cannot exercise because its largest lineage is
# unique. A port that leaves tie-breaking to hash order fails here and only
# here.
echo "Tie-heavy fixture (every lineage the same size):"
"$binary" --igblast "$root/tests/data/tie_heavy.tsv" \
    --dump-lineages "$work/th_lineages.tsv" \
    --dump-usage "$work/th_usage.tsv" \
    >"$work/th_stats.txt" 2>/dev/null
check "tie_heavy_stats.txt"    "$golden/tie_heavy_stats.txt"    "$work/th_stats.txt"
check "tie_heavy_lineages.tsv" "$golden/tie_heavy_lineages.tsv" "$work/th_lineages.tsv"
check "tie_heavy_usage.tsv"    "$golden/tie_heavy_usage.tsv"    "$work/th_usage.tsv"

# Thread count must not affect the result: blocks are independent and the
# lineage ordering is canonical, so the output is worker-count invariant.
echo "Worker-count invariance:"
for workers in 1 2 4 8; do
    "$binary" --igblast "$root/data/igblast_results.tsv" \
        --dump-lineages "$work/w$workers.tsv" >/dev/null 2>&1
done
check "identical across 1/2/4/8 workers" "$golden/lineages.tsv" "$work/w8.tsv"
for workers in 1 2 4; do
    cmp -s "$work/w8.tsv" "$work/w$workers.tsv" || { echo "  worker count $workers diverged"; status=1; }
done

[[ $status -eq 0 ]] && echo "All outputs identical to the frozen reference." \
                    || echo "EQUIVALENCE FAILED"
exit $status
