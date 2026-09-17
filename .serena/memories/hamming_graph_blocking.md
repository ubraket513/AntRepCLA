# Hamming graph: blocking-key optimization

## Result
`hamming_blocked.init()` replaces the O(n^2) pair scan in `hamming.init()`.
On the 10k-read dataset: **124s -> 0.28s (~390x)**, 49,995,000 pairs -> 33,151.

## Why it is exact
The edge predicate needs (a) equal `cdr_length`, (b) non-empty `v_call`
intersection, (c) non-empty `j_call` intersection, (d) hamming <= 10% length.
(a)-(c) are equality constraints, so they work as blocking keys: a pair can
only be joined if it shares a `(cdr_length, v_gene, j_gene)` triple. Sequences
with multiple calls are indexed under every triple they carry, so nothing is
missed. Only (d) is tested inside buckets.

## Two subtleties that must not be undone
1. **Union gene calls when collapsing reads to unique CDR3s.** 33 CDR3 groups
   in the real data are reported with *different* v_call sets on different
   reads. Taking the first read's calls would silently drop edges.
   Covered by `test_gene_calls_union_across_duplicate_reads`.
2. **Candidate pairs must be deduped.** A pair sharing several blocks would
   otherwise be compared repeatedly (1.37x inflation measured).

## Self-loop discrepancy (expected, is a fix)
Original emits 739 self-loops: it iterates over *reads*, so two reads with an
identical CDR3 produce `add_edge(x, x)`. Nodes are CDR3 strings, so these are
artifacts. Blocked version emits none. Genuine edges match exactly
(12033 old = 11294 new + 739 self-loops); all 623 components identical.

## Baseline for regression
`cache/HammingGraph.pkl` = 2477 nodes, 12033 edges (incl. self-loops),
623 components, largest 5 = [103, 99, 85, 65, 50], 240 singletons, 37 with >=10.
Data: 10000 reads, 2477 unique CDR3s (4x read redundancy).

## Env
No system pandas/networkx. Use `.venv` (uv, py3.12): `.venv/bin/python -m pytest tests/ -q`.
Real-data test takes ~107s because it runs the O(n^2) reference for comparison.
