# Design decisions

Why this code is shaped the way it is. Read this before changing the algorithm,
the output formats, or the dependency set — most of what looks arbitrary here
was measured, and several obvious-looking improvements were tried and reverted.

Last reviewed 2026-09-17.

## Current state

**The pipeline is C++ only.** The Python implementation was ported and then
deleted. `make && make check` builds it and proves it still reproduces the
published results exactly. `README.md` carries the science and the usage.

```bash
make            # build, uses every core
make check      # unit tests + byte-for-byte check against the frozen results
make run        # analyse the bundled dataset into out/
make igblast QUERY=reads.fasta
```

## Why deleting Python was safe

Two independent checks replaced it before it was removed. Do not weaken either.

1. **`build_edges_reference` in `src/hamming.cpp`** compares every pair directly
   — no blocking, no chunk filter, no canonical-block rule. It is the port of
   the Python `build_hamming_graph_reference` oracle and is **not dead code**.
   The fast path is three optimisations deep; each can only ever *drop* an edge,
   and a dropped edge silently splits a lineage with no other symptom. Tests
   require exact agreement on generated repertoires and on the real dataset
   (11,294 edges).
2. **`tests/golden/`** holds output frozen from the Python implementation.
   `tools/verify_output.sh` requires byte-identical results. These files are
   what remains of the Python code and are the reason it could go. Regenerate
   them only deliberately — there is no longer a second implementation to
   regenerate them *from*, so a change here is a change to the contract.

`tests/data/tie_heavy.tsv` exists solely because the real dataset cannot catch a
tie-breaking bug: its largest lineage is unique (103 vs 99), while all 40
lineages in the fixture have size 3. Lineage ordering is canonical —
`(-size, smallest member)` — precisely so a reimplementation can reproduce it.

## Verified against the published report

Every figure in `report/main.tex` is reproduced by the current code: 623
lineages, 501 sequences and 103 unique CDR3s in the largest, 272 lineages with
>= 10 sequences, 37 with >= 10 unique CDR3s, and the top six V genes
(979/964/769/768/732/719).

**One discrepancy in the report itself, not in the code:** its V gene table is
headed "Lineage Count" and the plot axis reads "number of clonal lineages formed
by each of V genes", but the values are **read counts**. 979 lineages is
impossible when there are only 623 in total. The code counts reads, the numbers
match, and `README.md` says "Reads". The `.tex` still says lineages.

## Performance (measured, same machine, best of three)

| | Python | C++ |
|---|---|---|
| 10,000 reads | 1.17 s, 129 MB | 0.03 s, 61 MB |
| 200,000 reads, 1 worker | 18.36 s, 591 MB | 0.72 s, 137 MB |
| 200,000 reads, 8 workers | 11.57 s, 628 MB | 0.47 s, 168 MB |

At 8 workers the run is ~47% parsing, 18% collapse, 25% scan, 9% union-find —
about three quarters serial, so more cores buy little. The scan is no longer the
bottleneck; parsing is. Anyone chasing this further should start there.

## Decisions worth not relitigating

- **Vendored headers, committed on purpose.** `csv.hpp` 5.3.0 (SHA-256 verified),
  ankerl/unordered_dense 5.0.1, doctest. This is what makes `make` work offline
  with no package manager. Note `unordered_dense` `main` is no longer a single
  header (needs `stl.h`), so the vendored copy is the v5.0.1 release directory.
- **No logging library.** spdlog was adopted and removed — 1.2 MB and 7.7 s of
  compile time per TU for ~8 log lines.
- **No graph library.** Union-find is 4.9% of the run; Boost.Graph or igraph
  would be a large dependency to speed up a twentieth of it.
- **No plotting library.** Matplot++ needs CMake to build and **Gnuplot at
  runtime**, which is the dependency class this port exists to remove. R and
  Julia were considered and rejected for the same reason, more so: a second
  language runtime plus its package manager, to draw two bar charts. Julia is
  not even in apt here, and its time-to-first-plot would exceed the entire
  analysis. `src/svg.cpp` writes the charts directly.
- **R is nevertheless available as an opt-in.** `tools/plot.R` + `make plot-r`
  read the TSV output and add a rank-abundance figure the binary does not
  produce. It is not a dependency of `run` or `check`. R 4.5.2 *is* installed on
  this machine (an earlier check in-session wrongly reported otherwise);
  ggplot2 is not, so the script's base-graphics fallback is the path that
  actually runs here, and it was tested.
- **No precompiled header.** Tried for `csv.hpp` and measured *worse*: a 187 MB
  `.gch` that costs more to load (8.3 s) than reparsing the header (7.3 s), plus
  a 5.9 s generation step that serialises the build. Clean build went 11.8 s ->
  35.3 s. It is included by one TU, so incremental builds never pay anyway.
- **IgBLAST being C++ is irrelevant** to any of this — it is invoked as a
  subprocess and never linked against.

## Environment

- `clangd` 22.1.6 at `~/.local/share/clangd_22.1.6`, symlinked into
  `~/.local/bin` (no sudo on this machine). `make compile-commands` writes the
  database it needs.
- `.serena/project.yml` is `languages: [cpp]`.
- `clangd --check` reporting "2 errors" on a clean tree is benign: it probes
  every refactoring action at every token and `ExtractFunction` refuses on
  `break`/`continue`.

## Open items

1. **Nothing regenerates `tests/golden/`.** If the algorithm ever legitimately
   changes, the goldens must be updated by hand and the change justified. That
   is the price of deleting the reference implementation, and it was accepted
   knowingly.
2. **The full BioProject scale-up is still not run.** 100 runs, 23.6M reads, 3
   donors (036: 8.7M, 122: 6.6M, 127: 8.3M). Group **per donor, pooling that
   donor's timepoints** — lineages are donor-private, so pooling across donors
   would invent biologically meaningless ones. Needs a run->donor manifest from
   `SampleName` (`036-S001` -> `036`). Refetch metadata with
   `curl "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc=PRJNA324093"`.
3. **`report/main.tex` still says "Lineage Count"** where it means reads (above).
4. **SIMD on the Hamming inner loop** is the remaining lever if the scan ever
   matters again — but parsing dominates now, so measure before starting.
