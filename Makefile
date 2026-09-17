# Build and run. There is nothing to install and nothing to download: every
# dependency is a vendored header under third_party/, and the only tools needed
# are g++ (C++20, OpenMP) and make.
#
#   make            build
#   make check      unit tests + byte-for-byte check against the frozen results
#   make run        analyse the bundled dataset, writing every output
#   make igblast    FASTA -> AIRR TSV, the step before the analysis
#   make plot-r     optional: extra figures via R, if you have it installed
#   make format     optional: apply .clang-format to src/ and tests/
#
# Builds use every core by default; override with `make -j1` or `NPROC=2 make`.

NPROC    ?= $(shell nproc 2>/dev/null || echo 4)
MAKEFLAGS += -j$(NPROC) --no-print-directory

CXX      ?= g++
CXXFLAGS ?= -std=c++20 -O2 -Wall -Wextra -Wpedantic -fopenmp
INCLUDES  = -Ithird_party -Isrc
LDFLAGS   = -fopenmp

BUILD   = build
BINARY  = $(BUILD)/antibody-repertoire
TESTS   = $(BUILD)/run-tests

SOURCES = $(wildcard src/*.cpp)
LIB_SOURCES = $(filter-out src/main.cpp,$(SOURCES))
LIB_OBJECTS = $(LIB_SOURCES:%.cpp=$(BUILD)/%.o)
TEST_OBJECTS = $(BUILD)/tests/test_main.o $(BUILD)/tests/doctest_main.o

# Inputs and outputs for `make run`.
IGBLAST_TSV ?= data/igblast_results.tsv
OUT         ?= out

.PHONY: all test check verify run igblast plot-r format compile-commands clean help
all: $(BINARY)

$(BINARY): $(LIB_OBJECTS) $(BUILD)/src/main.o
	@$(CXX) $^ -o $@ $(LDFLAGS)
	@echo "built $@"

$(TESTS): $(LIB_OBJECTS) $(TEST_OBJECTS)
	@$(CXX) $^ -o $@ $(LDFLAGS)

# No precompiled header for csv.hpp, though it is by far the slowest include.
# Measured: the .gch comes out at 187 MB and loading it costs *more* than
# reparsing the header (8.3 s vs 7.3 s), while also serialising the build
# behind a 5.9 s generation step. It is included by igblast.cpp alone, so
# incremental builds never pay for it anyway.

$(BUILD)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -MMD -MP -c $< -o $@

test: $(TESTS)
	@./$(TESTS)

# Unit tests, then the frozen-output check. This is the gate that matters.
check: test verify

# Split out from `check` so CI reports the two failure modes separately: a
# broken unit test and a changed pipeline result mean different things.
verify: $(BINARY)
	@./tools/verify_output.sh

run: $(BINARY)
	@mkdir -p $(OUT)
	@./$(BINARY) --igblast $(IGBLAST_TSV) \
		--dump-lineages $(OUT)/lineages.tsv \
		--dump-usage $(OUT)/v_gene_usage.tsv \
		--weblogo-query $(OUT)/weblogo_query.txt \
		--usage-plot $(OUT)/v_gene_usage.svg \
		--stats-plot $(OUT)/lineage_stats.svg
	@echo "outputs in $(OUT)/"

# Optional extra figures through R. Deliberately not a dependency of `run` or
# `check`: the binary's own SVG output never needs R, and this target is the
# only thing in the build that can fail for want of a tool.
plot-r:
	@command -v Rscript >/dev/null 2>&1 || { \
		echo "Rscript not found. This target is optional -- 'make run' already"; \
		echo "wrote $(OUT)/*.svg without it."; exit 1; }
	@Rscript tools/plot.R $(OUT) $(OUT)

# FASTA -> AIRR TSV. QUERY is required; OUT_TSV defaults to the analysis input.
igblast:
	@./tools/run_igblast.sh $(QUERY) $(if $(OUT_TSV),$(OUT_TSV),$(IGBLAST_TSV)) $(NPROC)

# Formatting is not enforced by CI: .clang-format is provided for editors and
# for this target, but the codebase predates it and has not been normalised
# against it, so gating on it would fail for reasons unrelated to any change.
format:
	@command -v clang-format >/dev/null 2>&1 || { \
		echo "clang-format not found; .clang-format is still used by editors."; exit 1; }
	@clang-format -i src/*.cpp src/*.hpp tests/*.cpp
	@echo "formatted src/ and tests/"

# clangd needs the include paths and the standard, which a Makefile does not
# publish. Generated here rather than requiring `bear`.
compile-commands:
	@printf '[\n' > compile_commands.json
	@first=1; for src in $(SOURCES) tests/test_main.cpp tests/doctest_main.cpp; do \
		if [ $$first -eq 0 ]; then printf ',\n' >> compile_commands.json; fi; first=0; \
		printf '  {"directory": "%s", "file": "%s", "command": "%s %s %s -c %s"}' \
			"$(CURDIR)" "$$src" "$(CXX)" "$(CXXFLAGS)" "$(INCLUDES)" "$$src" \
			>> compile_commands.json; \
	done
	@printf '\n]\n' >> compile_commands.json
	@echo "wrote compile_commands.json"

help:
	@awk '/^#/ { sub(/^# ?/, ""); print; next } { exit }' Makefile

clean:
	@rm -rf $(BUILD) compile_commands.json

-include $(shell find $(BUILD) -name '*.d' 2>/dev/null)
