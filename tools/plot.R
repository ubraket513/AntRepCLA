#!/usr/bin/env Rscript
#
# Publication figures from the pipeline's TSV output.
#
# This is OPTIONAL and deliberately outside the build: `make run` and
# `make check` never invoke it, so nothing here breaks when R is absent. The
# binary always emits serviceable SVG on its own; this is for when you want
# ggplot2 typography for a paper.
#
#   Rscript tools/plot.R [input_dir] [output_dir]     # defaults: out/ out/
#
# ggplot2 is used when installed and base graphics otherwise, so the script
# works against a bare r-base-core with no install.packages() step.

args <- commandArgs(trailingOnly = TRUE)
in_dir  <- if (length(args) >= 1) args[1] else "out"
out_dir <- if (length(args) >= 2) args[2] else "out"

usage_path   <- file.path(in_dir, "v_gene_usage.tsv")
lineage_path <- file.path(in_dir, "lineages.tsv")

for (p in c(usage_path, lineage_path)) {
  if (!file.exists(p)) {
    stop(sprintf("%s not found. Run `make run` first, or pass the directory: Rscript tools/plot.R <dir>", p))
  }
}
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

usage <- read.delim(usage_path, stringsAsFactors = FALSE)
usage <- usage[order(-usage$reads), ]
usage$v_gene <- factor(usage$v_gene, levels = usage$v_gene)

# Lineage sizes. n_unique_cdr3 counts distinct CDR3s, n_reads counts the reads
# collapsed onto them; the gap between the two is the sequencing-error
# tolerance at work, so both are worth showing.
lineages <- read.delim(lineage_path, stringsAsFactors = FALSE,
                       colClasses = c("integer", "integer", "integer", "character"))

has_ggplot <- requireNamespace("ggplot2", quietly = TRUE)

if (has_ggplot) {
  library(ggplot2)

  fill_scale <- if (requireNamespace("viridisLite", quietly = TRUE)) {
    scale_fill_viridis_d(guide = "none")
  } else {
    scale_fill_grey(guide = "none")
  }

  p1 <- ggplot(usage, aes(x = v_gene, y = reads, fill = v_gene)) +
    geom_col(colour = "grey20", linewidth = 0.2) +
    fill_scale +
    labs(title = "V Gene Usage Statistics", x = "V Genes", y = "Read Count") +
    theme_minimal(base_size = 12) +
    theme(
      plot.title  = element_text(face = "bold", hjust = 0.5, colour = "#141482"),
      axis.title  = element_text(face = "bold"),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                 face = "bold", colour = "#2f7d32", size = 7),
      axis.text.y = element_text(face = "bold", colour = "#14579e"),
      panel.grid.major.x = element_blank()
    )
  ggsave(file.path(out_dir, "v_gene_usage_ggplot.pdf"), p1, width = 14, height = 8)

  # Rank-abundance: lineage size against rank, both axes logarithmic. This is
  # the figure the binary does not produce, and it is the one that shows the
  # repertoire is dominated by a few large lineages.
  ranked <- data.frame(rank = seq_len(nrow(lineages)),
                       reads = lineages$n_reads,
                       cdr3s = lineages$n_unique_cdr3)
  p2 <- ggplot(ranked, aes(x = rank)) +
    geom_line(aes(y = reads, colour = "reads"), linewidth = 0.7) +
    geom_line(aes(y = cdr3s, colour = "unique CDR3s"), linewidth = 0.7) +
    scale_x_log10() + scale_y_log10() +
    scale_colour_manual(values = c("reads" = "#14579e", "unique CDR3s" = "#2f7d32"),
                        name = NULL) +
    labs(title = "Clonal Lineage Rank Abundance",
         x = "Lineage rank (largest first)", y = "Count") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, colour = "#141482"),
          axis.title = element_text(face = "bold"),
          legend.position = "top")
  ggsave(file.path(out_dir, "lineage_rank_abundance.pdf"), p2, width = 9, height = 6)

  message("wrote ", out_dir, "/v_gene_usage_ggplot.pdf and lineage_rank_abundance.pdf")
} else {
  message("ggplot2 not installed; falling back to base graphics.")

  pdf(file.path(out_dir, "v_gene_usage_base.pdf"), width = 14, height = 8)
  par(mar = c(9, 5, 4, 2))
  barplot(usage$reads,
          names.arg = as.character(usage$v_gene),
          las = 2, cex.names = 0.6,
          col = hcl.colors(nrow(usage), "viridis"),
          border = "grey20",
          main = "V Gene Usage Statistics",
          ylab = "Read Count")
  mtext("V Genes", side = 1, line = 7, font = 2)
  grid(nx = NA, ny = NULL, lty = 2)
  dev.off()

  pdf(file.path(out_dir, "lineage_rank_abundance_base.pdf"), width = 9, height = 6)
  plot(seq_len(nrow(lineages)), lineages$n_reads, log = "xy", type = "l",
       col = "#14579e", lwd = 2,
       main = "Clonal Lineage Rank Abundance",
       xlab = "Lineage rank (largest first)", ylab = "Count")
  lines(seq_len(nrow(lineages)), lineages$n_unique_cdr3, col = "#2f7d32", lwd = 2)
  legend("topright", legend = c("reads", "unique CDR3s"),
         col = c("#14579e", "#2f7d32"), lwd = 2, bty = "n")
  dev.off()

  message("wrote ", out_dir, "/v_gene_usage_base.pdf and lineage_rank_abundance_base.pdf")
}
