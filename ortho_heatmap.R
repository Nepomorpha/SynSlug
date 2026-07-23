#!/usr/bin/env Rscript
# This is a supplementary script aim to help you identify orthologous Chr for further analysis. You would firstly need to run export.py.
# Usage: 
#   Rscript heatmap.R A_B_matrix.tsv [out.png]
#
# Rows and columns are clustered so orthologous Chr arranges. 

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {cat("usage: Rscript heatmap.R A_B_matrix.tsv [out.png]\n") quit(status = 1)}

infile <- args[1]
out <- if (length(args) >= 2) args[2] else "heatmap.png"

m <- as.matrix(read.table(infile, header = TRUE, row.names = 1, sep = "\t"))

png(out, width = 700, height = 600)
heatmap(m, scale = "none", margins = c(6, 6))   # scale = "row" instead for per-chromosome normalisation
dev.off()
cat("saved", out, "\n")
