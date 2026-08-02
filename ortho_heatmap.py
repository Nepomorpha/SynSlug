#!/usr/bin/env python3
# This is a supplementary script aim to help you identify orthologous Chr for further analysis. You would firstly need to run export.py.
# Usage: 
#   python ortho_heatmap.py A_B_matrix.tsv [out.png]
#
# Rows and columns are clustered so orthologous Chr arranges. 

# Note: Unlike other scripts this one requires seaborn package.

# imports -------
import sys
import matplotlib
matplotlib.use("Agg")
import pandas as pd
import seaborn as sns

# args ----------
if len(sys.argv) < 2:
    print("usage: python ortho_heatmap.py links.tsv [out.png]")
    sys.exit(1)

link_file = sys.argv[1]
output_file = "heatmap.png"

if len(sys.argv) > 2:
    output_file = sys.argv[2]

# read and plot -------
links = pd.read_csv(link_file, sep="\t")
matrix = pd.crosstab(links["chr1"], links["chr2"])

# core function from seaborn package doing the job
heatmap = sns.clustermap(
    matrix,
    method="complete",
    metric="euclidean",
    cmap="viridis",
    figsize=(7, 6),
)

# save -----
heatmap.savefig(output_file, dpi=150)
print(f"saved {output_file}")