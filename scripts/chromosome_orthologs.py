#!/usr/bin/env python3
# Advanced version of ortho_heatmap.py.
# Aimed to facilitate choice of orthologous chromosomes for further visualisaiton/synteny analysis.
# Compare chromosomes from every species with the focal species.
#
# Usage:
#   python scripts/chromosome_orthologs.py config.toml [out.png]
#
# This uses the OrthoFinder and GFF files taht you would use for normal SynSlug run.
# Unlike the main SynSlug program, this script also needs seaborn and scipy.

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.spatial.distance import cdist

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from functions import read_config
from orthofinder import read_gff_positions, read_orthogroups, _resolve_chr_sizes


if len(sys.argv) < 2:
    print("usage: python scripts/chromosome_orthologs.py config.toml [out.png]")
    sys.exit(1)

config_path = sys.argv[1]
output_file = "orthologous_chromosomes.png"

if len(sys.argv) > 2:
    output_file = sys.argv[2]

cfg = read_config(config_path)
chain = cfg["chain"]
focal = cfg["focal"]

if not cfg["orthogroups"]:
    print("chromosome_orthologs.py needs orthogroups and GFF files in the config")
    sys.exit(1)

gene_og = read_orthogroups(cfg["orthogroups"])

positions = {}
chr_sizes = {}

for species in chain:
    positions[species], chr_sizes[species] = read_gff_positions(cfg["gff"][species])

chr_sizes = _resolve_chr_sizes(cfg, chain, chr_sizes)

rows = []

for species in chain:
    for gene, orthogroup in gene_og[species].items():
        if gene not in positions[species]:
            continue

        chromosome = positions[species][gene][0]

        if species == focal and cfg["focal_chrs"]:
            if chromosome not in cfg["focal_chrs"]:
                continue

        if cfg["min_chr_size"]:
            if chromosome not in chr_sizes[species]:
                continue
            if chr_sizes[species][chromosome] < cfg["min_chr_size"]:
                continue

        rows.append([species, chromosome, orthogroup])

genes = pd.DataFrame(
    rows,
    columns=["species", "chromosome", "orthogroup"],
)

genes["chromosome_id"] = genes["species"] + ":" + genes["chromosome"]

matrix = pd.crosstab(
    genes["chromosome_id"],
    genes["orthogroup"],
)

# Presence or absence stops large paralog families from dominating the comparison.
matrix = matrix > 0

focal_ids = []
other_ids = []

for species in chain:
    for chromosome in chr_sizes[species]:
        chromosome_id = species + ":" + chromosome

        if chromosome_id not in matrix.index:
            continue

        if species == focal:
            focal_ids.append(chromosome_id)
        else:
            other_ids.append(chromosome_id)

focal_matrix = matrix.loc[focal_ids]
other_matrix = matrix.loc[other_ids]

similarity = 1 - cdist(
    other_matrix,
    focal_matrix,
    metric="jaccard",
)

similarity = pd.DataFrame(
    similarity,
    index=other_ids,
    columns=focal_ids,
)

# Choose the focal chromosome that shares the most orthogroups with each row.
shared_orthogroups = other_matrix.astype(int).dot(focal_matrix.astype(int).T)

figure_height = len(other_ids) / 8
if figure_height < 6:
    figure_height = 6

plt.figure(figsize=(10, figure_height))
heatmap = sns.heatmap(
    similarity,
    cmap="Blues",
    vmin=0,
    xticklabels=True,
    yticklabels=True,
    cbar_kws={"label": "Jaccard similarity"},
)

# I also added a small circle that marks the Chr with maximum OGs. Helps see what will be missed if you'd filter for max OGs only without looking at the heatmap.
for row_number, chromosome_id in enumerate(shared_orthogroups.index):
    most_shared = shared_orthogroups.loc[chromosome_id].max()

    if most_shared == 0:
        continue

    for column_number, focal_chromosome in enumerate(shared_orthogroups.columns):
        if shared_orthogroups.loc[chromosome_id, focal_chromosome] == most_shared:
            heatmap.scatter(
                column_number + 0.5,
                row_number + 0.5,
                marker="o",
                facecolors="none",
                edgecolors="red",
                s=20,
            )

plt.xlabel(focal + " chromosome")
plt.ylabel("chromosome from another species")
plt.tight_layout()
plt.savefig(output_file, dpi=150)
plt.close()

print(f"saved {output_file}")
