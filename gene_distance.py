#!/usr/bin/env python3
"""
This script produces the plot showing the distance change between two orthogroups across all genomes from config. 
It is helpful for validating syntenically linked genes found by syn_blocks and similar scripts. Or just to validate what we visually see on ribbon plot. 

    python gene_distance.py config.toml GENE1 GENE2 [out.png]

"""

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import functions

MB = 1_000_000 # I prefer it displaying in MB

if len(sys.argv) < 4:
    print("usage: python gene_distance.py config.toml GENE1 GENE2 [out.png]")
    sys.exit(1)

config_path, gene1, gene2 = sys.argv[1], sys.argv[2], sys.argv[3]

out = "gene_distance.png"
if len(sys.argv) > 4:
    out = sys.argv[4]


cfg = functions.read_config(config_path)
chain = cfg["chain"]



if cfg["orthogroups"]:
    chr_sizes, links_by_pair = functions.load_orthofinder(cfg, chain)
else:
    chr_sizes = functions.read_chr_sizes(cfg["chr_sizes"])
    links_by_pair = functions.load_links(cfg, chr_sizes, chain)

sides = []



for (a, b), df in links_by_pair.items():
    sides.append(pd.DataFrame({"species": a, "gene": df["gene1"], "chr": df["chr1"], "pos": df["pos1"]}))
    sides.append(pd.DataFrame({"species": b, "gene": df["gene2"], "chr": df["chr2"], "pos": df["pos2"]}))

genes = pd.concat(sides).drop_duplicates(["species", "gene"]).set_index(["species", "gene"])

print(f"{'species':<16}{'chr ' + gene1:>20}{'chr ' + gene2:>20}{'distance (Mb)':>16}")

rows = []
points = []



for sp in chain:
    label = cfg["labels"].get(sp, sp)
    here1 = (sp, gene1) in genes.index
    here2 = (sp, gene2) in genes.index

    chr1 = genes.loc[(sp, gene1), "chr"] if here1 else "-"
    chr2 = genes.loc[(sp, gene2), "chr"] if here2 else "-"

    distance = "-"
    if here1 and here2 and chr1 == chr2:
        mb = abs(genes.loc[(sp, gene1), "pos"] - genes.loc[(sp, gene2), "pos"]) / MB
        distance = f"{mb:.2f}"
        points.append((mb, chr_sizes[sp][chr1] / MB, label))

    rows.append([label, chr1, chr2, distance])
    print(f"{label:<16}{chr1:>20}{chr2:>20}{distance:>16}")

fig, (ax_table, ax_scatter) = plt.subplots(
    1, 2, figsize=(11, 0.45 * len(chain) + 2.5), gridspec_kw={"width_ratios": [1, 1.3]})

ax_table.axis("off")
table = ax_table.table(cellText=rows, colLabels=["species", gene1, gene2, "distance (Mb)"],
                       loc="center", cellLoc="center")
table.auto_set_font_size(False)
table.set_fontsize(9)
table.scale(1, 1.4)
ax_table.set_title("chromosome of each gene", fontsize=10)

for distance, chr_size, label in points:
    ax_scatter.scatter(distance, chr_size, color="#377eb8", zorder=3)
    ax_scatter.annotate(label, (distance, chr_size), fontsize=8,xytext=(4, 4), textcoords="offset points")

ax_scatter.set_xlabel(f"distance {gene1}-{gene2} (Mb)")
ax_scatter.set_ylabel("chromosome size (Mb)")
ax_scatter.set_title("chromosome size vs. gene distance", fontsize=10)
ax_scatter.margins(0.2)

fig.tight_layout()
fig.savefig(out, dpi=150, bbox_inches="tight")
print(f"saved {out}")
