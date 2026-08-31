#!/usr/bin/env python3
"""Dot plot of orthologs between one chromosome of each species. For finding Ortho-Chr use heatmap.R. 
Usage:
    python dotplot.py config.toml SPECIES_A SPECIES_B CHR_A CHR_B [out.png]

Each dot is one orthologous gene pair, at its position on CHR_A (x) and CHR_B (y). There is also spearman correlation which suppose to indicate level of colinearity. 
"""

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import functions

if len(sys.argv) < 6:
    print("usage: python dotplot.py config.toml SPECIES_A SPECIES_B CHR_A CHR_B [out.png]")
    sys.exit(1)

config_path = sys.argv[1]
sp_a = sys.argv[2]
sp_b = sys.argv[3]
chr_a = sys.argv[4]
chr_b = sys.argv[5]
out = "dotplot.png"
if len(sys.argv) > 6:
    out = sys.argv[6]

cfg = functions.read_config(config_path)
# in case if Orthogroups/GFFs were obtained with miniprot/AGAT pipeline there was a problem with readability, including upper/lowercase in CDS. 
  # Thats how this block of code was introduced. 
if cfg["orthogroups"]:
    chr_sizes, links_by_pair = functions.load_orthofinder(cfg, cfg["chain"])
else:
    chr_sizes = functions.read_chr_sizes(cfg["chr_sizes"])
    links_by_pair = functions.load_links(cfg, chr_sizes, cfg["chain"])

if (sp_a, sp_b) not in links_by_pair:
    print(f"no data for {sp_a} vs {sp_b}; you should pass them in the very same order they are appearing in the config chain")
    sys.exit(1)
df = links_by_pair[(sp_a, sp_b)]

# keep only the orthologs that are on two chosen Chr
pair = df[(df["chr1"] == chr_a) & (df["chr2"] == chr_b)]

# colinearity
# +1 colinear, -1 inverted, ~0 scrambled
rho = pair["pos1"].rank().corr(pair["pos2"].rank())

fig, ax = plt.subplots(figsize=(7, 7))
ax.scatter(pair["pos1"], pair["pos2"], s=8, color="#377eb8")
ax.set_xlabel(f"{sp_a}  {chr_a}")
ax.set_ylabel(f"{sp_b}  {chr_b}")
ax.set_title(f"colinearity (Spearman) = {rho:+.2f}")
fig.savefig(out, dpi=150, bbox_inches="tight")
print(f"saved {out}  ({len(pair)} orthologs, colinearity {rho:+.2f})")
