#!/usr/bin/env python3

import os
import re
import sys
import colorsys
import tomllib

# Defaults for everything a config does not set.
DEFAULTS = {
    "links": {},                # [links] table: '<spA>_<spB>' = path, one per adjacent pair
    "chr_sizes": None,          # chromosome sizes TSV (optional in OrthoFinder mode)
    "orthogroups": None,        # OrthoFinder Orthogroups.tsv           (OrthoFinder mode)
    "gff": {},                  # species -> GFF path                   (OrthoFinder mode)
    "max_copies": 0,            # skip orthogroups with > this many located genes in a species (0 = off)
    "tracks": [],               # [[tracks]] list: extra bed/histogram lanes below/above a species
    "track_height": 0.18,       # height of one track lane (track units)
    "bedgraph": {},             # species -> .bedgraph path; colours chromosome bars by signal
    "bedgraph_cmap": "viridis",  # colormap for the bedgraph signal
    "bedgraph_vmin": None,      # signal value mapped to the low end of the colormap (default: data min)
    "bedgraph_vmax": None,      # signal value mapped to the high end (default: data max)
    "bedgraph_label": "",       # tiny legend title
    "sort": "input",            # input | name | name_reverse | size | synteny
    "min_links": 0,             # drop chromosomes with fewer coloured links (0 = keep all)
    "max_chrs": 0,              # keep only the N largest chromosomes per species (0 = all)
    "max_links": 2000,          # max links drawn per species pair
    "link_alpha": 0.2,
    "half_width": 5000,         # bp; each ribbon foot spans pos +/- half_width
    "min_ribbon": 0.0006,       # minimum ribbon foot half-width (axis units), keeps links visible
    "min_bar_width": 0.01,      # minimum chromosome bar width (axis units), keeps tiny chromosomes visible
    "default_color": "#e3e3e3",  # links that don't trace back to the focal species
    "drop_grey": False,         # if true, don't draw links that didn't reach the focal species
    "chr_label": False,
    "chr_label_style": "number",  # number | name | short
    "chr_label_size": 6,
    "chr_label_pos": "above",   # "above" the bar, or "inside" (centred on the bar)
    "chr_label_color": "black",  # chromosome label colour
    "chr_names_tsv": None,      # TSV/CSV: real_name -> display_name for chromosome labels
    "bar_color": "#bdbdbd",     # default chromosome bar colour
    "chr_colors_tsv": None,     # TSV: chromosome -> bar colour (overrides bar_color per chr)
    "scale_bar_mb": 0,          # draw a reference scale bar this many Mb long (0 = off)
    "highlight_chrs": [],       # highlight mode: colour only these chromosomes' synteny, rest grey
    "highlight_links_tsv": None,  # TSV: link id -> colour [-> label]; draws those links on top
    "highlight_legend": False,    # show a legend of the highlighted-link colours
    "highlight_legend_title": "",
    "grey_alpha": None,         # alpha for grey links (default: same as link_alpha)
    "bar_height": 0.05,         # thickness of a chromosome bar (track units)
    "space": 0.004,             # gap between chromosomes, as a fraction of the widest genome
    "fig_width": 12.0,          # figure width in inches
    "fig_height": None,         # figure height in inches (auto from species count if unset)
    "dpi": 300,
    "labels": {},
    "focal_chrs": [],
    "chr_order": {},
    "colors_tsv": None,
    "out": None,}

# Drawing

def main_ribbon_painter(top_chr_span, bot_chr_span, y_top, y_bot, rgba):
# see `drawing_studies.ipynb` for detailed info how I come to this function. 
    xa1, xa2 = top_chr_span
    xb1, xb2 = bot_chr_span
    y_mid = (y_top + y_bot) / 2.0
    verts = [(xa1, y_top), (xa1, y_mid), (xb1, y_mid), (xb1, y_bot), (xb2, y_bot), (xb2, y_mid), (xa2, y_mid), (xa2, y_top), (xa1, y_top)]
    codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.LINETO, Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.CLOSEPOLY]
    return PathPatch(Path(verts, codes), facecolor=rgba, edgecolor="none", linewidth=0)
