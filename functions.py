#!/usr/bin/env python3

import os
import re
import sys
import tomllib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.cm import ScalarMappable
from matplotlib.collections import PatchCollection
from matplotlib.colors import to_hex, to_rgba, Normalize
from matplotlib.path import Path
from matplotlib.patches import PathPatch, Rectangle, Patch, Arc, Polygon

from loaders import read_chr_sizes, read_colors, read_bedgraph, read_name_map, read_highlight_links, load_links, load_tracks
from orthofinder import load_orthofinder

# Defaults for everything a config does not set.
DEFAULTS = {
    "links": {},                # [links] table: '<spA>_<spB>' = path, one per adjacent pair
    "chr_sizes": None,          # chromosome sizes TSV (optional in OrthoFinder mode)
    "orthogroups": None,        # OrthoFinder Orthogroups.tsv
    "gff": {},                  # species -> GFF path                   (OrthoFinder mode)
    "max_copies": 0,            # skip orthogroups with > this many located genes in a species (0 = off)
    "tracks": [],               # [[tracks]] list: bed, hist, loop or tad lanes above/below a species
    "track_height": 0.18,       # height of one track lane (track units)
    "bedgraph": {},             # species -> .bedgraph path; colours chromosome bars by signal
    "bedgraph_cmap": "viridis",  # colormap for the bedgraph signal
    "bedgraph_vmin": None,      # signal value mapped to the low end of the colormap (default: data min)
    "bedgraph_vmax": None,      # signal value mapped to the high end (default: data max)
    "bedgraph_label": "",       # tiny legend title
    "sort": "input",            # input | name | name_reverse | size | synteny
    "min_links": 0,             # drop chromosomes with fewer coloured links (0 = keep all)
    "min_chr_size": 0,          # remove chromosomes shorter than this many bp (0 = keep all)
    "max_links": 2000,          # number of links kept from the group selected below
    "max_links_mode": "grey",  # all = cap all links; grey = cap only grey links
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

# Config
# just reading and processing config.
def read_config(path):
    with open(path, "rb") as fh:
        cfg = tomllib.load(fh)

    for key in ("chain", "focal"):
        if key not in cfg:
            raise ValueError(f"{path}: missing required key '{key}'")

    if cfg["focal"] not in cfg["chain"]:
        raise ValueError(f"{path}: focal '{cfg['focal']}' is not in chain {cfg['chain']}")

    if "orthogroups" in cfg:
        missing = []
        for species in cfg["chain"]:
            if species not in cfg.get("gff", {}):
                missing.append(species)
        if missing:
            raise ValueError(f"{path}: orthogroups mode needs a [gff] path for {missing}")
    else:
        for key in ("chr_sizes", "links"):
            if key not in cfg:
                raise ValueError(f"{path}: missing required key '{key}'")

    out = dict(DEFAULTS)
    out.update(cfg)

    if out["max_links_mode"] not in ("all", "grey"):
        raise ValueError(f"{path}: max_links_mode must be 'all' or 'grey'")

    for table in ("gff", "bedgraph", "links"):
        stray = []
        for key, value in out[table].items():
            if not isinstance(value, str):
                stray.append(key)
        if stray:
            raise ValueError(f"{path}: the [{table}] table contains non-path key(s) {stray}. Put normal settings before the TOML tables.")

    config_folder = os.path.dirname(os.path.abspath(path))

    for key in ("chr_sizes", "orthogroups", "colors_tsv", "chr_colors_tsv", "chr_names_tsv", "highlight_links_tsv"):
        value = out[key]
        if value and not os.path.isabs(value):
            out[key] = os.path.normpath(os.path.join(config_folder, value))

    for table in ("gff", "bedgraph", "links"):
        for key, value in out[table].items():
            if not os.path.isabs(value):
                out[table][key] = os.path.normpath(os.path.join(config_folder, value))

    for track in out["tracks"]:
        if "species" not in track or "file" not in track:
            raise ValueError(f"{path}: each [[tracks]] entry needs 'species' and 'file'")
        if not os.path.isabs(track["file"]):
            track["file"] = os.path.normpath(os.path.join(config_folder, track["file"]))

    out_paths = out["out"]
    if not out_paths:
        out_paths = [f"synteny_{'_'.join(out['chain'])}.png"]
    if isinstance(out_paths, str):
        out_paths = [out_paths]

    out["out"] = []
    for output_path in out_paths:
        if not os.path.isabs(output_path):
            output_path = os.path.normpath(os.path.join(config_folder, output_path))
        out["out"].append(output_path)

    return out

# Chromosome order

def natural_key(name):
    """Sort key so that 'group2' comes before 'group10'."""
    parts = []
    for text in re.split(r"(\d+)", name):
        if text.isdigit():
            parts.append(int(text))
        else:
            parts.append(text.lower())
    return parts


def _keep_only(items, wanted):
    kept = []
    for item in items:
        if item in wanted:
            kept.append(item)
    return kept


def _sort_chrs(chrs, sizes_sp, sort):
    if sort == "size":
        def chr_size(c):
            return sizes_sp[c]
        return sorted(chrs, key=chr_size, reverse=True)
    if sort == "name" or sort == "synteny" or sort == "synteny_partner":
        return sorted(chrs, key=natural_key)
    if sort == "name_reverse":
        return sorted(chrs, key=natural_key, reverse=True)
    return chrs                      # "input" and "synteny" keep the file order at this stage


def _apply_pin(chrs, pinned_list):
    pinned = _keep_only(pinned_list, chrs)
    unpinned = []
    for c in chrs:
        if c not in pinned:
            unpinned.append(c)
    return pinned + unpinned


def order_chromosomes(chr_sizes, chain, sort, chr_order, focal, focal_chrs):
    # note: sort="synteny" is handled later once the links are loaded - here it just
    # behaves like "name" for now, and gets re-sorted by sort_by_synteny afterward
    order = {}
    for sp in chain:
        chrs = list(chr_sizes.get(sp, {}).keys())
        chrs = _sort_chrs(chrs, chr_sizes[sp], sort)

        if sp in chr_order:
            chrs = _apply_pin(chrs, chr_order[sp])

        if sp == focal and focal_chrs:
            chrs = _keep_only(chrs, focal_chrs)

        order[sp] = chrs
    return order


def _links_to_upper(order_lower, upper_rank, df):
    # links_to[lower_chrom][upper_chrom] = n_links between them (upper side only counts
    # chromosomes still present in upper_rank)
    links_to = {}
    for chrom in order_lower:
        links_to[chrom] = {}
    for cu, cl in zip(df["chr1"], df["chr2"]):
        if cl in links_to and cu in upper_rank:
            links_to[cl][cu] = links_to[cl].get(cu, 0) + 1
    return links_to


def _synteny_sort_key(partners, upper_rank, mode):
    if mode == "partner":
        if partners:
            best = max(partners, key=partners.get)
            return (upper_rank[best], -partners[best])
        return (len(upper_rank), 0)

    n = sum(partners.values())                          # mode "center"
    if n:
        weighted_sum = 0
        for upper_chrom, link_count in partners.items():
            weighted_sum += upper_rank[upper_chrom] * link_count
        return (weighted_sum / n,)
    return (float(len(upper_rank)),)


def sort_by_synteny(chain, order, links_by_pair, mode):
    """Re-order every lower track to follow the track above. There is more detail in the notebook `synteny_sorting_study.ipynb`.
    mode "partner": each chromosome goes directly under the single partner it links to
        most above (most-linked grouped first) — cleanest for ~1:1 karyotypes.
    mode "center":  each chromosome goes to the average position of the chromosomes it
    links to above (barycenter) — usually fewest crossings when the correspondence is many-to-one.
    """
    new_order = {chain[0]: list(order[chain[0]])}
    for i in range(len(chain) - 1):
        upper, lower = chain[i], chain[i + 1]
        upper_rank = {}
        for rank, chrom in enumerate(new_order[upper]):
            upper_rank[chrom] = rank

        # links_to[lower_chr][upper_chr] = number of links between them
        links_to = _links_to_upper(order[lower], upper_rank, links_by_pair[(upper, lower)])

        ranked_chromosomes = []
        for old_index, chrom in enumerate(order[lower]):
            sort_key = _synteny_sort_key(links_to[chrom], upper_rank, mode)
            ranked_chromosomes.append((sort_key, old_index, chrom))

        ranked_chromosomes.sort()
        new_order[lower] = []
        for _sort_key, _old_index, chrom in ranked_chromosomes:
            new_order[lower].append(chrom)
    return new_order


# Layout

def layout(
# This is the main functions for positioning of elements defined by other functions. 
# Its directly needed for function draw()
    chain,
    order,
    chr_sizes,
    space,
    tracks_for,
    track_height,
    track_gap,
    bar_height,
):
    # Calculate the total width of each species.
    genome_widths = {}

    for species in chain:
        width = 0

        for chromosome in order[species]:
            width += chr_sizes[species][chromosome]

        genome_widths[species] = width

    widest_genome = max(genome_widths.values())
    chromosome_gap = space * widest_genome

    # add those gaps between Chr
    for species in chain:
        number_of_gaps = len(order[species]) - 1

        if number_of_gaps > 0:
            genome_widths[species] += chromosome_gap * number_of_gaps

    widest_track = max(genome_widths.values())
    scale = 1 / widest_track

    # Calculate vertical space that extra tracks need
    space_above = {}
    space_below = {}

    for species in chain:
        number_above = 0
        number_below = 0

        if species in tracks_for:
            for track in tracks_for[species]:
                if track["where"] == "above":
                    number_above += 1
                else:
                    number_below += 1

        space_above[species] = 0
        space_below[species] = 0

        if number_above > 0:
            space_above[species] = (
                number_above * (track_height + track_gap) + track_gap
            )

        if number_below > 0:
            space_below[species] = (
                number_below * (track_height + track_gap) + track_gap
            )

    # place the species from top to bottom.
    bar_y = {}
    bar_y[chain[0]] = 0

    for i in range(1, len(chain)):
        previous_species = chain[i - 1]
        species = chain[i]

        bar_y[species] = (
            bar_y[previous_species]
            - 1
            - space_below[previous_species]
            - space_above[species]
        )

    # Calculate the position of every chromosome and extra track.
    species_tracks = {}

    for species in chain:
        x_offset = (1 - genome_widths[species] * scale) / 2

        chromosome_spans = {}
        bp_before = {}
        cursor = 0

        for chromosome in order[species]:
            chromosome_size = chr_sizes[species][chromosome]

            chromosome_spans[chromosome] = (
                x_offset + cursor * scale,
                x_offset + (cursor + chromosome_size) * scale,
            )

            bp_before[chromosome] = cursor
            cursor += chromosome_size + chromosome_gap

        track_lanes = []
        number_above = 0
        number_below = 0

        if species in tracks_for:
            for track in tracks_for[species]:
                if track["where"] == "above":
                    y_center = (
                        bar_y[species]
                        + bar_height / 2
                        + track_gap
                        + track_height / 2
                        + number_above * (track_height + track_gap)
                    )
                    number_above += 1
                else:
                    y_center = (
                        bar_y[species]
                        - bar_height / 2
                        - track_gap
                        - track_height / 2
                        - number_below * (track_height + track_gap)
                    )
                    number_below += 1

                track_lanes.append((track, y_center))

        species_tracks[species] = {
            "y": bar_y[species],
            "x_offset": x_offset,
            "span": chromosome_spans,
            "bp_before": bp_before,
            "lanes": track_lanes,
        }

    y_top = (
        bar_y[chain[0]]
        + bar_height / 2
        + space_above[chain[0]]
    )

    y_bottom = (
        bar_y[chain[-1]]
        - bar_height / 2
        - space_below[chain[-1]]
    )

    return {
        "scale": scale,
        "track": species_tracks,
        "y_top": y_top,
        "y_bottom": y_bottom,
    }


# Drawing

def main_ribbon_painter(top_chr_span, bot_chr_span, y_top, y_bot, rgba):
# see `drawing_studies.ipynb` for detailed info how I come to this function. 
    xa1, xa2 = top_chr_span
    xb1, xb2 = bot_chr_span
    y_mid = (y_top + y_bot) / 2.0
    verts = [(xa1, y_top), (xa1, y_mid), (xb1, y_mid), (xb1, y_bot), (xb2, y_bot), (xb2, y_mid), (xa2, y_mid), (xa2, y_top), (xa1, y_top)]
    codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.LINETO, Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.CLOSEPOLY]
    return PathPatch(Path(verts, codes), facecolor=rgba, edgecolor="none", linewidth=0)


def _x_of(track, chrom, pos, scale):
    return track["x_offset"] + (track["bp_before"][chrom] + pos) * scale


def track_patches(trk, track, scale, in_view, lane_h, y_center):
    # BED boxes, BEDPE loop arcs, BED TAD triangles, or a bedgraph histogram.
    data, color = trk["data"], trk["color"]
    patches = []
    if trk["type"] == "bed":
        for chrom, start, end in zip(data.chrom, data.start, data.end):
            if chrom in in_view:
                x0, x1 = _x_of(track, chrom, start, scale), _x_of(track, chrom, end, scale)
                patches.append(Rectangle((x0, y_center - lane_h / 2), max(x1 - x0, 0), lane_h, facecolor=color, edgecolor="none"))
    elif trk["type"] == "loop":
        y = y_center - lane_h / 2
        a1, a2 = 0, 180
        if trk["where"] == "below":
            y = y_center + lane_h / 2
            a1, a2 = 180, 360
        for i in data.index:
            c = data.chr1[i]
            if c not in in_view:
                continue
            p1 = (data.start1[i] + data.end1[i]) / 2
            p2 = (data.start2[i] + data.end2[i]) / 2
            x1 = _x_of(track, c, p1, scale)
            x2 = _x_of(track, c, p2, scale)
            # Arc takes the whole ellipse height; we draw only half.
            patches.append(Arc(((x1 + x2) / 2, y), abs(x2 - x1), 2 * lane_h,
                               theta1=a1, theta2=a2, color=color, lw=trk["lw"], ls=trk["ls"]))
    elif trk["type"] == "tad":
        y = y_center - lane_h / 2
        tip = y + lane_h
        if trk["where"] == "below":
            y = y_center + lane_h / 2
            tip = y - lane_h
        for c, start, end in zip(data.chrom, data.start, data.end):
            if c not in in_view:
                continue
            x1 = _x_of(track, c, start, scale)
            x2 = _x_of(track, c, end, scale)
            patches.append(Polygon([(x1, y), ((x1 + x2) / 2, tip), (x2, y)],
                                   facecolor=color, edgecolor="none"))
    else:                                              # histogram from a bedgraph
        vmax = data["value"].max() or 1.0
        base = y_center - lane_h / 2
        for chrom, start, end, value in zip(data.chrom, data.start, data.end, data.value):
            if chrom in in_view:
                x0, x1 = _x_of(track, chrom, start, scale), _x_of(track, chrom, end, scale)
                patches.append(Rectangle((x0, base), max(x1 - x0, 0), max(value / vmax * lane_h, 0), facecolor=color, edgecolor="none"))
    return patches


def _foot(track, chrom, pos, half_width, scale, min_half):
    # foot is at least min_half wide (axis units) so ribbons stay visible even on huge
    # genomes where half_width bp is a tiny sliver of the track
    center = track["x_offset"] + (track["bp_before"][chrom] + pos) * scale
    half = max(half_width * scale, min_half)
    x0, x1 = track["span"][chrom]
    return max(x0, center - half), min(x1, center + half)


def bedgraph_patches(chrs, track, scale, bg_df, cmap, norm, y, bar_h):
    in_view = set(chrs)
    x_off = track["x_offset"]
    patches = []
    for chrom, start, end, value in zip(bg_df.chrom, bg_df.start, bg_df.end, bg_df.value):
        if chrom not in in_view:
            continue
        cx0, cx1 = track["span"][chrom]
        x0 = max(x_off + (track["bp_before"][chrom] + start) * scale, cx0)
        x1 = min(x_off + (track["bp_before"][chrom] + end) * scale, cx1)
        if x1 > x0:
            patches.append(Rectangle((x0, y - bar_h / 2), x1 - x0, bar_h, facecolor=cmap(norm(value)), edgecolor="none"))
    return patches


def chr_label_text(style, name, index):
    if style == "name":
        return name
    if style == "short":
        m = re.search(r"(\d+)(?:\.\d+)?$", name)
        return m.group(1) if m else name
    return str(index)        # "number"


def _figure_height(cfg, lay, n_species):
    if cfg["fig_height"]:
        return cfg["fig_height"]
    if cfg["tracks"]:
        return (lay["y_top"] - lay["y_bottom"]) * 0.9 + 0.8   # grow for track lanes
    return 0.9 * n_species + 0.8


def draw_chromosome_bars(ax, chain, order, lay, cfg, bedgraphs, norm, cmap, scale, left_edge):
    bar_h = cfg["bar_height"]
    min_w = cfg["min_bar_width"]
    chr_colors, chr_names = cfg["chr_colors"], cfg["chr_names"]
    for sp in chain:
        t = lay["track"][sp]
        y = t["y"]
        for index, c in enumerate(order[sp], start=1):
            x0, x1 = t["span"][c]
            if x1 - x0 < min_w:                   # widen just the drawn rectangle, not
                mid = (x0 + x1) / 2               # the span itself, so ribbon feet still
                x0, x1 = mid - min_w / 2, mid + min_w / 2   # attach at the true position
            ax.add_patch(Rectangle((x0, y - bar_h / 2), x1 - x0, bar_h, facecolor=chr_colors.get(c, cfg["bar_color"]), edgecolor="none"))
            if cfg["chr_label"]:
                label = chr_names.get(c) or chr_label_text(cfg["chr_label_style"], c, index)
                if cfg["chr_label_pos"] == "inside":
                    ax.text((x0 + x1) / 2, y, label, ha="center", va="center", zorder=5, fontsize=cfg["chr_label_size"], color=cfg["chr_label_color"])
                else:
                    ax.text((x0 + x1) / 2, y + bar_h / 2 + 0.05, label, ha="center", va="bottom", fontsize=cfg["chr_label_size"], color=cfg["chr_label_color"])
        if sp in bedgraphs:
            bins = bedgraph_patches(order[sp], t, scale, bedgraphs[sp], cmap, norm, y, bar_h)
            if bins:
                ax.add_collection(PatchCollection(bins, match_original=True, zorder=3))
        ax.text(left_edge - 0.015, y, cfg["labels"].get(sp, sp), ha="right", va="center", fontstyle="italic", fontsize=11)


def draw_track_lanes(ax, chain, order, lay, cfg, scale, left_edge):
    for sp in chain:
        t = lay["track"][sp]
        in_view = set(order[sp])
        for trk, yc in t["lanes"]:
            tp = track_patches(trk, t, scale, in_view, cfg["track_height"], yc)
            if tp:
                ax.add_collection(PatchCollection(tp, match_original=True, zorder=3))
            if trk["label"]:
                ax.text(left_edge - 0.015, yc, trk["label"], ha="right", va="center", fontsize=6, color=trk["color"])


def draw_ribbons(ax, chain, lay, links_by_pair, cfg, scale):
    bar_h = cfg["bar_height"]
    default_color = cfg["default_color"]
    skip_default = bool(cfg["focal_chrs"]) or cfg["drop_grey"]
    for i in range(len(chain) - 1):
        upper, lower = chain[i], chain[i + 1]
        tu, tl = lay["track"][upper], lay["track"][lower]
        y_top = tu["y"] - bar_h / 2
        y_bot = tl["y"] + bar_h / 2
        df = links_by_pair[(upper, lower)]
        patches = []
        for c1, p1, c2, p2, col in zip(df.chr1, df.pos1, df.chr2, df.pos2, df.color):
            if skip_default and col == default_color:
                continue
            if c1 not in tu["span"] or c2 not in tl["span"]:
                continue
            foot_top = _foot(tu, c1, p1, cfg["half_width"], scale, cfg["min_ribbon"])
            foot_bot = _foot(tl, c2, p2, cfg["half_width"], scale, cfg["min_ribbon"])
            grey = col == default_color and cfg["grey_alpha"] is not None
            alpha = cfg["grey_alpha"] if grey else cfg["link_alpha"]
            patches.append(main_ribbon_painter(foot_top, foot_bot, y_top, y_bot, to_rgba(col, alpha)))
        if patches:
            ax.add_collection(PatchCollection(patches, match_original=True))


def draw_highlighted_links(ax, chain, lay, highlights, cfg, scale):
    # drawn opaque and on top of the normal ribbons
    bar_h = cfg["bar_height"]
    for i in range(len(chain) - 1):
        df = highlights.get((chain[i], chain[i + 1]))
        if df is None or df.empty:
            continue
        tu, tl = lay["track"][chain[i]], lay["track"][chain[i + 1]]
        y_top, y_bot = tu["y"] - bar_h / 2, tl["y"] + bar_h / 2
        patches = []
        for c1, p1, c2, p2, hlc in zip(df.chr1, df.pos1, df.chr2, df.pos2, df.hl_color):
            if c1 in tu["span"] and c2 in tl["span"]:
                ft = _foot(tu, c1, p1, cfg["half_width"], scale, cfg["min_ribbon"] * 2)
                fb = _foot(tl, c2, p2, cfg["half_width"], scale, cfg["min_ribbon"] * 2)
                patches.append(main_ribbon_painter(ft, fb, y_top, y_bot, to_rgba(hlc, 0.9)))
        if patches:
            ax.add_collection(PatchCollection(patches, match_original=True, zorder=4))


def draw_scale_bar(ax, cfg, lay, left_edge, scale):
    # bottom-left reference bar. returns the y_bottom the axis limits should use
    # (lowered a bit to leave room for it), unchanged if scale_bar_mb is off
    y_bottom = lay["y_bottom"]
    if not cfg["scale_bar_mb"]:
        return y_bottom
    length = cfg["scale_bar_mb"] * 1_000_000 * scale
    y0 = y_bottom - 0.18
    ax.plot([left_edge, left_edge + length], [y0, y0], color="black", lw=1.5, solid_capstyle="butt")
    ax.text(left_edge + length / 2, y0 - 0.05, f"{cfg['scale_bar_mb']:g} Mb", ha="center", va="top", fontsize=7)
    return y0 - 0.05


def add_bedgraph_legend(fig, ax, bedgraphs, norm, cmap, cfg):
    if not bedgraphs:
        return
    bar = ScalarMappable(norm=norm, cmap=cmap)
    bar.set_array([])
    cax = ax.inset_axes([0.86, 0.02, 0.12, 0.03])
    cb = fig.colorbar(bar, cax=cax, orientation="horizontal")
    cb.set_ticks([norm.vmin, norm.vmax])
    cb.ax.tick_params(labelsize=5, length=2, width=0.4, pad=1)
    cb.outline.set_linewidth(0.4)
    if cfg["bedgraph_label"]:
        cb.set_label(cfg["bedgraph_label"], fontsize=5, labelpad=1)


def add_highlight_legend(ax, cfg):
    if not (cfg["highlight_legend"] and cfg["highlight_legend_entries"]):
        return
    handles = [Patch(facecolor=color, edgecolor="none", label=label) for label, color in cfg["highlight_legend_entries"]]
    ax.legend(handles=handles, loc="upper right", fontsize=7, frameon=True, framealpha=0.9, borderpad=0.5, labelspacing=0.3, handlelength=1.2, title=cfg["highlight_legend_title"] or None, title_fontsize=7)


def draw(chain, order, lay, links_by_pair, cfg, bedgraphs, norm, cmap, highlights):
    # bedgraphs is species -> DataFrame (colours the bars, adds a colorbar legend);
    # highlights is pair -> DataFrame with an 'hl_color' column, drawn opaque on top
    scale = lay["scale"]
    left_edge = min(t["x_offset"] for t in lay["track"].values())
    fig_height = _figure_height(cfg, lay, len(chain))
    fig, ax = plt.subplots(figsize=(cfg["fig_width"], fig_height))

    draw_chromosome_bars(ax, chain, order, lay, cfg, bedgraphs, norm, cmap, scale, left_edge)
    draw_track_lanes(ax, chain, order, lay, cfg, scale, left_edge)
    draw_ribbons(ax, chain, lay, links_by_pair, cfg, scale)
    draw_highlighted_links(ax, chain, lay, highlights, cfg, scale)
    y_bottom = draw_scale_bar(ax, cfg, lay, left_edge, scale)

    ax.set_xlim(left_edge - 0.22, 1.02)
    ax.set_ylim(y_bottom - 0.45, lay["y_top"] + 0.45)
    ax.axis("off")
    fig.tight_layout()

    add_bedgraph_legend(fig, ax, bedgraphs, norm, cmap, cfg)
    add_highlight_legend(ax, cfg)
    return fig


# colours 

def auto_colors(chromosomes):
    # uses matplotlib hsv colormap and gnerates hex colours for chromosomes. 
    colors = {}
    palette = matplotlib.colormaps["hsv"]
    number_of_chromosomes = len(chromosomes)

    for i, chromosome in enumerate(chromosomes):
        frac_for_colour = i / number_of_chromosomes
        colors[chromosome] = to_hex(palette(frac_for_colour))

    return colors


def setup_focal_colors(cfg, chr_sizes, focal):
    focal_chrs = list(cfg["focal_chrs"])

    if cfg["colors_tsv"]:
        focal_colors = read_colors(cfg["colors_tsv"])
    else:
        chromosomes = list(chr_sizes.get(focal, {}).keys())
        focal_colors = auto_colors(chromosomes)
        print(f"  auto-generated {len(focal_colors)} focal colours")

    if focal_chrs:
        selected_colors = {}
        for chromosome, color in focal_colors.items():
            if chromosome in focal_chrs:
                selected_colors[chromosome] = color
        focal_colors = selected_colors

    remove_grey_links = cfg["drop_grey"]
    if focal_chrs:
        remove_grey_links = True

    if not cfg["highlight_chrs"]:
        return focal_colors, focal_chrs, remove_grey_links

    highlight_chrs = list(cfg["highlight_chrs"])
    automatic_colors = auto_colors(highlight_chrs)
    highlight_colors = {}

    for chromosome in highlight_chrs:
        if chromosome in focal_colors:
            highlight_colors[chromosome] = focal_colors[chromosome]
        else:
            highlight_colors[chromosome] = automatic_colors[chromosome]

    cfg["focal_chrs"] = []
    cfg["drop_grey"] = False
    cfg["min_links"] = 0

    if cfg["grey_alpha"] is None:
        cfg["grey_alpha"] = 0.05

    return highlight_colors, [], False

## sweeping functions, pls check the 'rest_of_functions.ipynb' for more details about why they are needed.
def _sweep_right(chain, focal_i, focal, focal_colors, links_by_pair, known, default_color):
    # focal -> bottom of chain; the upper species is nearer the focal species
    for i in range(focal_i, len(chain) - 1):
        upper, lower = chain[i], chain[i + 1]
        links = links_by_pair[(upper, lower)]
        colors = []
        for chr1, gene1, chr2, gene2 in zip(links.chr1, links.gene1, links.chr2, links.gene2):
            if upper == focal:
                color = focal_colors.get(chr1, default_color)
            else:
                color = known[upper].get(str(gene1), default_color)
            gene2 = str(gene2)
            if gene2 in known[lower] and known[lower][gene2] != color:
                known[lower][gene2] = default_color
            else:
                known[lower][gene2] = color
            colors.append(color)
        links["color"] = colors

def _sweep_left(chain, focal_i, focal, focal_colors, links_by_pair, known, default_color):
    # focal -> top of chain; the lower species is nearer the focal species
    for i in range(focal_i, 0, -1):
        upper, lower = chain[i - 1], chain[i]
        links = links_by_pair[(upper, lower)]
        colors = []
        for chr1, gene1, chr2, gene2 in zip(links.chr1, links.gene1, links.chr2, links.gene2):
            if lower == focal:
                color = focal_colors.get(chr2, default_color)
            else:
                color = known[lower].get(str(gene2), default_color)
            gene1 = str(gene1)
            if gene1 in known[upper] and known[upper][gene1] != color:
                known[upper][gene1] = default_color
            else:
                known[upper][gene1] = color
            colors.append(color)
        links["color"] = colors


def assign_colors(chain, focal, focal_colors, links_by_pair, default_color):
    """Add a 'color' column to every link DataFrame, traced back to the focal species.

    Colours start on the focal chromosomes and spread further along the chain, step/chr by step/chr. The link is coloured by which of its two ends is closer the focal species.
    Chaining uses gene names, which are required in every link file.
    A link whose near end was never coloured (no ortholog back to the focal) gets default_color.
    """
    focal_i = chain.index(focal)
    known = {}      # known[sp][gene] -> colour already decided
    for sp in chain:
        known[sp] = {}

    _sweep_right(chain, focal_i, focal, focal_colors, links_by_pair, known, default_color)
    _sweep_left(chain, focal_i, focal, focal_colors, links_by_pair, known, default_color)


def filter_chromosomes(chain, focal, order, links_by_pair, default_color, min_links):
    count = {}

    for species in chain:
        count[species] = {}

    for pair, pair_links in links_by_pair.items():
        species_a, species_b = pair

        for chromosome_a, chromosome_b, color in zip(
            pair_links["chr1"],
            pair_links["chr2"],
            pair_links["color"],
        ):
            if color != default_color:
                # This coloured link adds one to both chromosomes.
                old_count = count[species_a].get(chromosome_a, 0)
                count[species_a][chromosome_a] = old_count + 1

                old_count = count[species_b].get(chromosome_b, 0)
                count[species_b][chromosome_b] = old_count + 1

    kept = {}

    for species in chain:
        if species == focal:
            kept[species] = list(order[species])
        else:
            kept[species] = []

            for chromosome in order[species]:
                number_of_links = count[species].get(chromosome, 0)

                if number_of_links >= min_links:
                    kept[species].append(chromosome)

    return kept


def sample_links(links, max_links, default_color, max_links_mode):
    if max_links_mode == "all":
        if len(links) > max_links:
            links = links.sample(n=max_links, random_state=42)
        return links

    colored = links[links["color"] != default_color]
    grey = links[links["color"] == default_color]

    if len(grey) > max_links:
        grey = grey.sample(n=max_links, random_state=42)

    return pd.concat([grey, colored])


def build_highlights(cfg, links_by_pair, order):
    # This function reads highlight_links_tsv and builds a dict of highlights and a legend.
    # Its very important one so It should clearly give an error if, lets say, requested gene ID is absent etc.
    # But keep in mind that gene IDs must be !unique enough! to identify the intended gene. Hence be carefull with miniprot outputs, bc they often comes as MP_123 etc. rename them!
    highlights = {}
    legend = []

    highlight_file = cfg["highlight_links_tsv"]
    if highlight_file is None or highlight_file == "":
        return highlights, legend

    id_color, id_label = read_highlight_links(highlight_file)
    found_ids = []
    ids_on_removed_chromosomes = []
    number_of_highlighted_links = 0

    for pair, links in links_by_pair.items():
        species_a, species_b = pair
        keep_row = []
        highlight_colors = []

        for chromosome_a, gene_a, chromosome_b, gene_b in zip(
            links["chr1"],
            links["gene1"],
            links["chr2"],
            links["gene2"],
        ):
            gene_a = str(gene_a).strip()
            gene_b = str(gene_b).strip()
            matched_ids = []

            if gene_a in id_color:
                matched_ids.append(gene_a)
            if gene_b in id_color and gene_b != gene_a:
                matched_ids.append(gene_b)

            if len(matched_ids) == 0:
                keep_row.append(False)
                continue

            for gene_id in matched_ids:
                if gene_id not in found_ids:
                    found_ids.append(gene_id)

            color = id_color[matched_ids[0]]
            for gene_id in matched_ids:
                if id_color[gene_id] != color:
                    raise ValueError(
                        f"{highlight_file}: link {gene_a} - {gene_b} "
                        f"has two different highlight colours"
                    )

            is_visible = chromosome_a in order[species_a] and chromosome_b in order[species_b]
            if is_visible:
                keep_row.append(True)
                highlight_colors.append(color)
            else:
                for gene_id in matched_ids:
                    if gene_id not in ids_on_removed_chromosomes:
                        ids_on_removed_chromosomes.append(gene_id)
                keep_row.append(False)

        highlighted = links[keep_row].copy()
        highlighted["hl_color"] = highlight_colors
        highlights[pair] = highlighted
        number_of_highlighted_links += len(highlighted)

    missing_ids = []
    for gene_id in id_color:
        if gene_id not in found_ids:
            missing_ids.append(gene_id)

    if len(missing_ids) > 0:
        raise ValueError(f"{highlight_file}: gene ID(s) not found in the links: {missing_ids}")

    if len(ids_on_removed_chromosomes) > 0:
        raise ValueError(
            f"{highlight_file}: highlighted gene ID(s) are on chromosomes "
            f"removed from the plot: {ids_on_removed_chromosomes}"
        )

    for gene_id in id_color:
        entry = (id_label[gene_id], id_color[gene_id])
        if entry not in legend:
            legend.append(entry)

    print(f"  highlighting {number_of_highlighted_links} links from {len(id_color)} ids")
    return highlights, legend


def setup_bedgraph(cfg):
    bedgraphs = {}
    all_values = []

    for species, path in cfg["bedgraph"].items():
        if species not in cfg["chain"]:
            raise ValueError(f"bedgraph species '{species}' is not in the chain")

        data = read_bedgraph(path)

        if len(data) == 0:
            raise ValueError(f"{path}: bedgraph file is empty")

        bedgraphs[species] = data

        for value in data["value"]:
            all_values.append(value)

    # Without bedgraph data there is no colour scale or colour map to make.
    if len(bedgraphs) == 0:
        norm = None
        cmap = None
        return bedgraphs, norm, cmap

    vmin = min(all_values)
    vmax = max(all_values)

    if cfg["bedgraph_vmin"] is not None:
        vmin = cfg["bedgraph_vmin"]

    if cfg["bedgraph_vmax"] is not None:
        vmax = cfg["bedgraph_vmax"]

    if vmax < vmin:
        raise ValueError("bedgraph_vmax cannot be smaller than bedgraph_vmin")

    if vmax == vmin:
        vmax = vmin + 1

    norm = Normalize(vmin=vmin, vmax=vmax)
    cmap = matplotlib.colormaps[cfg["bedgraph_cmap"]]

    return bedgraphs, norm, cmap


# This is our main plotting function, it is called from the main.py script. It takes a config file as input and produces a figure as output.
def plot(config_path):
    """Make one figure from one config file."""
    # First we read the toml config with, I decided to use tomllib package for that.
    cfg = read_config(config_path)
    # Then we put main config values into variables
    chain = cfg["chain"]
    focal = cfg["focal"]
    default_color = cfg["default_color"]

    # Load chromosome sizes and links, either from OrthoFinder or from the config's [links] table.
    ## 1. load chr colors
    cfg["chr_colors"] = {}
    if cfg["chr_colors_tsv"]: # if .tsv file is provided, we read it with read_colors() funciton.
        cfg["chr_colors"] = read_colors(cfg["chr_colors_tsv"]) # this comes from loaders.py, it reads a .tsv and returns a dict of chr + its color.
    ## 2. load chr names
    cfg["chr_names"] = {}
    if cfg["chr_names_tsv"]: # if .tsv file is provided, we read it with read_name_map() funciton. That tsv could be easily made with indexing with samtools faidx, and then editing the .fai file to have two columns: chr_name and display_name.
        cfg["chr_names"] = read_name_map(cfg["chr_names_tsv"])

    print(f"== {os.path.basename(config_path)} ==") # Prints the name of the config file being processed.
    ## 3. load orthogroups or links
    if cfg["orthogroups"]: # if orthofinder ouput is proveded and the path to the orthogroups file is provided, we load it with load_orthofinder() function.
        chr_sizes, links_by_pair = load_orthofinder(cfg, chain) # two outputs: chr_sizes and links_by_pair. chr_sizes is a dict of species -> {chr_name: chr_size}, links_by_pair is a dict of (species1, species2) -> DataFrame with columns: chr1, gene1, pos1, chr2, gene2, pos2.
    # Important note on paralogs; if a gene has many paralogs load_orthofinder() creates all possible pairs so we dont decide which is 'true' paralog. 
    # We then can exclude confusing families with many paralogs by setting max_copies in the config file. 
    # If max_copies is set to 0, we keep all paralogs. 
    # If max_copies is set to 1, we keep only single copy orthologs. 
    # If max_copies is set to 2, we keep only OGs with max 2 copies in any species, etc.

    # So if orthogroups are not provided, we switch to second mode, it was originally designed for PSYNT like outputs.
    else:
        chr_sizes = read_chr_sizes(cfg["chr_sizes"])
        links_by_pair = load_links(cfg, chr_sizes, chain)

    # Now we have all the data we need to make a figure, we can start processing it.
    # This is done by using 4 functions. 


    ## 1. we need to figure out colors for plotting, there is an focal_colors_study.ipynb file that explains how we do that.
    focal_colors, focal_chrs, remove_grey_links = setup_focal_colors(cfg, chr_sizes, focal)
    ## 2. another importatnt step is to figure out the order of chromosomes for plotting, this is done by order_chromosomes() function. It takes into account the sort mode, which can be 
    # - "input"
    # - "name"
    # -"name_reverse"
    # - "size"
    # - "synteny"
    #  "synteny_partner"
    sort_mode = cfg["sort"]
    order = order_chromosomes(
        chr_sizes,
        chain,
        sort_mode,
        cfg["chr_order"],
        focal,
        focal_chrs,
    )

    # Remove small chromosomes and scaffolds before they affect synteny sorting.
    for species in chain:
        chromosomes_to_keep = []

        for chromosome in order[species]:
            chromosome_size = chr_sizes[species][chromosome]

            if chromosome_size >= cfg["min_chr_size"]:
                chromosomes_to_keep.append(chromosome)

        order[species] = chromosomes_to_keep

    # if the sort mode is "synteny" or "synteny_partner", we need to re-order the chromosomes based on the synteny information. This is done by sort_by_synteny() function, we have a notebook explaining it.
    # difference between "synteny" and "synteny_partner" is in the way how they handle chr with multiple partners. 
    # Lets say you have 3 chr A, B, C in focal species and A1, B1, C1 in bottom one. 
    # A1 and A have 3 links shared 
    # A1 and B have 0 links shared
    # A1 and C have 2 links shared
    # In "synteny" mode A1 will be placed at approximate average between A and C, hence just under B. 
    # In "synteny_partner" mode A1 will be placed just under A, because it has the most links with A.
    # I didn't find the best way to handle this, although "synteny_partner" works much better for species maintaining macrosynteny. 
    # I therefore encourage users to try both modes and see which one works better for your data.

    if sort_mode == "synteny":
        order = sort_by_synteny(chain, order, links_by_pair, "center")

    if sort_mode == "synteny_partner":
        order = sort_by_synteny(chain, order, links_by_pair, "partner")

    # assign_colours() does a lot of work, so consult the focal_colors_study.ipynb notebook for details. 
    # In short, it assigns colors to links based on the focal species and its chromosomes.
    assign_colors(chain, focal, focal_colors, links_by_pair, default_color)

    # Count every coloured link before reducing the number of ribbons drawn.
    order = filter_chromosomes(chain, focal, order, links_by_pair, default_color, cfg["min_links"])

    sampled = {}
    for pair, pair_links in links_by_pair.items():
        if remove_grey_links:
            pair_links = pair_links[pair_links["color"] != default_color]

        sampled[pair] = sample_links(
            pair_links,
            cfg["max_links"],
            default_color,
            cfg["max_links_mode"],
        )

    highlights, cfg["highlight_legend_entries"] = build_highlights(cfg, links_by_pair, order)
    bedgraphs, norm, cmap = setup_bedgraph(cfg)

    tracks_for = load_tracks(cfg, chain)

    lay = layout(chain, order, chr_sizes, cfg["space"], tracks_for, cfg["track_height"], 0.04, cfg["bar_height"])
    fig = draw(chain, order, lay, sampled, cfg, bedgraphs, norm, cmap, highlights)

    for output_path in cfg["out"]:
        output_folder = os.path.dirname(output_path)

        if output_folder != "":
            os.makedirs(output_folder, exist_ok=True)

        if output_path.lower().endswith(".png"):
            fig.savefig(output_path, dpi=cfg["dpi"], bbox_inches="tight")
        else:
            fig.savefig(output_path, bbox_inches="tight")

        print(f"  saved {output_path}")

    plt.close(fig)
