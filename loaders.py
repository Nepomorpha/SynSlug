"""

Program     : SynSlug, synteny visualisation tool.
Language    : Python
Authors     : Yehor Tertyshnyk, Mirela Krol
Email       : egor.tertyshnyk@gmail.com
Github      : https://github.com/Nepomorpha/SynSlug
License     : MIT
Citation    : -TBS-

Description : A collection of functions for parsing/reading of input data etc

"""


import os
import pandas as pd

# READERS OF NECESSARY SUPPLEMENTARY FILES 
def read_chr_sizes(path):
    df = pd.read_csv(path, sep="\t")

    for col in ["species", "chr", "size"]:
        if col not in df.columns:
            raise ValueError(
                f"missing column '{col}' in {path}, got {list(df.columns)}. "
                f"check that the file is tab separated, "
                f"check the names of cols, should be 'species', 'chr', 'size' "
                f"I attached some recomendations how to easy generate such file "
            )

    sizes = {}
    for species, chrom, size in zip(df["species"], df["chr"], df["size"]):
        species = str(species)
        if species not in sizes:
            sizes[species] = {}
        sizes[species][str(chrom)] = int(size)

    return sizes


def read_colors(path):
    """read 'chr <tab> HEXcolour' table into {chr: color}. Header line is ignored."""
    df = pd.read_csv(path, sep="\t", header=None, dtype=str).dropna()
    df[1] = df[1].str.strip()

    if len(df) and not df[1].iloc[0].startswith("#"):
        df = df.iloc[1:]                      # first row only: could be a header row

    bad = df[~df[1].str.startswith("#")]
    if len(bad):
        raise ValueError(f"{path}: not a hex colour (missing '#'?): {bad.values.tolist()}")

    return dict(zip(df[0].astype(str), df[1].astype(str)))


def read_name_map(path):
    """Read a 2-column (real_name, display_name) TSV/CSV into {real: display}."""
    sep = "," if path.lower().endswith(".csv") else "\t"
    df = pd.read_csv(path, sep=sep, header=None, dtype=str).dropna()
    return dict(zip(df[0].astype(str), df[1].astype(str)))

# READERS OF MAIN INPUT FILES 

def read_highlight_links(path):
    """Read a highlight TSV (link id, colour, optional label) into (colour, label) dicts."""
    df = pd.read_csv(path, sep="\t", header=None, dtype=str)

    if df.shape[1] not in (2, 3):
        raise ValueError(f"{path}: highlight file needs 2 or 3 columns: gene_id, colour, optional label")

    first_color = str(df.iloc[0, 1]).strip().lower()
    if first_color in ("color", "colour", "hexcolor", "hexcolour"):
        df = df.iloc[1:].copy()

    if df.empty or df[0].isna().any() or df[1].isna().any():
        raise ValueError(f"{path}: every highlight row needs a gene ID and colour")

    df[0] = df[0].str.strip()
    df[1] = df[1].str.strip()

    if (df[0] == "").any():
        raise ValueError(f"{path}: highlight gene IDs cannot be empty")

    bad_colors = df[~df[1].str.fullmatch(r"#[0-9A-Fa-f]{6}([0-9A-Fa-f]{2})?")]
    if len(bad_colors) > 0:
        raise ValueError(f"{path}: highlight colours must be 6- or 8-digit hex colours")

    duplicate_ids = []
    seen_ids = set()
    for gene_id in df[0]:
        if gene_id in seen_ids:
            duplicate_ids.append(gene_id)
        seen_ids.add(gene_id)

    if len(duplicate_ids) > 0:
        raise ValueError(f"{path}: duplicate highlight gene ID(s): {duplicate_ids}")

    color = {}
    for gene_id, color_value in zip(df[0], df[1]):
        color[str(gene_id)] = str(color_value)

    label = {}
    if df.shape[1] >= 3:
        for gene_id, label_value in zip(df[0], df[2]):
            clean_label = str(label_value).strip()
            if pd.notna(label_value) and clean_label != "":
                label[str(gene_id)] = clean_label
            else:
                label[str(gene_id)] = str(gene_id)
    else:
        for gene_id in df[0]:
            label[str(gene_id)] = str(gene_id)
          
    return color, label


def _read_link_file(path, swap):
    df = pd.read_csv(path, sep="\t")
    ncol = len(df.columns)
    if ncol != 6:
        raise ValueError(f"{path}: link file needs 6 columns: chr1, gene1, pos1, chr2, gene2, pos2")

    df.columns = ["chr1", "gene1", "pos1", "chr2", "gene2", "pos2"]
    if df["gene1"].isna().any() or df["gene2"].isna().any():
        raise ValueError(f"{path}: gene1 and gene2 cannot be empty")

    df["gene1"] = df["gene1"].astype(str).str.strip()
    df["gene2"] = df["gene2"].astype(str).str.strip()
    if (df["gene1"] == "").any() or (df["gene2"] == "").any():
        raise ValueError(f"{path}: gene1 and gene2 cannot be empty")

    if swap:
        df = df.rename(columns={"chr1": "chr2", "gene1": "gene2", "pos1": "pos2",
                                "chr2": "chr1", "gene2": "gene1", "pos2": "pos1"})
    df = df[["chr1", "gene1", "pos1", "chr2", "gene2", "pos2"]].copy()
    df["chr1"] = df["chr1"].astype(str)
    df["chr2"] = df["chr2"].astype(str)
    df["pos1"] = df["pos1"].astype(int)
    df["pos2"] = df["pos2"].astype(int)
    return df


def read_link_path(path, a, b, chr_sizes):
    """Load one link file and orient it so side 1 = a, side 2 = b.

    Which species is in which column is detected from the chromosome names (matched
    against a's and b's chromosomes in chr_sizes), so the column order in the file does
    not matter.
    """
    df = _read_link_file(path, swap=False)
    a_chrs, b_chrs = set(chr_sizes.get(a, {})), set(chr_sizes.get(b, {}))
    side1, side2 = set(df["chr1"]), set(df["chr2"])

    as_is = len(side1 & a_chrs) + len(side2 & b_chrs)
    swapped = len(side1 & b_chrs) + len(side2 & a_chrs)
    if swapped > as_is:
        df = _read_link_file(path, swap=True)
    return df


def link_path_for_pair(links, a, b):
    """Look up the [links] entry for pair (a, b) by the key '<a>_<b>' (either order).

    The key is built from the chain species and matched exactly, so a typo just means the
    pair isn't found (a clear error in plot), never a silently wrong file.
    """
    for key in (f"{a}_{b}", f"{b}_{a}"):
        if key in links:
            return links[key]
    return None


def load_links(cfg, chr_sizes, chain):
    """Read one link table per adjacent pair in the chain (link-file mode)."""
    links_by_pair = {}
    for i in range(len(chain) - 1):
        a, b = chain[i], chain[i + 1]
        path = link_path_for_pair(cfg["links"], a, b)
        if path is None:
            raise ValueError(f'no link file for pair {a}-{b}: add  {a}_{b} = "..."  under '
                             f"[links] (have: {sorted(cfg['links'])})")
        if not os.path.exists(path):
            raise FileNotFoundError(f"the link file for pair {a}-{b} not found: {path}")
        links_by_pair[(a, b)] = read_link_path(path, a, b, chr_sizes)
        print(f"  {a}-{b}: {len(links_by_pair[(a, b)]):,} links")
      
    return links_by_pair


##################################################################
# ADDITIONAL TRACKS READERS
# For adding additional tracks and data to your chr you need to load them, those are functions helping doing that 
def _n_header_lines(path):
    """some bed/bedgraph files may start with one 'track ...' or 'browser ...' line."""
    with open(path) as fh:
        first_line = fh.readline()
    if first_line.startswith("track") or first_line.startswith("browser"):
        return 1
    return 0


def read_bedgraph(path):
    """read a 4 column .bedgraph (chrom, start, end, value)."""
    df = pd.read_csv(path, sep="\t", header=None, comment="#", dtype=str, skiprows=_n_header_lines(path), names=["chrom", "start", "end", "value"])
    df["chrom"] = df["chrom"].astype(str) 
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    df["value"] = df["value"].astype(float) 
    return df


def read_bed(path):
    """Read  first 3 columns of a BED, for costum .bed extra column is ignored, just ensure order it correct."""
    df = pd.read_csv(path, sep="\t", header=None, comment="#", dtype=str, skiprows=_n_header_lines(path))
    df = df.rename(columns={0: "chrom", 1: "start", 2: "end"})
    df["chrom"] = df["chrom"].astype(str)
    df["start"] = df["start"].astype(int) 
    df["end"] = df["end"].astype(int)
    return df[["chrom", "start", "end"]] 


def read_bedpe(path):
    """Read the first 6 BEDPE columns, one same-chromosome loop per row."""
    df = pd.read_csv(path, sep="\t", header=None, comment="#", dtype=str, skiprows=_n_header_lines(path))
    if df.shape[1] < 6:
        raise ValueError(f"{path}: loops need 6 columns: chr1, start1, end1, chr2, start2, end2")
    df = df.iloc[:, :6].copy()
    df.columns = ["chr1", "start1", "end1", "chr2", "start2", "end2"]
    if df.isna().any().any():
        raise ValueError(f"{path}: loop anchors cannot be empty")
    for col in ["start1", "end1", "start2", "end2"]:
        df[col] = df[col].astype(int)
    if (df.chr1 != df.chr2).any():
        raise ValueError(f"{path}: both loop anchors must be on the same chromosome")
    if (df.start1 < 0).any() or (df.start2 < 0).any() or (df.end1 <= df.start1).any() or (df.end2 <= df.start2).any():
        raise ValueError(f"{path}: loop intervals need 0 <= start < end")
    return df


def load_tracks(cfg, chain):
    tracks_for = {}

    for track in cfg["tracks"]:
        species = track["species"]
        path = track["file"]

        if species not in chain:
            raise ValueError(f"track species '{species}' is not in the chain")

        if "type" in track:
            track_type = track["type"]
        elif path.lower().endswith(".bedpe"):
            track_type = "loop"
        elif path.lower().endswith(".bed"):
            track_type = "bed"
        else:
            track_type = "hist"

        if track_type not in ("bed", "hist", "loop", "tad"):
            raise ValueError("track type must be 'bed', 'hist', 'loop' or 'tad'")

        if "where" in track:
            where = track["where"]
        else:
            where = "below"

        if where != "above" and where != "below":
            raise ValueError("track position must be 'above' or 'below'")

        if "color" in track:
            color = track["color"]
        else:
            color = "#555555"

        if "label" in track:
            label = track["label"]
        else:
            label = ""

        if track_type == "bed" or track_type == "tad":
            data = read_bed(path)
        elif track_type == "loop":
            data = read_bedpe(path)
        else:
            data = read_bedgraph(path)

        if species not in tracks_for:
            tracks_for[species] = []

        tracks_for[species].append({
            "type": track_type,
            "color": color,
            "label": label,
            "where": where,
            "data": data,
            "lw": track.get("lw", 1.0),
            "ls": track.get("ls", "-"),
        })

    return tracks_for
