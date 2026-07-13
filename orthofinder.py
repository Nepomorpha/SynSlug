"""Builds synteny links from OrthoFinder output (Orthogroups.tsv + per-species GFFs)"""

import pandas as pd
from loaders import read_chr_sizes


def _check_species_columns(gene_og, chain, orthogroups_path):
    for sp in chain:
        if sp not in gene_og:
            raise ValueError(f"species '{sp}' is not a column in {orthogroups_path}; columns are {sorted(gene_og)}")


def _resolve_chr_sizes(cfg, chain, gff_sizes):
    # real GFF sizes
    if cfg["chr_sizes"]:
        return read_chr_sizes(cfg["chr_sizes"]) 

    no_size = [sp for sp in chain if not gff_sizes[sp]]
    if no_size:
        raise ValueError(f"no chr_sizes for {no_size}: here its either their GFF has no '##sequence-region' lines (typical for miniprot output), so real chromosome sizes are not in file. Add chr_sizes = \"...\" to the config, for these species.")
    return gff_sizes

def load_orthofinder(cfg, chain):
    """Build chr_sizes + per pair link tables from orthofinder output and GFFs"""
    print("  reading orthogroups ...")
    gene_og = read_orthogroups(cfg["orthogroups"])
    _check_species_columns(gene_og, chain, cfg["orthogroups"])

    positions, gff_sizes = {}, {}
    for sp in chain:
        positions[sp], gff_sizes[sp] = read_gff_positions(cfg["gff"][sp])
        print(f"  {sp}: {len(positions[sp]):,} proteins located from gff")

    links_by_pair = build_orthofinder_links(chain, gene_og, positions, cfg["max_copies"])
    for (a, b), df in links_by_pair.items():
        print(f"  {a}-{b}: {len(df):,} gene-pair links")

    chr_sizes = _resolve_chr_sizes(cfg, chain, gff_sizes)
    return chr_sizes, links_by_pair 


def read_orthogroups(path):
    """Orthogroups.tsv -> {species: {gene_id: orthogroup}}."""
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    species = [c for c in df.columns if c != "Orthogroup"]     # rest of the columns are species
    gene_og = {sp: {} for sp in species} 
    for row in df.itertuples(index=False):
        og = row[0] 
        for sp, cell in zip(species, row[1:]):
            if cell:                                            # cell is a ", "-separated gene list 
                for gene in cell.split(", "):
                    gene_og[sp][gene] = og
    return gene_og


def _gene_positions_from_cds(cds):
    grp = cds.groupby("gene")
    chrom = grp["seqid"].first()
    midpoint = (grp["start"].min() + grp["end"].max()) // 2
    positions = {}
    for gene in chrom.index:
        positions[gene] = (chrom[gene], int(midpoint[gene]))
    return positions


def _sequence_region_sizes(path):
    # NCBI has ##sequence-region lines through the file (one before each sequence's
    # features), so we just scan the whole thing rather than only the header block
    sizes = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("##sequence-region"):
                _, name, _start, end = line.split()
                sizes[name] = int(end)
    return sizes


def read_gff_positions(path):
    """GFF (NCBI or miniprot) -> ({gene_id: (chr, midpoint)}, {chr: size}).Position is the midpoint of the CDS span.

    Sizes come ONLY from ##sequence-region lines. note: miniprot GFFs don't have them -> diff handling is needed. (applied?)
    """
    cols = ["seqid", "source", "type", "start", "end", "attr"]
    gff = pd.read_csv(path, sep="\t", comment="#", header=None, usecols=[0, 1, 2, 3, 4, 8], names=cols, dtype=str, on_bad_lines="skip")        # tolerate odd/tral lines
    gff = gff.dropna(subset=["seqid", "type", "start", "end"])
    gff["start"] = gff["start"].astype(int)
    gff["end"] = gff["end"].astype(int)

    miniprot = (gff["source"] == "miniprot").any()
    cds = gff[gff["type"] == "CDS"].copy()
    cds["gene"] = cds["attr"].str.extract(r"Parent=([^;]+)" if miniprot else r"protein_id=([^;]+)")
    cds = cds.dropna(subset=["gene"]) 

    positions = _gene_positions_from_cds(cds)
    sizes = _sequence_region_sizes(path)
    return positions, sizes


def _cross_join_orthogroup(og_genes_a, og_genes_b):
    rows = []
    for gene_a, chr_a, pos_a in og_genes_a: 
        for gene_b, chr_b, pos_b in og_genes_b:
            rows.append([chr_a, gene_a, pos_a, chr_b, gene_b, pos_b])
    return rows


def _pair_links(genes_by_og_a, genes_by_og_b, max_copies):
    rows = []
    for og in genes_by_og_a:
        if og not in genes_by_og_b:
            continue
        genes_a = genes_by_og_a[og] 
        genes_b = genes_by_og_b[og] 
        if max_copies and (len(genes_a) > max_copies or len(genes_b) > max_copies):
            continue 
        rows.extend(_cross_join_orthogroup(genes_a, genes_b))
    return rows


def build_orthofinder_links(chain, gene_og, positions, max_copies):
    # every gene in species  A gets linked to gene in species B that shares an orthogroup.
    genes_by_og = {} 
    for sp in chain:
        genes_by_og[sp] = genes_by_orthogroup(gene_og[sp], positions[sp])

    links_by_pair = {}
    for i in range(len(chain) - 1):
        a, b = chain[i], chain[i + 1]
        rows = _pair_links(genes_by_og[a], genes_by_og[b], max_copies)
        links_by_pair[(a, b)] = pd.DataFrame(rows, columns=["chr1", "gene1", "pos1", "chr2", "gene2", "pos2"]) 
    return links_by_pair


def genes_by_orthogroup(gene_og, positions):
    by_og = {}
    for gene, og in gene_og.items():
        if gene in positions:
            if og not in by_og:
                by_og[og] = []
            chrom, pos = positions[gene]
            by_og[og].append((gene, chrom, pos))
    return by_og
