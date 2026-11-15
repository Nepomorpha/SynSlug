"""Builds synteny links from OrthoFinder output (Orthogroups.tsv + per-species GFFs)"""

import gzip
import pandas as pd
import pysam
from loaders import read_chr_sizes


def _check_species_columns(gene_og, chain, orthogroups_path):
    for sp in chain:
        if sp not in gene_og:
            raise ValueError(f"species '{sp}' is not a column in {orthogroups_path}; columns are {sorted(gene_og)}")


def _resolve_chr_sizes(cfg, chain, gff_sizes):
    sizes = dict(gff_sizes)
    for sp, path in cfg.get("genomes", {}).items():
        sizes[sp] = {}
        with pysam.FastxFile(path) as fasta:
            for seq in fasta:
                sizes[sp][seq.name] = len(seq.sequence)
    if cfg.get("chr_sizes"):
        sizes.update(read_chr_sizes(cfg["chr_sizes"]))

    no_size = []
    for sp in chain:
        if not sizes.get(sp):
            no_size.append(sp)
    if no_size:
        raise ValueError(f"no chromosome lengths for {no_size}: add their genome FASTA paths under [genomes], or provide a chr_sizes TSV")
    return sizes

def load_orthofinder(cfg, chain):
    """Build chr_sizes + per pair link tables from orthofinder output and GFFs"""
    print("  reading orthogroups ...")
    gene_og = read_orthogroups(cfg["orthogroups"])
    _check_species_columns(gene_og, chain, cfg["orthogroups"])

    positions, gff_sizes = {}, {}
    for sp in chain:
        positions[sp], gff_sizes[sp] = read_gff_positions(cfg["gff"][sp])
        n = 0
        for gene in gene_og[sp]:
            if gene in positions[sp]:
                n += 1
        print(f"  {sp}: {n:,} of {len(gene_og[sp]):,} orthogroup proteins matched to GFF")
        if n == 0:
            raise ValueError(f"{sp}: no OrthoFinder protein IDs match the GFF. Check that the proteins and annotation are from the same gene set.")

    links_by_pair = build_orthofinder_links(chain, gene_og, positions, cfg["max_copies"])
    for (a, b), df in links_by_pair.items():
        print(f"  {a}-{b}: {len(df):,} gene-pair links")
        if df.empty:
            raise ValueError(f"no gene-pair links for {a}-{b}: check the protein IDs and max_copies setting")

    chr_sizes = _resolve_chr_sizes(cfg, chain, gff_sizes)
    return chr_sizes, links_by_pair 


def read_orthogroups(path):
    """Orthogroups.tsv -> {species: {gene_id: orthogroup}}."""
    df = pd.read_csv(path, sep="\t", dtype=str)
    df = df.fillna("")
    species = [c for c in df.columns if c != "Orthogroup"]     # rest of the columns are species
    gene_og = {sp: {} for sp in species} 
    for row in df.itertuples(index=False):
        og = row[0] 
        for sp, cell in zip(species, row[1:]):
            if cell:                                            # cell is a ", "-separated gene list 
                for gene in cell.split(", "):
                    gene_og[sp][gene] = og
    return gene_og


def read_gff_positions(path):
    """gff -> ids, positions, chr sizes"""
    path = str(path)
    cds = {}
    versions = {}
    sizes = {}

    if path.endswith((".gz", ".bgz")):
        fh = gzip.open(path, "rt")
    else:
        fh = open(path)

    with fh:
        for line in fh:
            if line.startswith("##sequence-region"):
                row = line.split()
                sizes[row[1]] = int(row[3])
                continue
            if line.startswith("##FASTA"):
                break
            if line.startswith("#"):
                continue

            row = line.rstrip().split("\t")
            if len(row) != 9:
                continue
            if row[2].upper() != "CDS":
                continue

            chrom = row[0]
            start = int(row[3])
            end = int(row[4])

            # last col, protein_id=P1;version=1
            attrs = {}
            for entry in row[8].split(";"):
                if "=" in entry:
                    key, value = entry.split("=", 1)
                    attrs[key.strip()] = value.strip()

            if attrs.get("Target"):
                if "rank" in attrs and attrs["rank"] != "1":
                    continue
                gene = attrs["Target"].split()[0]
            elif attrs.get("protein_id"):
                gene = attrs["protein_id"]
                if attrs.get("version"):
                    version = attrs["version"]
                    if not gene.endswith("." + version):
                        versions[gene + "." + version] = gene
            elif attrs.get("Parent"):
                gene = attrs["Parent"]
            else:
                continue

            # CDS span
            if gene not in cds:
                cds[gene] = [chrom, start, end]
            else:
                cds[gene][1] = min(cds[gene][1], start)
                cds[gene][2] = max(cds[gene][2], end)

    positions = {}
    for gene, span in cds.items():
        chrom, start, end = span
        positions[gene] = (chrom, (start + end) // 2)

    # P1 and P1.1
    for name, gene in versions.items():
        positions[name] = positions[gene]

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
