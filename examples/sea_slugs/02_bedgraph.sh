#!/bin/bash
set -e
set -o pipefail

# from repo root
mkdir -p example_data/sea_slugs/work

# chr sizes from the header
gzip -dc example_data/sea_slugs/downloads/Berghia_stephanieae.gff.gz |
    awk 'BEGIN {OFS="\t"} $1 == "##sequence-region" && $3 == 1 && $4 >= 1000000 {print $2, $4}' > example_data/sea_slugs/work/sizes.tsv

# gene centres
gzip -dc example_data/sea_slugs/downloads/Berghia_stephanieae.gff.gz |
    awk 'BEGIN {OFS="\t"} $3 == "gene" {p = int(($4 + $5) / 2) - 1; print $1, p, p + 1}' > example_data/sea_slugs/work/genes.bed

bedtools makewindows -g example_data/sea_slugs/work/sizes.tsv -w 1000000 > example_data/sea_slugs/work/windows.bed
bedtools coverage -a example_data/sea_slugs/work/windows.bed -b example_data/sea_slugs/work/genes.bed -counts > example_data/sea_slugs/work/counts.tsv

# genes / Mb
awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4 * 1000000 / ($3 - $2)}' example_data/sea_slugs/work/counts.tsv > example_data/sea_slugs/gene_density.bedgraph
