#!/bin/bash
set -e

# from repo root
mkdir -p example_data/sea_slugs/downloads

curl -fL --retry 3 https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/034/508/935/GCF_034508935.2_UCSD_Bste_1.3/GCF_034508935.2_UCSD_Bste_1.3_genomic.gff.gz -o example_data/sea_slugs/downloads/Berghia_stephanieae.gff.gz
curl -fL --retry 3 https://ftp.ebi.ac.uk/pub/ensemblorganisms/GCA/964/266/695/1/ensembl/2025_03/geneset/genes.gff3.gz -o example_data/sea_slugs/downloads/Phyllidia_flava.gff.gz
curl -fL --retry 3 https://genome.senckenberg.de/download/etim/braker.noc.gff3.gz -o example_data/sea_slugs/downloads/Elysia_timida.gff.gz
curl -fL --retry 3 https://ftp.ebi.ac.uk/pub/ensemblorganisms/GCA/963/854/125/1/ensembl/2024_05/geneset/genes.gff3.gz -o example_data/sea_slugs/downloads/Elysia_crispata.gff.gz
curl -fL --retry 3 https://ftp.ebi.ac.uk/pub/ensemblorganisms/GCA/963/931/925/1/ensembl/2024_05/geneset/genes.gff3.gz -o example_data/sea_slugs/downloads/Onchidella_celtica.gff.gz

# true sizes for E. timida
curl -fL --retry 3 https://genome.senckenberg.de/download/etim/ely-polyC.yahs_scaffolds_final.FINAL.sort.nomt.noc.fa.gz -o example_data/sea_slugs/downloads/Elysia_timida.fa.gz
