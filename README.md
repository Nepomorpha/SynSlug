![SynSlug](synslug_baner.png.webp)

**Step 1. Clone SynSlug**

Let's get SynSlug onto your machine! 
Navigate to your project folder and run the following commands to clone the repository:

```bash
mkdir -p software
git clone https://github.com/Nepomorpha/SynSlug.git software/SynSlug
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r software/SynSlug/requirements.txt
```



**Step 2. Get your genomes and proteomes**

Let's say we want to visualise synteny across five sea slug species. 
Where do we start? 
First, we need to collect the genome and proteome files for each species.

In this guide, we'll use data from NCBI. To make things easier, let's walk through the process using Berghia stephanieae as our first example.

We'll download both the genome assembly and its corresponding proteome, which we'll need for the downstream synteny analysis.

For example, let's download the required files for Berghia stephanieae. We'll need three files:

Protein sequences (.faa.gz) – for the OrthoFinder analysis.
Gene annotations (.gff.gz) – to determine gene positions on chromosomes.
Genome sequences (.fa.gz) – to retrieve chromosome or scaffold lengths when needed.

You can download all three directly from NCBI using curl:

```bash
curl -fL --retry 3 https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/034/508/935/GCF_034508935.2_UCSD_Bste_1.3/GCF_034508935.2_UCSD_Bste_1.3_protein.faa.gz -o "$work/downloads/Berghia_stephanieae.faa.gz"
curl -fL --retry 3 https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/034/508/935/GCF_034508935.2_UCSD_Bste_1.3/GCF_034508935.2_UCSD_Bste_1.3_genomic.gff.gz -o "$work/downloads/Berghia_stephanieae.gff.gz"
curl -fL --retry 3 https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/034/508/935/GCF_034508935.2_UCSD_Bste_1.3/GCF_034508935.2_UCSD_Bste_1.3_genomic.fna.gz -o "$work/downloads/Berghia_stephanieae.fa.gz"
```

You can also compress the files using gzip if needed:

```bash
gzip -dc "$WORK/downloads/Berghia_stephanieae.faa.gz" > "$WORK/proteins/Berghia_stephanieae.fa"
gzip -dc "$WORK/downloads/Phyllidia_flava.faa.gz" > "$WORK/proteins/Phyllidia_flava.fa"
...
```



**Step 3. Run OrthoFinder**

Next, it's time to run OrthoFinder (ref). You can run it in whichever way works best for you.

Although Guide.me supports other input formats, such as links to files, OrthoFinder is probably the easiest option to start with, as it requires the least amount of file preparation.


```bash
OrthoFinder-2.5.5/orthofinder.py \
    -f /proteins \
    -o /results \
    -S diamond -t 12 -a 4 -og
```

For now, the file we're most interested in is `/results/Results_Nov15/Orthogroups/Orthogroups.tsv`

This file contains 19,312 orthogroups (plus one header line). We'll use it as the starting point for the next steps of our synteny visualisation workflow.

```bash
wc -l /results/Results_Nov15/Orthogroups/Orthogroups.tsv
19313 /results/Results_Nov15/Orthogroups/Orthogroups.tsv 
   ```



**Step 4. Create and configure config.toml**

Now that we have our orthogroups and annotation files, it's time to give SynSlug the instructions how to use them. We do this by creating a config.toml file.

A few key settings control how the synteny plot will look:

`chain` defines the order in which species will appear in the plot.

`focal` selects the reference species. Its chromosomes determine the colours used in the plot.

`[gff]` specifies the GFF annotation file for each species. Chromosome positions are inferred from these files.



For our five sea slugs, the complete configuration looks like this:


```toml
orthogroups = "results/Results_Nov15/Orthogroups/Orthogroups.tsv"
chain = ["Berghia_stephanieae", "Phyllidia_flava", "Elysia_timida", "Elysia_crispata", "Onchidella_celtica"]
focal = "Berghia_stephanieae"

out = "results/sea_slugs.png"

min_chr_size = 1_000_000
max_links = 2000
max_links_mode = "all"
dpi = 150

[gff]
Berghia_stephanieae = "downloads/Berghia_stephanieae.gff.gz"
Phyllidia_flava = "downloads/Phyllidia_flava.gff.gz"
Elysia_timida = "downloads/Elysia_timida.gff.gz"
Elysia_crispata = "downloads/Elysia_crispata.gff.gz"
Onchidella_celtica = "downloads/Onchidella_celtica.gff.gz"

[genomes]
Elysia_timida = "downloads/Elysia_timida.fa.gz"
```
As you can see you don't need to decompress your files first – SynSlug can read the original compressed GFFs and genome files directly.

A couple of basic settings are worth highlighting:

`min_chr_size` = 1_000_000 skips scaffolds shorter than 1 Mb, helping to keep the plot manageable.

`max_links` = 2000 limits the number of links drawn for each species pair, so the first plot renders faster. You can increase this value if you want to see more connections.

`max_links_mode` = "all" applies this limit across all links.

`dpi` = 150 sets the resolution of the output image.


Once your configuration file is ready, run SynSlug from the project folder:

```bash
source .venv/bin/activate
python software/SynViser/synviser.py config.toml
```

This saves `results/sea_slugs.png`.


**Want to customise your plot?**

What we've covered so far is just the basic functionality of SynSlug! The tool is highly customisable, allowing you to tweak many aspects of your visualisation, including colours, link opacity, gene highlighting and even loop visualisation.

There's much more to explore! You'll find instructions for all these additional features in the following sections of the guide.
