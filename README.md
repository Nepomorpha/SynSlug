# SynViser
A simple matplotlib based program for easy and adjustable plotting of Ribbon-plots to explore synteny. 

We suggest you to try our [Jupyter Notebook](01_genome_analysis_pipeline.ipynb) for step by step analysis of your dataset.

Example configs and plots are in [examples](examples/).
The smaller notebooks are in [studies](studies/), and extra scripts are in [scripts](scripts/).

Jump to relevant section:
- (getting started)
  - (cloning and env)   
  - (preparing inputs) 
  - (setting up config)
- (running whole genome and single chr)
- modifying appearance of plots 
- additional analysis 

## Installation

SynViser requires Python 3.11 or newer.

```
git clone https://github.com/Nepomorpha/SynViser.git
cd SynViser
pip install -r requirements.txt
```

## Quick start
Run one of examples:

```
python synviser.py examples/my_config.toml
```

This creates `examples/my_plot.png`.

To start with your data you should either have an Orhofinder output or a tab separated files of gene pairs in exactly this order::

```text
chr1  gene1  pos1  chr2  gene2  pos2
```

Every row needs a non empty gene identifier on both sides. For an intermediate
species, use the same gene IDs in both adjacent link files so colours can continue
through the species chain.

If you happen to have 4-column file you can convert it with:

```bash
python scripts/add_gene_names.py < old_links.tsv > links.tsv
```

The generated ID is `chromosome:position`. Use real gene IDs when available.
