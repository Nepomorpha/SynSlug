import os, sys
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import functions

DATA = os.path.join(os.path.dirname(__file__), "data")
DEFAULT = "#e3e3e3"


def load(name):
    return pd.read_csv(os.path.join(DATA, name), sep="\t")


def show(links):
    for (a, b), df in links.items():
        print(a, b)
        print(df[["chr1", "gene1", "chr2", "gene2", "color"]])


chain = ["X", "Y", "Z", "W"]
links = {("X", "Y"): load("ac1a.tsv"), ("Y", "Z"): load("ac1b.tsv"), ("Z", "W"): load("ac1c.tsv")}
functions.assign_colors(chain, "Y", {"Y1": "#ff0000", "Y2": "#0000ff"}, links, DEFAULT)
show(links)
print()

chain = ["X", "Y", "Z"]
links = {("X", "Y"): load("ac2a.tsv"), ("Y", "Z"): load("ac2b.tsv")}
functions.assign_colors(chain, "X", {"X1": "#ff0000", "X2": "#00ff00"}, links, DEFAULT)
show(links)
print()

chain = ["X", "Y", "Z"]
links = {("X", "Y"): load("ac3a.tsv"), ("Y", "Z"): load("ac3b.tsv")}
functions.assign_colors(chain, "X", {"X1": "#ff0000"}, links, DEFAULT)
show(links)
