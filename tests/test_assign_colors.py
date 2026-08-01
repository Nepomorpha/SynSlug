import os, sys
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import functions

DATA = os.path.join(os.path.dirname(__file__), "data")
DEFAULT = "#e3e3e3"


def load(name):
    return pd.read_csv(os.path.join(DATA, name), sep="\t")


chain = ["X", "Y", "Z", "W"]
links = {("X", "Y"): load("ac1a.tsv"), ("Y", "Z"): load("ac1b.tsv"), ("Z", "W"): load("ac1c.tsv")}
functions.assign_colors(chain, "Y", {"Y1": "#ff0000", "Y2": "#0000ff"}, links, DEFAULT)
print("X-Y")
print(links[("X", "Y")])
print("Y-Z")
print(links[("Y", "Z")])
print("Z-W")
print(links[("Z", "W")])

chain = ["X", "Y", "Z"]
links = {("X", "Y"): load("ac3a.tsv"), ("Y", "Z"): load("ac3b.tsv")}
functions.assign_colors(chain, "X", {"X1": "#ff0000"}, links, DEFAULT)
print("X-Y")
print(links[("X", "Y")])
print("Y-Z")
print(links[("Y", "Z")])
