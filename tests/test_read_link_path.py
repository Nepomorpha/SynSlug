import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import loaders

DATA = os.path.join(os.path.dirname(__file__), "data")

df = loaders.read_link_path(os.path.join(DATA, "lp1.tsv"), "A", "B", {"A": {"a1": 1000}, "B": {"b1": 1000}})
print(set(df["chr1"]))

df = loaders.read_link_path(os.path.join(DATA, "lp2.tsv"), "A", "B", {"A": {"a1": 1000}, "B": {"b1": 1000}})
print(set(df["chr1"]))

chr_sizes = {"A": {"1": 1000, "2": 2000}, "B": {"1": 1000, "4": 4000}}
df = loaders.read_link_path(os.path.join(DATA, "lp3.tsv"), "A", "B", chr_sizes)
print(set(df["chr1"]))
