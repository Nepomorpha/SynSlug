import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import loaders

DATA = os.path.join(os.path.dirname(__file__), "data")
print(loaders.read_chr_sizes(os.path.join(DATA, "cs1.tsv")))
