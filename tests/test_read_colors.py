import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import loaders

DATA = os.path.join(os.path.dirname(__file__), "data")
try:
    print(loaders.read_colors(os.path.join(DATA, "rc1.tsv")))
except ValueError as e:
    print(e)
print(loaders.read_colors(os.path.join(DATA, "rc2.tsv")))
