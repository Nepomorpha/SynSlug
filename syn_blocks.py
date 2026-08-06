#!/usr/bin/env python3
import argparse
import pandas as pd
from scipy.spatial.distance import pdist
from sklearn.cluster import AgglomerativeClustering
 
parser = argparse.ArgumentParser()
parser.add_argument("--mlt", required=True)
parser.add_argument("--output", required=True)
parser.add_argument("--maximum-mlt-distance", type=float, default=0.1)
parser.add_argument("--minimum-anchors", type=int, default=3)
parser.add_argument("--minimum-span", type=int, default=20_000)
parser.add_argument("--maximum-span", type=int, default=3_000_000)
parser.add_argument("--maximum-reference-gap", type=int, default=500_000)
args = parser.parse_args()
