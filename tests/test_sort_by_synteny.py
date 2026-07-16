import os, sys
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import functions

DATA = os.path.join(os.path.dirname(__file__), "data")


def load(name):
    return pd.read_csv(os.path.join(DATA, name), sep="\t")


order = {"A": ["A1", "A2"], "B": ["B1", "B2", "B_ORPHAN"]}
links = load("sb1.tsv")
print(functions.sort_by_synteny(["A", "B"], order, {("A", "B"): links}, "partner")["B"])
print(functions.sort_by_synteny(["A", "B"], order, {("A", "B"): links}, "center")["B"])
print()

order = {"A": ["A1", "A2"], "B": ["B3", "B1", "B2"]}
links = load("sb2.tsv")
print(functions.sort_by_synteny(["A", "B"], order, {("A", "B"): links}, "partner")["B"])
print(functions.sort_by_synteny(["A", "B"], order, {("A", "B"): links}, "center")["B"])
print()

order = {"A": ["A1", "A2", "A3", "A4"], "B": ["B1", "B2", "B3", "B4", "B_TIE"]}
links = load("sb3.tsv")
print(functions.sort_by_synteny(["A", "B"], order, {("A", "B"): links}, "partner")["B"])
print(functions.sort_by_synteny(["A", "B"], order, {("A", "B"): links}, "center")["B"])
