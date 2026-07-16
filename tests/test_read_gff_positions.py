import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import orthofinder as of

DATA = os.path.join(os.path.dirname(__file__), "data")
print(of.read_gff_positions(os.path.join(DATA, "gp1.gff")))
print(of.read_gff_positions(os.path.join(DATA, "gp2.gff")))
