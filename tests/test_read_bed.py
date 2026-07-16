import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import loaders

DATA = os.path.join(os.path.dirname(__file__), "data")
print(loaders.read_bed(os.path.join(DATA, "rb1.bed")))
try:
    print(loaders.read_bed(os.path.join(DATA, "rb2.bed")))
except Exception as e:
    print(type(e).__name__, e)
