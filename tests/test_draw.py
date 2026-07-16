import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import functions

DATA = os.path.join(os.path.dirname(__file__), "data")

for name in ["dr1.toml", "dr2.toml", "dr3.toml", "dr4.toml"]:
    functions.plot(os.path.join(DATA, name))
    print("saved", name.replace(".toml", ".png"))
