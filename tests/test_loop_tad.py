from pathlib import Path
from tempfile import TemporaryDirectory
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from functions import read_config, layout, track_patches
from loaders import load_tracks, read_bedpe


ROOT = Path(__file__).resolve().parents[1]
cfg = read_config(ROOT / "examples/loops_tads.toml")
tracks = load_tracks(cfg, cfg["chain"])
loops, tads = tracks["spA"]
assert loops["data"].chr1.tolist() == ["A1", "A2"]
assert tads["type"] == "tad"

# A2 comes after A1, so its loop must include that offset.
sizes = {"spA": {"A1": 20_000_000, "A2": 18_000_000}, "spB": {"B1": 20_000_000}}
order = {"spA": ["A1", "A2"], "spB": ["B1"]}
lay = layout(["spA", "spB"], order, sizes, 0.01, tracks, 0.18, 0.04, 0.05)
t = lay["track"]["spA"]
scale = lay["scale"]

for trk, yc in t["lanes"]:
    patches = track_patches(trk, t, scale, set(order["spA"]), 0.18, yc)
    assert len(patches) == 2
    box = patches[1].get_path().get_extents(patches[1].get_transform())
    left = t["span"]["A2"][0] + 2_200_000 * scale
    right = t["span"]["A2"][0] + 12_800_000 * scale
    assert np.allclose([box.xmin, box.xmax], [left, right])
    assert np.allclose([box.ymin, box.ymax], [yc - 0.09, yc + 0.09])
    hidden = track_patches(trk, t, scale, {"A1"}, 0.18, yc)
    assert len(hidden) == 1

# Also try a loop below and a TAD above.
loops["where"] = "below"
tads["where"] = "above"
loops["ls"] = "--"
for trk in [loops, tads]:
    p = track_patches(trk, t, scale, {"A1"}, 0.18, 0)[0]
    box = p.get_path().get_extents(p.get_transform())
    assert np.allclose([box.ymin, box.ymax], [-0.09, 0.09])
assert p.get_xy()[1, 1] == 0.09  # triangle tip above
p = track_patches(loops, t, scale, {"A1"}, 0.18, 0)[0]
assert p.theta1 == 180 and p.theta2 == 360
assert p.get_linestyle() == "--"

with TemporaryDirectory() as folder:
    path = Path(folder) / "loops.bedpe"
    # Extra BEDPE columns and reversed anchor order are allowed.
    path.write_text("track name=loops\n1\t90\t100\t1\t10\t20\tloop1\n")
    df = read_bedpe(path)
    assert df.chr1[0] == "1"
    loops["data"] = df
    t1 = {"x_offset": 0, "bp_before": {"1": 0}}
    p = track_patches(loops, t1, 0.01, {"1"}, 0.18, 0)[0]
    assert np.isclose(p.width, 0.8)

    for row in ["A1\t10\t20\tA2\t30\t40\n",
                "A1\t-1\t20\tA1\t30\t40\n",
                "A1\t20\t10\tA1\t30\t40\n",
                "A1\t10\t20\tA1\t30\n"]:
        path.write_text(row)
        try:
            read_bedpe(path)
        except ValueError:
            pass
        else:
            raise AssertionError("invalid loop was accepted")

print("Loop and TAD checks passed")
