"""
A synthetic input for processor-cytometry-umap, shaped like processor-flowsom's
`cells` output and fcs-prep's `metadata`:

  <root>/linked/cells/cells.csv.gz     sample_id, cluster, metacluster, label, CD3, CD4, ...
  <root>/linked/cells/cells.json       {"clustering_markers": [...], ...}
  <root>/linked/metadata/samples.csv   sample_id, file, batch, condition, n_events
  <root>/flat/...                      the same files side by side, nested one level

Four populations (CD4 T, CD8 T, B cells, NK) in four samples, two batches,
stimulated samples with higher CD69 on T cells; arcsinh-like values.

  python -m tests.make_fixture data/fixture [--cells 300]
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

CLUSTERING = ["CD3", "CD4", "CD8", "CD19", "CD56", "CD16"]
POPS = {
    "CD4 T": {"CD3": 2.5, "CD4": 2.5},
    "CD8 T": {"CD3": 2.5, "CD8": 3.0},
    "B cells": {"CD19": 2.5},
    "NK CD56dim CD16+": {"CD56": 2.5, "CD16": 2.0},
}
SAMPLES = [("b1_unstim", "batch01", "unstimulated"), ("b1_stim", "batch01", "stimulated"),
           ("b2_unstim", "batch02", "unstimulated"), ("b2_stim", "batch02", "stimulated")]


def cells_table(n_per_sample: int, seed: int = 1) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    rows = []
    for sid, batch, cond in SAMPLES:
        for k, (label, mu) in enumerate(POPS.items(), start=1):
            n = n_per_sample // len(POPS)
            d = {"sample_id": sid, "cluster": k * 3, "metacluster": k, "label": label}
            df = pd.DataFrame(d, index=range(n))
            for m in CLUSTERING:
                df[m] = rng.normal(mu.get(m, 0.0), 0.25, n) + (0.3 if batch == "batch02" else 0.0)
            df["CD69"] = rng.normal(2.0 if cond == "stimulated" and "T" in label else 0.0, 0.3, n)
            rows.append(df)
    return pd.concat(rows, ignore_index=True).round(4)


def make_fixture(root: Path, n_per_sample: int = 300) -> dict[str, Path]:
    cells = cells_table(n_per_sample)
    meta = {"clustering_markers": CLUSTERING, "state_markers": ["CD69"], "markers": CLUSTERING + ["CD69"],
            "cells_per_sample": n_per_sample, "seed": 42}
    sheet = pd.DataFrame([{"sample_id": s, "file": f"{s}.fcs", "batch": b, "condition": c, "n_events": n_per_sample}
                          for s, b, c in SAMPLES])
    linked = root / "linked"
    flat = root / "flat" / "nested"
    for cells_dir, meta_dir in ((linked / "cells", linked / "metadata"), (flat, flat)):
        cells_dir.mkdir(parents=True, exist_ok=True)
        meta_dir.mkdir(parents=True, exist_ok=True)
        cells.to_csv(cells_dir / "cells.csv.gz", index=False)
        (cells_dir / "cells.json").write_text(json.dumps(meta, indent=2))
        sheet.to_csv(meta_dir / "samples.csv", index=False)
    return {"linked": linked, "flat": root / "flat"}


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("root", type=Path)
    ap.add_argument("--cells", type=int, default=300, help="cells per sample")
    a = ap.parse_args()
    for k, v in make_fixture(a.root, a.cells).items():
        print(f"{k}: {v}")
