"""
processor-cytometry-umap: a 2-D UMAP embedding (McInnes et al.; Becht et al. 2018
for single-cell data) of the cell subsample that processor-flowsom writes, with
figures coloured by FlowSOM label, condition, batch and marker expression.

  INPUT_DIR/cells/cells.csv.gz      sample_id, cluster, metacluster, label, <marker>...
  INPUT_DIR/cells/cells.json        {"clustering_markers": [...]} (optional)
  INPUT_DIR/metadata/samples.csv    sample_id, condition, batch (optional)

  OUTPUT_DIR/embedding.csv.gz       sample_id, metacluster, label, umap1, umap2, condition, batch
  OUTPUT_DIR/umap_by_label.png, umap_by_condition.png, umap_by_batch.png, umap_markers.png
  OUTPUT_DIR/provenance.json

Files are found in their port folder when linked, else anywhere under INPUT_DIR
by name (symlinks followed).
"""
from __future__ import annotations

import json
import logging
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
log = logging.getLogger("processor-cytometry-umap")
# The Lambda runtime configures the root logger at WARNING before basicConfig can.
log.setLevel(logging.INFO)
for noisy in ("matplotlib", "PIL", "numba", "fontTools"):
    logging.getLogger(noisy).setLevel(logging.WARNING)

APP_NAME = "processor-cytometry-umap"
PARAMETERS = ("N_NEIGHBORS", "MIN_DIST", "MAX_CELLS", "MARKERS", "SEED")
DEFAULTS = {"N_NEIGHBORS": "15", "MIN_DIST": "0.3", "MAX_CELLS": "100000", "MARKERS": "", "SEED": "42"}
ID_COLUMNS = ("sample_id", "cluster", "metacluster", "label")
REPO_ROOT = Path(__file__).resolve().parent.parent


class InputError(Exception):
    pass


def params() -> dict:
    raw = {k: (os.environ.get(k) or DEFAULTS[k]).strip() for k in PARAMETERS}
    try:
        p = {
            "N_NEIGHBORS": int(raw["N_NEIGHBORS"]),
            "MIN_DIST": float(raw["MIN_DIST"]),
            "MAX_CELLS": int(raw["MAX_CELLS"]),
            "MARKERS": raw["MARKERS"],
            "SEED": int(raw["SEED"]),
        }
    except ValueError as e:
        raise InputError(f"bad parameter: {e}") from e
    if p["N_NEIGHBORS"] < 2:
        raise InputError("N_NEIGHBORS must be at least 2")
    if not 0.0 <= p["MIN_DIST"] <= 1.0:
        raise InputError("MIN_DIST must be between 0 and 1")
    if p["MAX_CELLS"] < 10:
        raise InputError("MAX_CELLS must be at least 10")
    return p


def find(input_dir: Path, port: str, name: str) -> Path | None:
    """The file in INPUT_DIR/<port>/ when linked, else anywhere under INPUT_DIR."""
    for root in (input_dir / port, input_dir):
        if not root.is_dir():
            continue
        hits = sorted(Path(d) / name for d, _, files in os.walk(root, followlinks=True) if name in files)
        if hits:
            return hits[0]
    return None


def rel(path: Path | None, root: Path) -> str | None:
    if path is None:
        return None
    try:
        return str(Path(path).relative_to(root))
    except ValueError:
        return str(path)


def norm(name: str) -> str:
    return "".join(c for c in name.upper() if c.isalnum())


def choose_markers(columns: list[str], spec: str, recorded: list[str] | None) -> list[str]:
    available = [c for c in columns if c not in ID_COLUMNS]
    by_norm = {norm(c): c for c in available}
    if spec.strip():
        wanted = [w.strip() for w in spec.split(",") if w.strip()]
        source = "MARKERS"
    elif recorded:
        wanted, source = recorded, "cells.json clustering_markers"
    else:
        return available
    missing = [w for w in wanted if norm(w) not in by_norm]
    if missing:
        raise InputError(f"{source}: not in cells.csv.gz: {', '.join(missing)} (have {', '.join(available)})")
    return list(dict.fromkeys(by_norm[norm(w)] for w in wanted))


def subsample(cells: pd.DataFrame, max_cells: int, seed: int) -> pd.DataFrame:
    """At most max_cells, the same number from every sample where possible (seeded)."""
    if len(cells) <= max_cells:
        return cells.reset_index(drop=True)
    rng = np.random.default_rng(seed)
    groups = {s: idx for s, idx in cells.groupby("sample_id", sort=True).indices.items()}
    # Water-filling: equal shares, small samples give their leftover to the others.
    quota = {s: 0 for s in groups}
    left = max_cells
    open_ = sorted(groups, key=lambda s: len(groups[s]))
    while left > 0 and open_:
        share = left // len(open_)
        if share == 0:
            break
        nxt = []
        for s in open_:
            add = min(share, len(groups[s]) - quota[s])
            quota[s] += add
            left -= add
            if quota[s] < len(groups[s]):
                nxt.append(s)
        open_ = nxt
    # The remainder (fewer cells than samples left), one each, in sample order.
    for s in sorted(open_):
        if left == 0:
            break
        quota[s] += 1
        left -= 1
    take = np.concatenate([np.sort(rng.choice(groups[s], quota[s], replace=False)) for s in sorted(groups) if quota[s]])
    return cells.iloc[np.sort(take)].reset_index(drop=True)


def main() -> int:
    t0 = time.time()
    started = datetime.now(timezone.utc).isoformat(timespec="seconds")
    input_dir = Path(os.environ.get("INPUT_DIR") or "")
    output_dir = Path(os.environ.get("OUTPUT_DIR") or "")
    version = (REPO_ROOT / "VERSION").read_text().strip() if (REPO_ROOT / "VERSION").exists() else "unknown"
    try:
        if not str(input_dir) or not str(output_dir) or str(input_dir) == "." or str(output_dir) == ".":
            raise InputError("INPUT_DIR and OUTPUT_DIR must be set")
        p = params()
        output_dir.mkdir(parents=True, exist_ok=True)
        cells_path = find(input_dir, "cells", "cells.csv.gz")
        if cells_path is None:
            raise InputError(f"no cells.csv.gz under {input_dir}")
        meta_path = find(input_dir, "cells", "cells.json")
        samples_path = find(input_dir, "metadata", "samples.csv")
        recorded = None
        if meta_path is not None:
            recorded = json.loads(meta_path.read_text()).get("clustering_markers")
        log.info("%s %s: cells %s, cells.json %s, samples.csv %s", APP_NAME, version,
                 rel(cells_path, input_dir), rel(meta_path, input_dir), rel(samples_path, input_dir))

        cells = pd.read_csv(cells_path, dtype={"sample_id": str, "label": str})
        if "sample_id" not in cells.columns:
            raise InputError("cells.csv.gz has no sample_id column")
        for col in ("metacluster", "label"):
            if col not in cells.columns:
                cells[col] = "NA" if col == "label" else -1
        markers = choose_markers(list(cells.columns), p["MARKERS"], recorded)
        if len(markers) < 2:
            raise InputError("need at least two markers for UMAP")
        n_in = len(cells)
        cells = cells.dropna(subset=markers)
        cells = subsample(cells, p["MAX_CELLS"], p["SEED"])
        if len(cells) <= p["N_NEIGHBORS"]:
            raise InputError(f"only {len(cells)} cells, need more than N_NEIGHBORS={p['N_NEIGHBORS']}")

        cond = pd.Series("NA", index=cells.index)
        batch = pd.Series("NA", index=cells.index)
        if samples_path is not None:
            sheet = pd.read_csv(samples_path, dtype=str).drop_duplicates("sample_id").set_index("sample_id")
            for col, target in (("condition", "cond"), ("batch", "batch")):
                if col in sheet.columns:
                    vals = cells["sample_id"].map(sheet[col]).fillna("NA")
                    if target == "cond":
                        cond = vals
                    else:
                        batch = vals
            unknown = sorted(set(cells["sample_id"]) - set(sheet.index))
            if unknown:
                log.warning("samples not in samples.csv: %s", ", ".join(unknown))
        per_sample = cells["sample_id"].value_counts().sort_index()
        for s, n in per_sample.items():
            log.info("%s: %d cells", s, n)
        log.info("UMAP on %d of %d cells from %d samples, %d markers (%s), n_neighbors=%d min_dist=%g seed=%d",
                 len(cells), n_in, len(per_sample), len(markers), ", ".join(markers),
                 p["N_NEIGHBORS"], p["MIN_DIST"], p["SEED"])

        t1 = time.time()
        import umap  # slow import (numba); after the cheap checks

        model = umap.UMAP(n_neighbors=p["N_NEIGHBORS"], min_dist=p["MIN_DIST"], n_components=2,
                          metric="euclidean", random_state=p["SEED"], n_jobs=1)
        emb = model.fit_transform(cells[markers].to_numpy(dtype=np.float32))
        t2 = time.time()
        log.info("UMAP done in %.1fs", t2 - t1)

        out = pd.DataFrame({
            "sample_id": cells["sample_id"].to_numpy(),
            "metacluster": cells["metacluster"].to_numpy(),
            "label": cells["label"].fillna("NA").to_numpy(),
            "umap1": np.round(emb[:, 0], 4),
            "umap2": np.round(emb[:, 1], 4),
            "condition": cond.to_numpy(),
            "batch": batch.to_numpy(),
        })
        out.to_csv(output_dir / "embedding.csv.gz", index=False, compression={"method": "gzip", "mtime": 0})

        from processor import figures

        figs = {}
        figs["umap_by_label.png"] = figures.categorical(out, "label", "FlowSOM metacluster label", output_dir / "umap_by_label.png", direct_labels=True)
        figs["umap_by_condition.png"] = figures.categorical(out, "condition", "Condition", output_dir / "umap_by_condition.png")
        figs["umap_by_batch.png"] = figures.categorical(out, "batch", "Batch", output_dir / "umap_by_batch.png")
        others = [c for c in cells.columns if c not in ID_COLUMNS and c not in markers]
        figs["umap_markers.png"] = figures.marker_grid(out, cells, markers, others, output_dir / "umap_markers.png")
        t3 = time.time()

        import matplotlib
        import numba
        import pynndescent
        import sklearn

        provenance = {
            "app": APP_NAME,
            "version": version,
            "tools": {"python": sys.version.split()[0], "umap-learn": _version("umap-learn"),
                      "pynndescent": _version("pynndescent"), "numba": numba.__version__, "numpy": np.__version__,
                      "pandas": pd.__version__, "scikit-learn": sklearn.__version__, "matplotlib": matplotlib.__version__},
            "parameters": {**p, "MARKERS_used": markers},
            "method": ("umap.UMAP(n_components=2, metric='euclidean') on the marker values as given "
                       "(arcsinh-transformed upstream, not rescaled); random_state=SEED, which makes the "
                       "embedding reproducible and runs the optimisation single-threaded"),
            "inputs": [r for r in (rel(cells_path, input_dir), rel(meta_path, input_dir), rel(samples_path, input_dir)) if r],
            "cells_in": n_in,
            "cells_embedded": len(out),
            "samples": [{"sample_id": s, "cells": int(n), "status": "ok"} for s, n in per_sample.items()],
            "figures": figs,
            "timings": {"read_seconds": round(t1 - t0, 2), "umap_seconds": round(t2 - t1, 2),
                        "figures_seconds": round(t3 - t2, 2), "total_seconds": round(time.time() - t0, 2)},
            "startedAt": started,
            "finishedAt": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        }
        with open(output_dir / "provenance.json", "w") as fh:
            json.dump(provenance, fh, indent=2)
        log.info("done: %d cells embedded, %.1fs", len(out), time.time() - t0)
        return 0
    except InputError as e:
        log.error("%s", e)
        return 1
    except Exception:  # noqa: BLE001
        log.exception("failed")
        return 1


def _version(dist: str) -> str:
    from importlib.metadata import PackageNotFoundError, version

    try:
        return version(dist)
    except PackageNotFoundError:
        return "unknown"


if __name__ == "__main__":
    sys.exit(main())
