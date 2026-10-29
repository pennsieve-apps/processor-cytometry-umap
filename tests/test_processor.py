"""End-to-end on the synthetic fixture, both input layouts, plus the Lambda handler and helpers."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from processor import handler
from processor import main as app
from tests.make_fixture import CLUSTERING, make_fixture

FIGURES = ("umap_by_label.png", "umap_by_condition.png", "umap_by_batch.png", "umap_markers.png")


@pytest.fixture(scope="session")
def fixture_root(tmp_path_factory) -> dict[str, Path]:
    return make_fixture(tmp_path_factory.mktemp("fixture"))


def run(input_dir: Path, output_dir: Path, monkeypatch, **params) -> int:
    monkeypatch.setenv("INPUT_DIR", str(input_dir))
    monkeypatch.setenv("OUTPUT_DIR", str(output_dir))
    for name in app.PARAMETERS:
        monkeypatch.delenv(name, raising=False)
    for name, value in params.items():
        monkeypatch.setenv(name, str(value))
    return app.main()


@pytest.mark.parametrize("layout", ["linked", "flat"])
def test_layouts(fixture_root, tmp_path, monkeypatch, layout):
    out = tmp_path / "out"
    assert run(fixture_root[layout], out, monkeypatch) == 0
    emb = pd.read_csv(out / "embedding.csv.gz")
    assert list(emb.columns) == ["sample_id", "metacluster", "label", "umap1", "umap2", "condition", "batch"]
    assert len(emb) == 1200
    assert set(emb["condition"]) == {"stimulated", "unstimulated"}
    assert set(emb["batch"]) == {"batch01", "batch02"}
    for f in FIGURES:
        assert (out / f).stat().st_size > 10_000, f
    prov = json.loads((out / "provenance.json").read_text())
    assert prov["app"] == "processor-cytometry-umap" and prov["version"] == "0.1.0"
    assert prov["parameters"]["MARKERS_used"] == CLUSTERING  # from cells.json
    assert prov["parameters"]["MIN_DIST"] == 0.3 and prov["parameters"]["N_NEIGHBORS"] == 15
    assert all(not p.startswith("/") for p in prov["inputs"])
    assert len(prov["samples"]) == 4


def test_populations_separate(fixture_root, tmp_path, monkeypatch):
    out = tmp_path / "out"
    assert run(fixture_root["linked"], out, monkeypatch) == 0
    emb = pd.read_csv(out / "embedding.csv.gz")
    centres = emb.groupby("label")[["umap1", "umap2"]].median()
    spread = emb.groupby("label")[["umap1", "umap2"]].std().max().max()
    d = np.linalg.norm(centres.to_numpy()[:, None] - centres.to_numpy()[None], axis=-1)
    assert d[np.triu_indices(len(centres), 1)].min() > 2 * spread


def test_reproducible_and_max_cells(fixture_root, tmp_path, monkeypatch):
    a, b = tmp_path / "a", tmp_path / "b"
    assert run(fixture_root["linked"], a, monkeypatch, MAX_CELLS=400, SEED=7) == 0
    assert run(fixture_root["linked"], b, monkeypatch, MAX_CELLS=400, SEED=7) == 0
    ea, eb = pd.read_csv(a / "embedding.csv.gz"), pd.read_csv(b / "embedding.csv.gz")
    assert len(ea) == 400
    assert ea["sample_id"].value_counts().tolist() == [100] * 4  # equal share per sample
    pd.testing.assert_frame_equal(ea, eb)


def test_markers_override_and_errors(fixture_root, tmp_path, monkeypatch):
    out = tmp_path / "out"
    assert run(fixture_root["linked"], out, monkeypatch, MARKERS="cd3, CD19,CD56") == 0
    prov = json.loads((out / "provenance.json").read_text())
    assert prov["parameters"]["MARKERS_used"] == ["CD3", "CD19", "CD56"]
    assert run(fixture_root["linked"], tmp_path / "bad", monkeypatch, MARKERS="CD3,NOPE") == 1
    assert run(fixture_root["linked"], tmp_path / "bad2", monkeypatch, MIN_DIST="2") == 1
    assert run(tmp_path / "nothing-here", tmp_path / "bad3", monkeypatch) == 1


def test_without_metadata_or_cells_json(tmp_path, monkeypatch):
    root = make_fixture(tmp_path / "fx")["linked"]
    (root / "metadata" / "samples.csv").unlink()
    (root / "cells" / "cells.json").unlink()
    out = tmp_path / "out"
    assert run(root, out, monkeypatch) == 0
    emb = pd.read_csv(out / "embedding.csv.gz", keep_default_na=False)
    assert set(emb["condition"]) == {"NA"} and set(emb["batch"]) == {"NA"}
    prov = json.loads((out / "provenance.json").read_text())
    assert prov["parameters"]["MARKERS_used"] == CLUSTERING + ["CD69"]  # every marker column


def test_symlinked_inputs(fixture_root, tmp_path, monkeypatch):
    linked = tmp_path / "in"
    linked.mkdir()
    for port in ("cells", "metadata"):
        (linked / port).symlink_to(fixture_root["linked"] / port, target_is_directory=True)
    assert run(linked, tmp_path / "out", monkeypatch) == 0


def test_subsample_water_filling():
    df = pd.DataFrame({"sample_id": ["a"] * 10 + ["b"] * 100 + ["c"] * 100, "x": range(210)})
    sub = app.subsample(df, 110, 1)
    assert sub["sample_id"].value_counts().to_dict() == {"a": 10, "b": 50, "c": 50}
    assert len(app.subsample(df, 103, 1)) == 103  # remainder handed out too


def test_lambda_handler(fixture_root, tmp_path, monkeypatch):
    for name in app.PARAMETERS:
        monkeypatch.delenv(name, raising=False)
    out = tmp_path / "out"
    res = handler.handler({"inputDir": str(fixture_root["linked"]), "outputDir": str(out), "N_NEIGHBORS": 10}, None)
    assert res["status"] == "ok"
    assert json.loads((out / "provenance.json").read_text())["parameters"]["N_NEIGHBORS"] == 10
    with pytest.raises(RuntimeError):
        handler.handler({"inputDir": str(tmp_path / "none"), "outputDir": str(tmp_path / "o2")}, None)
