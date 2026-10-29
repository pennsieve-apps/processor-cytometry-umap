# processor-cytometry-umap

A two-dimensional **UMAP** map of single cells (McInnes, Healy & Melville 2018; applied to
single-cell and cytometry data by Becht et al. 2018), the visualisation step of the Immune
Profiling & Differential Abundance workflow (`fcs-prep → peacoqc → cleanet → flowsom →
{cytometry-umap, diffcyt} → cytometry-report`).

It takes the per-sample cell subsample that **processor-flowsom** writes, embeds it on the
same markers FlowSOM clustered on, and draws the map coloured by FlowSOM metacluster label,
by condition, by batch, and by the expression of every marker. For an immunologist this is
the familiar "islands" picture: populations as separate islands, activation or
differentiation as gradients within them, and a quick visual check whether batches mix
(they should) and conditions separate (where biology differs). Distances between islands
are not quantitative; use the flowsom and diffcyt tables for numbers.

## Ports

| Port | Kind | Files | Notes |
|---|---|---|---|
| `cells` | input (subfolder) | `cells.csv.gz`, `cells.json` | flowsom's `cells` output: `sample_id, cluster, metacluster, label`, then one column per marker (arcsinh). `cells.json` names the `clustering_markers`. Found in `INPUT_DIR/cells/`, else anywhere under `INPUT_DIR` (symlinks followed). |
| `metadata` | input (subfolder, optional) | `samples.csv` | fcs-prep's sample table; its `condition` and `batch` columns are joined by `sample_id`. Without it both are `NA`. |
| `embedding` | output (the whole `OUTPUT_DIR`) | see below | |

## Parameters

| Parameter | Type | Default | What it does |
|---|---|---|---|
| `N_NEIGHBORS` | integer | `15` | UMAP `n_neighbors` (umap-learn's default): the neighbourhood size. Larger values keep more global structure, smaller ones more local detail. |
| `MIN_DIST` | number | `0.3` | UMAP `min_dist`: how tightly points may pack. umap-learn's default, 0.1, packs dense populations into tight blobs; with 100,000 cells that hides gradients inside a population (naive → memory, resting → activated) under over-plotted points. 0.3 spreads each island enough to see them while keeping populations clearly apart. Lower it for crisper islands, raise it (0.5) for more continuum. |
| `MAX_CELLS` | integer | `100000` | At most this many cells are embedded: a seeded subsample with an equal share per sample (a sample with fewer cells gives its share to the others), so every sample is represented alike. flowsom's default gives 5,000 per sample (120,000 for 24 samples). |
| `MARKERS` | string | empty | Comma list of marker columns to embed on (case and punctuation ignored). Empty: the clustering markers flowsom recorded in `cells.json`; without `cells.json`, every marker column. |
| `SEED` | integer | `42` | Seeds the subsample and UMAP (`random_state`). Same input + same seed = same embedding; it also makes umap-learn run its optimisation single-threaded. |

The marker values are embedded as given (already arcsinh-transformed upstream), with a
Euclidean metric and no rescaling, the same space FlowSOM clustered in.

## Outputs

| File | Content |
|---|---|
| `embedding.csv.gz` | `sample_id, metacluster, label, umap1, umap2, condition, batch`, one row per embedded cell. |
| `umap_by_label.png` | Cells coloured by FlowSOM label; every label is also written at its median position. Eight colours at most: the seven most common labels get one, the rest are grey "other" (still named on the map). |
| `umap_by_condition.png`, `umap_by_batch.png` | The same map by condition and by batch (same colour rule). |
| `umap_markers.png` | One panel per marker, cells coloured light → dark blue by expression, each marker scaled between its 1st and 99th percentile; the markers the embedding used come first, the rest are marked "not used for UMAP". |
| `provenance.json` | App and tool versions, parameters (and the markers used), inputs seen, cells in / embedded, per-sample cell counts, timings. |

## Results on the demo data

flowsom's `cells` output for the full CytoBatchFlagR demo data (24 samples × 5,000 cells),
fcs-prep's `samples.csv`, default parameters: 100,000 cells (4,167 or 4,166 per sample) on
14 markers. The map has four large islands: B cells; CD4 T cells (CD45RA+ and CD45RA−
halves, with the CD25+FOXP3+ cells at one edge); CD8 and DN T cells; and NK cells
(CD56bright, CD56dim CD16+ and CD16int) next to the small monocyte/DC-like groups.
Stimulated and unstimulated cells occupy different parts of the T and NK islands
(activation). Batches mostly overlap, but batch06 and batch01 form sub-regions of their own
in several islands: the batch effects this dataset was published to demonstrate (no
normalisation upstream in this workflow).

## Runtime

`computeTypes: [standard, lambda]`, 2 vCPU / 4 GiB, 15 min timeout.

Measured on the demo data above (100,000 cells × 14 markers) in Docker, **amd64 emulated
on an Apple M-series Mac**: **111–119 s** wall time (UMAP 100–108 s including numba's
just-in-time compilation, figures 10 s), **peak memory 1.33–1.59 GB** (cgroup
`memory.peak`, two runs). That is well inside a Lambda function's limits (10 GB, 15 min),
so the app also declares `lambda`; the image serves both runtimes (`processor/entry.py`
hands over to the Lambda runtime interface client when `AWS_LAMBDA_RUNTIME_API` is set;
`processor/handler.py` maps the payload's `inputDir`, `outputDir` and parameters onto the
environment). numba's cache and matplotlib's config go to `/tmp`, the only writable
directory on Lambda. Raising `MAX_CELLS` well beyond 100,000 should stay on `standard`.

The image is about 600 MB unpacked (195 MB compressed).

## Running it locally

```sh
python3 -m venv .venv && .venv/bin/pip install -r requirements-dev.txt
make test PYTHON=.venv/bin/python   # pytest: both input layouts, symlinks, parameters, Lambda handler (no Docker)
make run PYTHON=.venv/bin/python    # synthetic fixture through the Docker image into data/output/
```

By hand, on flowsom and fcs-prep outputs:

```sh
docker build --platform linux/amd64 -t processor-cytometry-umap .
docker run --rm --platform linux/amd64 -e INPUT_DIR=/data/input -e OUTPUT_DIR=/data/output \
  -v "$PWD/flowsom/cells:/data/input/cells:ro" -v "$PWD/prep/metadata:/data/input/metadata:ro" \
  -v "$PWD/data/output:/data/output" processor-cytometry-umap
```

## Limitations

- UMAP shows neighbourhoods, not distances or densities: island sizes and gaps between
  islands are not measurements, and a different seed gives a differently arranged map.
- The map is of a subsample (flowsom's `CELLS_FOR_PLOT` per sample, then `MAX_CELLS`); rare
  populations may be thin.
- Seeding makes umap-learn single-threaded; that is the price of reproducibility.

## Licences

Apache-2.0 for our code; every bundled component is permissively licensed (`THIRD_PARTY.md`).

## Citations

> Becht E, McInnes L, Healy J, Dutertre CA, Kwok IWH, Ng LG, Ginhoux F, Newell EW.
> **Dimensionality reduction for visualizing single-cell data using UMAP.**
> Nat Biotechnol 37(1):38-44 (2019; online 2018). https://doi.org/10.1038/nbt.4314, PMID 30531897

> McInnes L, Healy J, Melville J. **UMAP: Uniform Manifold Approximation and Projection for
> Dimension Reduction.** arXiv:1802.03426 (2018). https://doi.org/10.48550/arXiv.1802.03426
> (not in PubMed; checked against the DataCite record of the DOI)

FlowSOM (the clusters and labels shown): Van Gassen S et al. Cytometry A 87(7):636-645
(2015), https://doi.org/10.1002/cyto.a.22625, PMID 25573116.

Demo data: Koenig Z, Andorf S, Tilburgs T, Eswar S. **CytoBatchFlagR: Spectral flow
cytometry dataset.** Zenodo (2025). https://doi.org/10.5281/zenodo.15388817 (CC BY 4.0)
