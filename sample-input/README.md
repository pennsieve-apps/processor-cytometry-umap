# sample-input

processor-cytometry-umap embeds the cell subsample that processor-flowsom writes to its
`cells` output (`cells.csv.gz` + `cells.json`), optionally with fcs-prep's `metadata`
(`samples.csv`) for condition and batch. There is no public dataset to download for it
directly: run flowsom on the demo data (see processor-flowsom's `sample-input/README.md`)
and point `INPUT_DIR` at a folder holding `cells/` (flowsom's output) and `metadata/`
(fcs-prep's), or at a folder where they sit side by side.

For a local run without upstream outputs, generate the synthetic fixture:

```sh
make fixture PYTHON=.venv/bin/python   # data/fixture/linked/ and data/fixture/flat/
make run PYTHON=.venv/bin/python       # runs the Docker image over data/fixture/linked/
```

`tests/make_fixture.py` writes four populations (CD4 T, CD8 T, B, NK) in four samples, two
batches, two conditions, with higher CD69 on stimulated T cells.
