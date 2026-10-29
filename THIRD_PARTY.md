# Third-party software

This app is Apache-2.0 (see `LICENSE`). The image it builds contains the following
third-party components, all under permissive licences; their licences apply to the image.

| Component | Version | What | Licence | Source |
|---|---|---|---|---|
| umap-learn | 0.5.9.post2 | UMAP embedding | BSD-3-Clause | https://github.com/lmcinnes/umap |
| pynndescent | 0.5.13 | Approximate nearest neighbours for UMAP | BSD-2-Clause | https://github.com/lmcinnes/pynndescent |
| numba / llvmlite | 0.62.1 / 0.45.1 | JIT compiler used by umap-learn (bundles LLVM, Apache-2.0 with LLVM exceptions) | BSD-2-Clause | https://numba.pydata.org |
| scikit-learn | 1.7.2 | Used by umap-learn | BSD-3-Clause | https://scikit-learn.org |
| NumPy, SciPy, pandas | 2.3.4, 1.16.3, 2.3.3 | Arrays, CSV | BSD-3-Clause (NumPy/SciPy bundle OpenBLAS, BSD-3, and gfortran runtime, GPL-3 with the GCC Runtime Library Exception) | PyPI |
| matplotlib | 3.10.7 | Figures | Matplotlib licence (PSF-based, BSD-compatible) | https://matplotlib.org |
| joblib, threadpoolctl, contourpy, cycler, kiwisolver, pyparsing, python-dateutil, pytz, tzdata, six, packaging, tqdm, cloudpickle, fonttools, pillow | see `requirements.txt` | Transitive dependencies | BSD-3 / BSD-2 / MIT / Apache-2.0 / MPL-2.0 (tqdm) / MIT-CMU (pillow) | PyPI |
| awslambdaric | 3.1.1 (+ simplejson 4.2.0, snapshot-restore-py 1.0.0) | Lambda runtime interface client | Apache-2.0 (simplejson MIT) | https://github.com/aws/aws-lambda-python-runtime-interface-client |
| python:3.12-slim-bookworm | 3.12.16 | Base image (Debian 12 + CPython) | PSF; Debian packages under their own free licences | https://hub.docker.com/_/python |

The colour palette in the figures is a colour-vision-deficiency-checked eight-colour set
chosen for this app; no third-party assets are bundled.

## Citations

See README.md, Citations.
