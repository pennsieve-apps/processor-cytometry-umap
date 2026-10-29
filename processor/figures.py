"""UMAP scatter figures (matplotlib, Agg, PNG)."""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.patheffects as pe  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

# Eight categorical slots in fixed order (colour-vision-deficiency checked set);
# a ninth category is never a new hue: the smallest ones fold into grey "other".
CATEGORICAL = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
OTHER = "#b4b2ab"
INK = "#0b0b0b"
INK_2 = "#52514e"
SURFACE = "#fcfcfb"
# One hue, light -> dark, for expression.
SEQUENTIAL = LinearSegmentedColormap.from_list(
    "blue", ["#e4eefa", "#9ec5f4", "#5598e7", "#2a78d6", "#1c5cab", "#104281", "#0d366b"])
DPI = 110


def _point_size(n: int) -> float:
    return float(np.clip(4000.0 / max(n, 1), 0.3, 6.0))


def _style(ax) -> None:
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_facecolor(SURFACE)


def categorical(df: pd.DataFrame, column: str, title: str, path: Path, direct_labels: bool = False) -> dict:
    values = df[column].astype(str)
    counts = values.value_counts()
    if len(counts) <= len(CATEGORICAL):
        shown = list(counts.index)
    else:
        shown = list(counts.index[: len(CATEGORICAL) - 1])
    # Stable hue per category: alphabetical among the shown ones, so the same
    # label keeps its colour between runs with similar abundance.
    shown_sorted = sorted(shown)
    colour = {c: CATEGORICAL[i] for i, c in enumerate(shown_sorted)}
    folded = [c for c in counts.index if c not in colour]
    cols = np.array([colour.get(v, OTHER) for v in values])
    # Draw in a seeded random order so no category hides under another.
    order = np.random.default_rng(0).permutation(len(df))
    fig, ax = plt.subplots(figsize=(9, 7.5), dpi=DPI, facecolor=SURFACE)
    ax.scatter(df["umap1"].to_numpy()[order], df["umap2"].to_numpy()[order], c=cols[order],
               s=_point_size(len(df)), linewidths=0, rasterized=True)
    _style(ax)
    handles = [plt.Line2D([], [], marker="o", ls="", color=colour[c], markersize=8,
                          label=f"{c} ({counts[c]:,})") for c in shown_sorted]
    if folded:
        handles.append(plt.Line2D([], [], marker="o", ls="", color=OTHER, markersize=8,
                                  label=f"{len(folded)} other ({int(counts[folded].sum()):,})"))
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.0, 1.0), frameon=False,
              fontsize=8, labelcolor=INK_2, title=title, title_fontsize=9)
    if direct_labels:
        # Every category is named at its median position, folded ones included,
        # so identity never rests on colour alone.
        for c, g in df.groupby(values):
            if len(g) < max(20, 0.002 * len(df)):
                continue
            ax.text(g["umap1"].median(), g["umap2"].median(), c, fontsize=7, color=INK, ha="center", va="center",
                    path_effects=[pe.withStroke(linewidth=2.5, foreground=SURFACE)])
    ax.set_title(f"UMAP of {len(df):,} cells by {title.lower()}", loc="left", fontsize=11, color=INK)
    fig.tight_layout()
    fig.savefig(path, facecolor=SURFACE)
    plt.close(fig)
    return {"column": column, "categories": int(len(counts)), "coloured": shown_sorted, "folded_into_other": folded}


def marker_grid(df: pd.DataFrame, cells: pd.DataFrame, used: list[str], others: list[str], path: Path) -> dict:
    markers = list(used) + list(others)
    ncol = 6 if len(markers) > 12 else 4
    nrow = int(np.ceil(len(markers) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(2.6 * ncol, 2.5 * nrow + 0.6), dpi=DPI, facecolor=SURFACE, squeeze=False)
    x = df["umap1"].to_numpy()
    y = df["umap2"].to_numpy()
    size = _point_size(len(df)) * 0.5
    for ax, m in zip(axes.ravel(), markers):
        v = cells[m].to_numpy(dtype=float)
        lo, hi = np.nanpercentile(v, [1, 99])
        if hi <= lo:
            hi = lo + 1e-6
        order = np.argsort(v)  # high values on top
        ax.scatter(x[order], y[order], c=np.clip(v[order], lo, hi), cmap=SEQUENTIAL, vmin=lo, vmax=hi,
                   s=size, linewidths=0, rasterized=True)
        _style(ax)
        ax.set_title(f"{m}{'' if m in used else ' (not used for UMAP)'}", fontsize=8,
                     color=INK if m in used else INK_2, loc="left")
        ax.text(0.99, 0.01, f"{lo:.1f} to {hi:.1f}", transform=ax.transAxes, fontsize=6, color=INK_2, ha="right", va="bottom")
    for ax in axes.ravel()[len(markers):]:
        ax.axis("off")
    fig.suptitle("Marker expression (arcsinh), colour scale per marker from its 1st to 99th percentile, light = low",
                 x=0.01, ha="left", fontsize=10, color=INK)
    fig.tight_layout()
    fig.savefig(path, facecolor=SURFACE)
    plt.close(fig)
    return {"markers": markers}
