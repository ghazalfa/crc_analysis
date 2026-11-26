#!/usr/bin/env python3
"""PLOT_DE: figures for the differential-expression stage.

  volcano_plots.png        one volcano per contrast, top genes labelled
  top_de_heatmap.png       z-scored expression of the top DE genes, tumours grouped by cluster
  expression_pca.png       PCA of the most variable genes, coloured by mutation cluster
  de_gene_counts.png       number of significant up / down genes per contrast
"""
import argparse

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
from matplotlib.patches import Patch
from scipy.cluster.hierarchy import leaves_list, linkage

from plot_style import (BLUE_RAMP, GRID, INK, INK2, MUTED, NEUTRAL_MID, RED_RAMP, SERIES, SURFACE,
                        cluster_colors, cluster_markers)

UP, DOWN, OTHER = SERIES[1], SERIES[0], "#d3d1ca"
N_LABEL = 8


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--de", nargs="+", required=True, help="de_<A>_vs_<B>.tsv files")
    p.add_argument("--summary", required=True)
    p.add_argument("--top-genes", required=True)
    p.add_argument("--metadata", required=True)
    p.add_argument("--pca", required=True)
    p.add_argument("--pca-variance", required=True)
    p.add_argument("--alpha", type=float, default=0.05)
    p.add_argument("--min-lfc", type=float, default=1.0)
    return p.parse_args()


def contrast_name(path):
    return path.split("/")[-1].removeprefix("de_").removesuffix(".tsv")


def volcano(paths, alpha, min_lfc):
    paths = sorted(p for p in paths if "_vs_" in contrast_name(p))
    fig, axes = plt.subplots(1, len(paths), figsize=(5.2 * len(paths), 4.8), squeeze=False)
    for ax, path in zip(axes[0], paths):
        name = contrast_name(path)
        test, ref = name.split("_vs_")
        r = pd.read_csv(path, sep="\t", index_col=0).dropna(subset=["padj"])
        y = -np.log10(r["padj"].clip(lower=1e-300))
        x = r["log2FoldChange"].clip(-8, 8)
        up = r["significant"] & (r["log2FoldChange"] > 0)
        down = r["significant"] & (r["log2FoldChange"] < 0)
        ax.scatter(x[~(up | down)], y[~(up | down)], s=4, color=OTHER, linewidths=0, rasterized=True)
        ax.scatter(x[up], y[up], s=7, color=UP, linewidths=0, rasterized=True, label=f"Up in {test} ({up.sum()})")
        ax.scatter(x[down], y[down], s=7, color=DOWN, linewidths=0, rasterized=True,
                   label=f"Down in {test} ({down.sum()})")
        for v in (-min_lfc, min_lfc):
            ax.axvline(v, color=MUTED, lw=0.8, ls="--")
        ax.axhline(-np.log10(alpha), color=MUTED, lw=0.8, ls="--")
        # label the strongest genes on each side (score = -log10 padj x |log2FC|), stacked so they don't overlap
        ymax = y.max() * 1.05
        ax.set_ylim(0, ymax)
        xlo, xhi = ax.get_xlim()
        w = xhi - xlo
        xlo, xhi = xlo - 0.15 * w, xhi + 0.15 * w   # room for the stacked labels
        ax.set_xlim(xlo, xhi)
        for mask, side in ((up, 1), (down, -1)):
            top = r[mask].assign(score=y[mask] * r.loc[mask, "log2FoldChange"].abs()) \
                .sort_values("score", ascending=False).head(N_LABEL)
            top = top.assign(px=x[top.index], py=y[top.index]).sort_values("py", ascending=False)
            gap, last = 0.045 * ymax, None
            for gid, row in top.iterrows():
                ty = row["py"] if last is None else min(row["py"], last - gap)
                last = ty
                tx = row["px"] + side * 0.04 * (xhi - xlo)
                ax.annotate(str(row["symbol"]), (row["px"], row["py"]), xytext=(tx, ty),
                            ha="left" if side > 0 else "right", va="center", fontsize=7, color=INK,
                            fontstyle="italic",
                            arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.5, shrinkA=0, shrinkB=2))
        ax.set_xlabel(f"log2 fold change ({test} / {ref})")
        ax.set_ylabel("-log10 adjusted p")
        ax.set_title(f"{test} vs {ref}", loc="left", fontsize=11, pad=22)
        ax.grid(color=GRID, lw=0.5); ax.set_axisbelow(True)
        ax.legend(frameon=False, fontsize=8, loc="lower left", bbox_to_anchor=(0, 1.0), ncol=2,
                  markerscale=2, borderaxespad=0.2)
    fig.text(0.01, 0.005, f"Significant: BH-adjusted p < {alpha} and |log2 fold change| >= {min_lfc}. "
             "\nFold changes beyond +/-8 drawn at the edge.", fontsize=7.5, color=MUTED)
    fig.tight_layout(rect=(0, 0.03, 1, 1)); fig.savefig("volcano_plots.png", dpi=200); plt.close(fig)


def heatmap(top_path, meta):
    top = pd.read_csv(top_path, sep="\t", index_col=0)
    if top.empty:
        return
    symbols = top.pop("symbol").fillna(pd.Series(top.index, index=top.index))
    z = top.sub(top.mean(axis=1), axis=0).div(top.std(axis=1).replace(0, 1), axis=0).clip(-3, 3)

    # tumours grouped by cluster, ordered within cluster by similarity; genes clustered
    order = []
    for c in sorted(meta["cluster"].unique()):
        cols = meta.index[meta["cluster"] == c].intersection(z.columns)
        sub = z[cols]
        order += list(sub.columns[leaves_list(linkage(sub.T.values, "average", metric="correlation"))]) \
            if len(cols) > 2 else list(cols)
    rows = z.index[leaves_list(linkage(z.values, "average", metric="correlation"))]
    z = z.loc[rows, order]

    cmap = LinearSegmentedColormap.from_list("div", BLUE_RAMP[::-1][1:5] + [NEUTRAL_MID] + RED_RAMP)
    colors = cluster_colors(meta["cluster"].unique())
    fig = plt.figure(figsize=(11, 0.13 * len(z) + 2.2))
    gs = fig.add_gridspec(2, 2, height_ratios=[0.25, max(len(z) * 0.13, 1)], width_ratios=[1, 0.02],
                          hspace=0.02, wspace=0.02, left=0.1, right=0.92, top=0.93, bottom=0.04)
    ax = fig.add_subplot(gs[0, 0])
    keys = sorted(colors)
    codes = np.array([keys.index(c) for c in meta.loc[order, "cluster"]])
    ax.imshow(codes[None, :], aspect="auto", cmap=ListedColormap([colors[k] for k in keys]),
              vmin=0, vmax=len(keys) - 1, interpolation="nearest")
    ax.set_xticks([]); ax.set_yticks([0]); ax.set_yticklabels(["Cluster"]); ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_title("Top differentially expressed genes (z-scored log2 normalised counts)", loc="left", fontsize=11)

    ax = fig.add_subplot(gs[1, 0])
    im = ax.imshow(z.values, aspect="auto", cmap=cmap, vmin=-3, vmax=3, interpolation="nearest")
    bounds = np.flatnonzero(codes[1:] != codes[:-1]) + 0.5
    for b in bounds:
        ax.axvline(b, color=SURFACE, lw=2)
    ax.set_xticks([]); ax.set_yticks(range(len(z)))
    ax.set_yticklabels(symbols.loc[rows], fontsize=6, fontstyle="italic"); ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    cax = fig.add_subplot(gs[1, 1])
    cb = fig.colorbar(im, cax=cax); cb.set_label("z-score"); cb.outline.set_visible(False)
    present = [k for k in keys if k in set(meta.loc[order, "cluster"])]
    fig.legend(handles=[Patch(color=colors[k], label=f"{k} (n={(meta['cluster'] == k).sum()})") for k in present],
               loc="upper right", bbox_to_anchor=(0.92, 0.995), ncol=len(present), frameon=False, fontsize=8)
    fig.savefig("top_de_heatmap.png", dpi=200); plt.close(fig)


def pca_plot(pca_path, var_path):
    pcs = pd.read_csv(pca_path, sep="\t", index_col=0)
    var = pd.read_csv(var_path, sep="\t").set_index("pc")["variance_explained"]
    colors, markers = cluster_colors(pcs["cluster"].unique()), cluster_markers(pcs["cluster"].unique())
    fig, ax = plt.subplots(figsize=(6.2, 5))
    for c in sorted(pcs["cluster"].unique()):
        sub = pcs[pcs["cluster"] == c]
        ax.scatter(sub["PC1"], sub["PC2"], s=22, color=colors[c], marker=markers[c], alpha=0.85,
                   edgecolor=SURFACE, linewidths=0.6, label=f"{c} (n={len(sub)})")
    ax.set_xlabel(f"PC1 ({100 * var['PC1']:.1f}% of variance)")
    ax.set_ylabel(f"PC2 ({100 * var['PC2']:.1f}% of variance)")
    ax.grid(color=GRID, lw=0.5); ax.set_axisbelow(True)
    ax.legend(frameon=False, fontsize=8, title="Mutation cluster", title_fontsize=8.5)
    ax.set_title("Expression PCA (most variable genes)", loc="left", fontsize=11)
    fig.tight_layout(); fig.savefig("expression_pca.png", dpi=200); plt.close(fig)


def counts_plot(summary_path):
    s = pd.read_csv(summary_path, sep="\t")
    y = np.arange(len(s))
    fig, ax = plt.subplots(figsize=(6.5, 0.6 * len(s) + 1.4))
    ax.barh(y, s["up"], color=UP, height=0.55, label="Up")
    ax.barh(y, -s["down"], color=DOWN, height=0.55, label="Down")
    for yi, (u, d) in enumerate(zip(s["up"], s["down"])):
        ax.text(u, yi, f" {u}", va="center", fontsize=8, color=INK2)
        ax.text(-d, yi, f"{d} ", va="center", ha="right", fontsize=8, color=INK2)
    ax.axvline(0, color=INK2, lw=0.8)
    ax.set_yticks(y); ax.set_yticklabels([c.replace("_vs_", " vs ") for c in s["contrast"]])
    lim = max(s["up"].max(), s["down"].max()) * 1.25 + 1
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{abs(v):.0f}"))
    ax.set_xlim(-lim, lim); ax.set_xlabel("Significant genes (down | up in the first cluster)")
    ax.grid(axis="x", color=GRID, lw=0.5); ax.set_axisbelow(True)
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    ax.set_title("Differentially expressed genes per contrast", loc="left", fontsize=11)
    fig.tight_layout(); fig.savefig("de_gene_counts.png", dpi=200); plt.close(fig)


def main():
    a = parse_args()
    meta = pd.read_csv(a.metadata, sep="\t", index_col=0)
    volcano(a.de, a.alpha, a.min_lfc)
    heatmap(a.top_genes, meta)
    pca_plot(a.pca, a.pca_variance)
    counts_plot(a.summary)


if __name__ == "__main__":
    main()
