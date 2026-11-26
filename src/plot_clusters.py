#!/usr/bin/env python3
"""PLOT_CLUSTERS: figures for visual inspection of the mutation clustering.

  oncoprint.png / .pdf          driver-gene oncoprint ordered by the dendrogram, with TMB and
                                cluster / hypermutation / side / stage tracks
  k_selection.png               mean silhouette for each k (chosen k highlighted)
  cluster_gene_frequency.png    % of tumours mutated per gene in each cluster (* = enriched, q < 0.05)
  cluster_characteristics.png   per-cluster hypermutation, READ share and bootstrap stability
"""
import argparse
import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from matplotlib.patches import Patch
from scipy.cluster.hierarchy import dendrogram

from plot_style import (AXIS, BLUE_RAMP, EMPTY, GRID, INK, INK2, MUTED, NO_DRIVER, SERIES,  # noqa: F401
                        SURFACE, cluster_colors)

MUT_COLORS = {"Truncating": SERIES[0], "Missense": SERIES[1], "In-frame indel": SERIES[2], "Multi-hit": INK}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for name in ["types", "annotation", "assignments", "linkage", "linkage-samples", "k-selection",
                 "enrichment", "summary", "cluster-info", "genes"]:
        p.add_argument(f"--{name}", required=True)
    p.add_argument("--hypermut-tmb", type=float, default=10.0)
    return p.parse_args()


def track(ax, values, palette, label):
    keys = list(palette)
    codes = np.array([keys.index(v) if v in palette else -1 for v in values])
    cmap = ListedColormap([palette[k] for k in keys] + [EMPTY])
    codes = np.where(codes < 0, len(keys), codes)
    ax.imshow(codes[None, :], aspect="auto", cmap=cmap, vmin=0, vmax=len(keys), interpolation="nearest")
    ax.set_xticks([]); ax.set_yticks([0]); ax.set_yticklabels([label])
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)


def oncoprint(types, ann, assign, Z, linked, genes_meta, info, tmb_cutoff):
    # column order: dendrogram leaves, then tumours without a driver mutation
    leaves = dendrogram(Z, no_plot=True)["leaves"]
    order = [linked[i] for i in leaves] + [s for s in types.index if assign[s] == NO_DRIVER]
    genes = (types != "").mean().sort_values(ascending=False).index
    T = types.loc[order, genes].fillna("")
    a = ann.loc[order]
    n_genes = len(genes)

    tracks = [
        ("Cluster", assign[order].values, cluster_colors(assign.unique())),
        ("Hypermutated", a["hypermutated"].values, {"Hypermutated": INK, "Non-hypermutated": "#d8d6cf"}),
        # anatomical order proximal -> distal, so an ordered grey ramp
        ("Side", a["side"].values, {"Right colon": "#3a3936", "Left colon": "#898781", "Rectum": "#c9c7c0"}),
        ("Stage", a["stage"].values, {"I": BLUE_RAMP[1], "II": BLUE_RAMP[3], "III": BLUE_RAMP[4], "IV": BLUE_RAMP[6]}),
    ]
    heights = [1.6, 1.0] + [0.28] * len(tracks) + [0.15] + [0.3] * n_genes
    fig = plt.figure(figsize=(15, 2.6 + 0.3 * (n_genes + len(tracks)) + 1.6))
    gs = fig.add_gridspec(len(heights), 2, height_ratios=heights, width_ratios=[1, 0.07],
                          hspace=0.08, wspace=0.01, left=0.08, right=0.86, top=0.95, bottom=0.03)

    # dendrogram (only over clustered tumours; the no-driver block sits to its right)
    n_link = len(linked)
    ax = fig.add_subplot(gs[0, 0])
    dendrogram(Z, ax=ax, no_labels=True, color_threshold=0, above_threshold_color=INK2)
    ax.set_xlim(0, len(order) * 10)
    ax.axis("off")
    ax.set_title(f"Driver-mutation clustering of TCGA COAD + READ primary tumours "
                 f"({info['distance'].capitalize()} distance, {info['linkage']} linkage, k = {info['k']})",
                 loc="left", fontsize=11, color=INK)

    # TMB bar (log scale)
    ax = fig.add_subplot(gs[1, 0])
    ax.bar(np.arange(len(order)), a["tmb"].clip(lower=0.05), width=1.0, color=INK2, linewidth=0)
    ax.set_yscale("log"); ax.set_xlim(-0.5, len(order) - 0.5); ax.set_xticks([])
    ax.axhline(tmb_cutoff, color=INK, lw=1, ls="--")
    ax.set_ylabel("TMB\n(mut/Mb)", rotation=0, ha="right", va="center")
    ax.grid(axis="y", color=GRID, lw=0.5); ax.set_axisbelow(True)

    for i, (label, vals, pal) in enumerate(tracks):
        track(fig.add_subplot(gs[2 + i, 0]), vals, pal, label)

    # cluster boundaries
    bounds = np.flatnonzero(assign[order].values[1:] != assign[order].values[:-1]) + 0.5

    # gene rows
    first = 3 + len(tracks)
    freq = (T != "").mean()
    for j, g in enumerate(genes):
        ax = fig.add_subplot(gs[first + j, 0])
        row = T[g].values
        keys = list(MUT_COLORS)
        codes = np.array([keys.index(v) if v else len(keys) for v in row])
        ax.imshow(codes[None, :], aspect="auto", interpolation="nearest",
                  cmap=ListedColormap([MUT_COLORS[k] for k in keys] + [EMPTY]), vmin=0, vmax=len(keys))
        for b in bounds:
            ax.axvline(b, color=SURFACE, lw=1.5)
        ax.set_xticks([]); ax.set_yticks([0]); ax.set_yticklabels([g], fontstyle="italic")
        ax.tick_params(length=0)
        for s in ax.spines.values():
            s.set_visible(False)
        fax = fig.add_subplot(gs[first + j, 1])
        fax.barh([0], [freq[g]], color=INK2, height=0.6)
        fax.set_xlim(0, 1); fax.set_ylim(-0.5, 0.5); fax.axis("off")
        fax.text(freq[g] + 0.03, 0, f"{100 * freq[g]:.0f}%", va="center", fontsize=7.5, color=INK2)

    # legends
    handles = [Patch(color=c, label=k) for k, c in MUT_COLORS.items()] + [Patch(color=EMPTY, label="Not mutated")]
    fig.legend(handles=handles, title="Mutation", loc="upper left", bbox_to_anchor=(0.875, 0.55),
               frameon=False, fontsize=8, title_fontsize=8.5)
    y = 0.95
    for label, _, pal in tracks:
        h = [Patch(color=c, label=k) for k, c in pal.items()]
        leg = fig.legend(handles=h, title=label, loc="upper left", bbox_to_anchor=(0.875, y),
                         frameon=False, fontsize=7.5, title_fontsize=8.5, ncol=2 if len(h) > 3 else 1)
        y -= 0.025 + 0.018 * ((len(h) + 1) // 2 if len(h) > 3 else len(h)) + 0.012
        fig.add_artist(leg)
    fig.text(0.08, 0.005, f"{n_link} clustered tumours; {len(order) - n_link} tumours with no mutation in "
             f"these genes shown at right (grey cluster). Dashed line: hypermutation threshold, {tmb_cutoff:g} mut/Mb.",
             fontsize=7.5, color=MUTED)
    fig.savefig("oncoprint.png", dpi=200)
    fig.savefig("oncoprint.pdf")
    plt.close(fig)


def k_selection_plot(ksel):
    fig, ax = plt.subplots(figsize=(5.5, 3.6))
    ax.plot(ksel["k"], ksel["mean_silhouette"], color=SERIES[0], lw=2, marker="o", ms=7, zorder=3)
    ch = ksel[ksel["chosen"]]
    ax.scatter(ch["k"], ch["mean_silhouette"], s=140, facecolor="none", edgecolor=INK, lw=1.5, zorder=4)
    for _, r in ch.iterrows():
        ax.annotate(f"chosen k = {r['k']}", (r["k"], r["mean_silhouette"]), xytext=(8, 6),
                    textcoords="offset points", color=INK, fontsize=8.5)
    inel = ksel[~ksel["eligible"]]
    if len(inel):
        ax.scatter(inel["k"], inel["mean_silhouette"], marker="x", color=MUTED, s=50, zorder=5,
                   label="smallest cluster below minimum size")
        ax.legend(frameon=False, fontsize=8)
    ax.set_xlabel("Number of clusters (k)"); ax.set_ylabel("Mean silhouette width")
    ax.set_xticks(ksel["k"]); ax.grid(axis="y", color=GRID, lw=0.5); ax.set_axisbelow(True)
    ax.set_title("Choosing k", loc="left", fontsize=11)
    fig.tight_layout(); fig.savefig("k_selection.png", dpi=200); plt.close(fig)


def gene_frequency_plot(enr, summary):
    freq = enr.pivot(index="gene", columns="cluster", values="freq_in_cluster")
    sig = enr.assign(s=(enr["q_value"] < 0.05) & (enr["direction"] == "enriched")) \
        .pivot(index="gene", columns="cluster", values="s")
    freq = freq.loc[freq.max(axis=1).sort_values(ascending=False).index]
    sig = sig.loc[freq.index, freq.columns]
    sizes = summary.set_index("cluster")["n"]
    fig, ax = plt.subplots(figsize=(1.4 + 1.0 * freq.shape[1], 0.32 * len(freq) + 1.6))
    cmap = LinearSegmentedColormap.from_list("blue", ["#f4f8fd"] + BLUE_RAMP)
    im = ax.imshow(freq.values, aspect="auto", cmap=cmap, vmin=0, vmax=1)
    for i in range(freq.shape[0]):
        for j in range(freq.shape[1]):
            v = freq.values[i, j]
            ax.text(j, i, f"{100 * v:.0f}{'*' if sig.values[i, j] else ''}", ha="center", va="center",
                    fontsize=7.5, color="#ffffff" if v > 0.5 else INK)
    ax.set_xticks(range(freq.shape[1]))
    ax.set_xticklabels([f"{c}\nn={sizes.get(c, '')}" for c in freq.columns], color=INK2)
    ax.set_yticks(range(len(freq))); ax.set_yticklabels(freq.index, fontstyle="italic")
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    cb = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.02)
    cb.set_label("Fraction of tumours mutated"); cb.outline.set_visible(False)
    ax.set_title("Mutation frequency by cluster (%)\n* enriched vs other clusters, Fisher BH q < 0.05",
                 loc="left", fontsize=10)
    fig.tight_layout(); fig.savefig("cluster_gene_frequency.png", dpi=200); plt.close(fig)


def characteristics_plot(summary):
    s = summary.sort_values("cluster")
    colors = cluster_colors(s["cluster"])
    panels = [("pct_hypermutated", "% hypermutated", 100), ("pct_READ", "% rectal (READ)", 100),
              ("stability", "Bootstrap stability\n(mean Jaccard)", 1)]
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.4), sharey=True)
    y = np.arange(len(s))
    for ax, (col, title, xmax) in zip(axes, panels):
        vals = s[col].values.astype(float)
        ax.barh(y, np.nan_to_num(vals), color=[colors[c] for c in s["cluster"]], height=0.6)
        for yi, v in zip(y, vals):
            ax.text((0 if np.isnan(v) else v) + xmax * 0.02, yi, "n/a" if np.isnan(v) else
                    (f"{v:.2f}" if xmax == 1 else f"{v:.0f}"), va="center", fontsize=8, color=INK2,
                    bbox=dict(facecolor=SURFACE, edgecolor="none", pad=1))
        ax.set_xlim(0, xmax * 1.15); ax.set_title(title, loc="left", fontsize=9.5)
        ax.grid(axis="x", color=GRID, lw=0.5); ax.set_axisbelow(True)
        if col == "stability":
            ax.axvline(0.75, color=MUTED, lw=1, ls="--")
            ax.axvline(0.5, color=MUTED, lw=1, ls=":")
    axes[0].set_yticks(y)
    axes[0].set_yticklabels([f"{c}  (n={n})" for c, n in zip(s["cluster"], s["n"])])
    axes[0].invert_yaxis()
    fig.suptitle("Cluster characteristics", x=0.01, ha="left", fontsize=11)
    fig.text(0.01, 0.01, "Stability: dashed 0.75 = stable, dotted 0.5 = dissolves under resampling. "
             "C0 = tumours with no mutation in the panel genes (not clustered).", fontsize=7.5, color=MUTED)
    fig.tight_layout(rect=(0, 0.05, 1, 0.95)); fig.savefig("cluster_characteristics.png", dpi=200); plt.close(fig)


def main():
    a = parse_args()
    types = pd.read_csv(a.types, sep="\t", index_col=0, keep_default_na=False)
    ann = pd.read_csv(a.annotation, sep="\t", index_col=0, keep_default_na=False)
    ann["tmb"] = ann["tmb"].astype(float)
    assign = pd.read_csv(a.assignments, sep="\t", index_col=0)["cluster"]
    Z = pd.read_csv(a.linkage, sep="\t").values.astype(float)
    linked = pd.read_csv(a.linkage_samples, sep="\t")["sample"].tolist()
    ksel = pd.read_csv(a.k_selection, sep="\t")
    enr = pd.read_csv(a.enrichment, sep="\t")
    summary = pd.read_csv(a.summary, sep="\t")
    genes_meta = pd.read_csv(a.genes, sep="\t")
    with open(a.cluster_info) as f:
        info = json.load(f)

    oncoprint(types, ann, assign, Z, linked, genes_meta, info, a.hypermut_tmb)
    k_selection_plot(ksel)
    gene_frequency_plot(enr, summary)
    characteristics_plot(summary)


if __name__ == "__main__":
    main()
