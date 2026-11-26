#!/usr/bin/env python3
"""DIFFERENTIAL_EXPRESSION: compare the mutation clusters with DESeq2 (PyDESeq2).

  * One negative-binomial GLM is fitted to all clustered tumours with design
    ~ batch + cluster (batch = TCGA tissue source site), then every pair of clusters is
    tested with a Wald test (C2 vs C1, C3 vs C1, C3 vs C2 for k = 3).
  * Significance: BH-adjusted p < --alpha and |log2 fold change| >= --min-lfc.
  * For the figures, log2(normalised count + 1) is written for the top DE genes, and an
    expression PCA is computed on the most variable genes so the clusters can be seen in
    expression space.

Outputs: de_<A>_vs_<B>.tsv per contrast, de_summary.tsv, top_genes_log2norm.tsv,
expression_pca.tsv, expression_pca_variance.tsv
"""
import argparse
import itertools

import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from sklearn.decomposition import PCA


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--counts", required=True)
    p.add_argument("--metadata", required=True)
    p.add_argument("--genes", required=True)
    p.add_argument("--alpha", type=float, default=0.05)
    p.add_argument("--min-lfc", type=float, default=1.0)
    p.add_argument("--covariates", default="batch", help="comma-separated metadata columns, or 'none'")
    p.add_argument("--n-top", type=int, default=25, help="top DE genes per contrast for the heatmap")
    p.add_argument("--n-pca-genes", type=int, default=2000)
    p.add_argument("--cpus", type=int, default=1)
    return p.parse_args()


def main():
    a = parse_args()
    counts = pd.read_csv(a.counts, sep="\t", index_col=0)
    meta = pd.read_csv(a.metadata, sep="\t", index_col=0).loc[counts.columns]
    symbols = pd.read_csv(a.genes, sep="\t", index_col=0)["symbol"]

    covariates = [] if a.covariates == "none" else [c for c in a.covariates.split(",") if c]
    clusters = sorted(meta["cluster"].unique())
    design = "~" + " + ".join(covariates + ["cluster"])
    print(f"[de] {counts.shape[1]} tumours, {counts.shape[0]} genes, design {design}", flush=True)

    dds = DeseqDataSet(counts=counts.T, metadata=meta[covariates + ["cluster"]].astype(str),
                       design=design, n_cpus=a.cpus, quiet=True)
    dds.deseq2()

    summary, top = [], []
    for ref, test in itertools.combinations(clusters, 2):  # e.g. (C1, C2) -> C2 vs C1
        name = f"{test}_vs_{ref}"
        stats = DeseqStats(dds, contrast=["cluster", test, ref], alpha=a.alpha, n_cpus=a.cpus, quiet=True)
        stats.summary()
        res = stats.results_df.copy()
        res.insert(0, "symbol", symbols.reindex(res.index))
        res["significant"] = (res["padj"] < a.alpha) & (res["log2FoldChange"].abs() >= a.min_lfc)
        res = res.sort_values(["padj", "pvalue"])
        res.index.name = "gene_id"
        res.to_csv(f"de_{name}.tsv", sep="\t")

        sig = res[res["significant"]]
        summary.append({"contrast": name, "genes_tested": int(res["padj"].notna().sum()),
                        "padj_below_alpha": int((res["padj"] < a.alpha).sum()),
                        "significant": len(sig),
                        f"up_in_{test}": int((sig["log2FoldChange"] > 0).sum()),
                        f"down_in_{test}": int((sig["log2FoldChange"] < 0).sum())})
        top.extend(sig.index[: a.n_top])
        print(f"[de] {name}: {len(sig)} significant genes", flush=True)

    summary = pd.DataFrame(summary)
    # keep the up/down columns readable in one table
    up = summary.filter(like="up_in_").bfill(axis=1).iloc[:, 0]
    down = summary.filter(like="down_in_").bfill(axis=1).iloc[:, 0]
    summary = summary[["contrast", "genes_tested", "padj_below_alpha", "significant"]].assign(
        up=up.astype(int), down=down.astype(int))
    summary.to_csv("de_summary.tsv", sep="\t", index=False)

    # log2 normalised counts for plotting
    norm = pd.DataFrame(dds.layers["normed_counts"], index=counts.columns, columns=counts.index).T
    log2n = np.log2(norm + 1)
    top = list(dict.fromkeys(top))
    out = log2n.loc[top]
    out.insert(0, "symbol", symbols.reindex(top))
    out.to_csv("top_genes_log2norm.tsv", sep="\t")

    hv = log2n.var(axis=1).sort_values(ascending=False).index[: a.n_pca_genes]
    X = log2n.loc[hv].T
    pca = PCA(n_components=10, random_state=0).fit(X - X.mean())
    pcs = pd.DataFrame(pca.transform(X - X.mean())[:, :4], index=X.index, columns=["PC1", "PC2", "PC3", "PC4"])
    pcs.join(meta[["cluster", "batch", "side", "hypermutated"]]).to_csv("expression_pca.tsv", sep="\t")
    pd.DataFrame({"pc": [f"PC{i + 1}" for i in range(10)],
                  "variance_explained": pca.explained_variance_ratio_.round(4)}) \
        .to_csv("expression_pca_variance.tsv", sep="\t", index=False)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
