#!/usr/bin/env python3
"""PREPARE_COUNTS: build the raw count matrix and sample table for differential expression.

  * UCSC Xena stores GDC STAR gene counts as log2(count + 1); they are converted back to
    integer counts, which is what DESeq2 models.
  * Only tumours in the mutation clusters are kept (the no-driver group C0 is dropped).
  * Genes are kept if they have >= --min-count reads in at least as many tumours as the
    smallest cluster (the DESeq2 vignette's pre-filter), so a gene expressed in only one
    cluster can still be tested.
  * The tissue source site (TCGA barcode characters 6-7, e.g. TCGA-AA-xxxx -> AA) is a strong
    batch effect in TCGA; sites with fewer than --min-site tumours are pooled as "other".

Outputs: counts.tsv.gz (genes x tumours), sample_metadata.tsv, gene_annotation.tsv,
prepare_counts_summary.json
"""
import argparse
import json

import numpy as np
import pandas as pd

NO_DRIVER = "C0_no_driver"


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--expression", nargs="+", required=True)
    p.add_argument("--probemap", required=True)
    p.add_argument("--assignments", required=True)
    p.add_argument("--annotation", required=True)
    p.add_argument("--min-count", type=int, default=10)
    p.add_argument("--min-site", type=int, default=10)
    return p.parse_args()


def main():
    a = parse_args()
    assign = pd.read_csv(a.assignments, sep="\t", index_col=0)["cluster"]
    assign = assign[assign != NO_DRIVER]
    ann = pd.read_csv(a.annotation, sep="\t", index_col=0, keep_default_na=False)

    parts = []
    for path in a.expression:
        header = pd.read_csv(path, sep="\t", nrows=0).columns
        wanted = ["Ensembl_ID"] + [s for s in header[1:] if s in assign.index]
        parts.append(pd.read_csv(path, sep="\t", usecols=wanted, index_col=0))
    log2 = pd.concat(parts, axis=1)
    log2 = log2.loc[~log2.index.str.endswith("_PAR_Y"), assign.index]  # duplicated chrY copies
    counts = np.rint(np.exp2(log2) - 1).clip(lower=0).astype(np.int64)

    smallest = int(assign.value_counts().min())
    expressed = (counts >= a.min_count).sum(axis=1) >= smallest
    counts = counts[expressed]

    site = pd.Series(assign.index.str[5:7], index=assign.index)
    big = site.value_counts()
    site = site.where(site.map(big) >= a.min_site, "other")

    meta = pd.DataFrame({
        "cluster": assign,
        "batch": site,
        "project": ann["project"].reindex(assign.index),
        "side": ann["side"].reindex(assign.index),
        "hypermutated": ann["hypermutated"].reindex(assign.index),
        "tmb": ann["tmb"].reindex(assign.index),
        "library_size": counts.sum(axis=0),
    })
    meta.index.name = "sample"

    probemap = pd.read_csv(a.probemap, sep="\t", index_col=0)
    genes = probemap.reindex(counts.index)[["gene", "chrom"]].rename(columns={"gene": "symbol"})
    genes.index.name = "gene_id"

    counts.index.name = "gene_id"
    counts.to_csv("counts.tsv.gz", sep="\t", compression="gzip")
    meta.to_csv("sample_metadata.tsv", sep="\t")
    genes.to_csv("gene_annotation.tsv", sep="\t")
    summary = {
        "tumours": int(counts.shape[1]),
        "tumours_per_cluster": assign.value_counts().sort_index().to_dict(),
        "genes_in_file": int(len(log2)),
        "genes_kept": int(expressed.sum()),
        "filter": f">= {a.min_count} counts in >= {smallest} tumours (smallest cluster)",
        "batch_levels": int(meta["batch"].nunique()),
        "median_library_size": int(meta["library_size"].median()),
    }
    with open("prepare_counts_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
