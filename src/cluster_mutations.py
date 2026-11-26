#!/usr/bin/env python3
"""CLUSTER_MUTATIONS: hierarchical clustering of tumours on their driver-mutation profiles.

  * Distance: Jaccard on the binary matrix (shared mutated genes / genes mutated in either tumour).
    Tumours with no mutation in any kept gene have no defined Jaccard distance; they are set aside
    as group "C0_no_driver" and not clustered.
  * Linkage: Ward by default (average / complete available).
  * k: chosen as the k in [k-min, k-max] with the highest mean silhouette, among cuts whose
    smallest cluster has >= --min-cluster-size tumours (so each cluster is usable for expression
    analysis). --k overrides the choice. Assignments for every k are written as well.
  * Stability: each cluster's mean Jaccard recovery over --n-boot bootstrap resamples of tumours
    (Hennig 2007, as in fpc::clusterboot): >0.75 stable, 0.6-0.75 weak, <0.5 not reproducible.
  * Characterisation: Fisher exact test per gene (in cluster vs other clustered tumours, BH-adjusted),
    and chi-square tests of cluster vs hypermutation, tumour side, project and stage.

Outputs: cluster_assignments.tsv, cluster_assignments_all_k.tsv, k_selection.tsv, linkage.tsv,
cluster_stability.tsv, cluster_gene_enrichment.tsv, cluster_annotation_tests.tsv,
cluster_summary.tsv, cluster_summary.json
"""
import argparse
import json

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
from scipy.stats import chi2_contingency, fisher_exact
from sklearn.metrics import silhouette_score
from statsmodels.stats.multitest import multipletests

NO_DRIVER = "C0_no_driver"


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--matrix", required=True)
    p.add_argument("--annotation", required=True)
    p.add_argument("--linkage", default="ward", choices=["ward", "average", "complete"])
    p.add_argument("--k-min", type=int, default=2)
    p.add_argument("--k-max", type=int, default=8)
    p.add_argument("--k", type=int, default=0, help="force this k (0 = choose by silhouette)")
    p.add_argument("--min-cluster-size", type=int, default=20)
    p.add_argument("--n-boot", type=int, default=200)
    p.add_argument("--seed", type=int, default=501)
    return p.parse_args()


def jaccard_distance(X):
    X = X.astype(np.int32)
    inter = X @ X.T
    n = X.sum(axis=1)
    union = n[:, None] + n[None, :] - inter
    D = 1.0 - inter / union
    np.fill_diagonal(D, 0.0)
    return D


def cut(D, method, k):
    Z = linkage(squareform(D, checks=False), method=method)
    return Z, fcluster(Z, k, criterion="maxclust")


def relabel_by_size(labels):
    """Renumber clusters 1..k from largest to smallest so names are stable and readable."""
    order = pd.Series(labels).value_counts().index
    mapping = {old: new for new, old in enumerate(order, start=1)}
    return np.array([mapping[x] for x in labels])


def bootstrap_stability(D, labels, method, k, n_boot, rng):
    n = len(labels)
    clusters = np.unique(labels)
    scores = {c: [] for c in clusters}
    for _ in range(n_boot):
        idx = np.unique(rng.integers(0, n, n))
        _, boot = cut(D[np.ix_(idx, idx)], method, k)
        orig = labels[idx]
        for c in clusters:
            a = orig == c
            if not a.any():
                continue
            best = max(np.sum(a & (boot == b)) / np.sum(a | (boot == b)) for b in np.unique(boot))
            scores[c].append(best)
    return {c: float(np.mean(v)) if v else float("nan") for c, v in scores.items()}


def main():
    a = parse_args()
    rng = np.random.default_rng(a.seed)
    M = pd.read_csv(a.matrix, sep="\t", index_col=0)
    # "NA" stays a string (unknown) so it can be excluded explicitly below
    ann = pd.read_csv(a.annotation, sep="\t", index_col=0, keep_default_na=False).loc[M.index]

    has_driver = M.sum(axis=1) > 0
    Mc = M[has_driver]
    D = jaccard_distance(Mc.values)

    # k selection
    rows, all_k = [], pd.DataFrame(index=M.index)
    for k in range(a.k_min, a.k_max + 1):
        _, lab = cut(D, a.linkage, k)
        lab = relabel_by_size(lab)
        sizes = np.bincount(lab)[1:]
        sil = silhouette_score(D, lab, metric="precomputed") if len(sizes) > 1 else np.nan
        rows.append({"k": k, "mean_silhouette": round(float(sil), 4), "min_cluster_size": int(sizes.min()),
                     "cluster_sizes": ",".join(map(str, sizes)),
                     "eligible": bool(sizes.min() >= a.min_cluster_size and len(sizes) == k)})
        all_k[f"k{k}"] = pd.Series([f"C{x}" for x in lab], index=Mc.index)
    all_k = all_k.fillna(NO_DRIVER)
    ksel = pd.DataFrame(rows)

    if a.k:
        k = a.k
        reason = "set by --k"
    elif ksel["eligible"].any():
        k = int(ksel.loc[ksel[ksel["eligible"]]["mean_silhouette"].idxmax(), "k"])
        reason = f"highest mean silhouette with every cluster >= {a.min_cluster_size} tumours"
    else:
        k = int(ksel.loc[ksel["mean_silhouette"].idxmax(), "k"])
        reason = "highest mean silhouette (no k met the minimum cluster size)"
    ksel["chosen"] = ksel["k"] == k

    Z, lab = cut(D, a.linkage, k)
    lab = relabel_by_size(lab)
    stability = bootstrap_stability(D, lab, a.linkage, k, a.n_boot, rng)

    assign = pd.Series(NO_DRIVER, index=M.index, name="cluster")
    assign[Mc.index] = [f"C{x}" for x in lab]

    # gene enrichment: each cluster vs the other clustered tumours
    enr = []
    for c in sorted(set(assign[Mc.index])):
        inc = (assign[Mc.index] == c).values
        for g in M.columns:
            x = Mc[g].values.astype(bool)
            table = [[np.sum(x & inc), np.sum(~x & inc)], [np.sum(x & ~inc), np.sum(~x & ~inc)]]
            odds, p = fisher_exact(table)
            enr.append({"cluster": c, "gene": g, "freq_in_cluster": table[0][0] / inc.sum(),
                        "freq_in_others": table[1][0] / (~inc).sum(), "odds_ratio": odds, "p_value": p})
    enr = pd.DataFrame(enr)
    enr["q_value"] = multipletests(enr["p_value"], method="fdr_bh")[1]
    enr["direction"] = np.where(enr["freq_in_cluster"] > enr["freq_in_others"], "enriched", "depleted")
    enr = enr.round({"freq_in_cluster": 4, "freq_in_others": 4, "odds_ratio": 3})

    # clinical / genomic annotation vs cluster (clustered tumours, known values only)
    tests = []
    for col in ["hypermutated", "side", "project", "stage"]:
        sub = ann.loc[Mc.index, col].astype(str)
        known = sub != "NA"
        tab = pd.crosstab(assign[Mc.index][known], sub[known])
        chi2, p, dof, _ = chi2_contingency(tab)
        tests.append({"variable": col, "n": int(known.sum()), "chi2": round(chi2, 2), "dof": dof, "p_value": p})
    tests = pd.DataFrame(tests)

    # per-cluster summary
    full = ann.join(assign)
    summ = []
    for c, grp in full.groupby("cluster"):
        top = enr[(enr["cluster"] == c) & (enr["direction"] == "enriched") & (enr["q_value"] < 0.05)]
        top = top.sort_values("q_value")["gene"].head(5).tolist()
        summ.append({"cluster": c, "n": len(grp),
                     "median_tmb": round(grp["tmb"].median(), 2),
                     "pct_hypermutated": round(100 * (grp["hypermutated"] == "Hypermutated").mean(), 1),
                     "pct_right_colon": round(100 * (grp.loc[grp["side"] != "NA", "side"] == "Right colon").mean(), 1),
                     "pct_READ": round(100 * (grp["project"] == "READ").mean(), 1),
                     "stability": round(stability.get(int(c[1:]), np.nan), 3) if c != NO_DRIVER else np.nan,
                     "enriched_genes_q05": ";".join(top)})
    summ = pd.DataFrame(summ)

    assign.to_frame().to_csv("cluster_assignments.tsv", sep="\t")
    all_k.to_csv("cluster_assignments_all_k.tsv", sep="\t")
    ksel.to_csv("k_selection.tsv", sep="\t", index=False)
    pd.DataFrame(Z, columns=["left", "right", "height", "size"]).to_csv("linkage.tsv", sep="\t", index=False)
    pd.Series(Mc.index, name="sample").to_csv("linkage_samples.tsv", sep="\t", index=False)
    pd.DataFrame({"cluster": [f"C{c}" for c in stability], "mean_jaccard_recovery": stability.values()}) \
        .to_csv("cluster_stability.tsv", sep="\t", index=False)
    enr.to_csv("cluster_gene_enrichment.tsv", sep="\t", index=False)
    tests.to_csv("cluster_annotation_tests.tsv", sep="\t", index=False)
    summ.to_csv("cluster_summary.tsv", sep="\t", index=False)
    info = {"linkage": a.linkage, "distance": "jaccard", "k": k, "k_reason": reason,
            "n_clustered": int(has_driver.sum()), "n_no_driver": int((~has_driver).sum()),
            "n_boot": a.n_boot, "seed": a.seed}
    with open("cluster_summary.json", "w") as f:
        json.dump(info, f, indent=2)
    print(json.dumps(info, indent=2))
    print(summ.to_string(index=False))


if __name__ == "__main__":
    main()
