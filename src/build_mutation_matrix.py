#!/usr/bin/env python3
"""BUILD_MUTATION_MATRIX: turn TCGA masked somatic mutation calls into a sample x driver-gene matrix.

Steps
  1. Pool COAD + READ mutation calls (GDC MuTect2/MuSE/VarScan2 consensus, masked MAF via UCSC Xena).
  2. Keep primary tumours (sample type 01), one aliquot per patient, and, by default, only
     tumours that also have RNA-seq so every clustered sample can be used for expression analysis.
  3. Keep protein-altering mutations (missense, truncating, splice-site, in-frame indels).
  4. Compute tumour mutational burden (TMB, protein-altering mutations per Mb of exome) and flag
     hypermutated tumours (MSI / POLE-like); they are removed if --exclude-hypermutated is true.
  5. Restrict to the CRC driver panel and keep genes mutated in >= --min-freq of tumours.

Outputs
  mutation_matrix.tsv     binary 0/1, samples x genes (input to clustering)
  mutation_types.tsv      mutation class per sample x gene (for the oncoprint colours)
  sample_annotation.tsv   per-sample project, TMB, hypermutation flag, stage, tumour side
  gene_summary.tsv        per panel gene: mutation frequency and whether it was kept
  build_summary.json      sample / gene counts at each filter
"""
import argparse
import gzip
import json
import re

import numpy as np
import pandas as pd

# Sequence Ontology terms that change the protein, grouped into oncoprint classes.
CLASS_OF_TERM = {
    "frameshift_variant": "Truncating",
    "stop_gained": "Truncating",
    "splice_acceptor_variant": "Truncating",
    "splice_donor_variant": "Truncating",
    "start_lost": "Truncating",
    "missense_variant": "Missense",
    "stop_lost": "Missense",
    "protein_altering_variant": "Missense",
    "inframe_deletion": "In-frame indel",
    "inframe_insertion": "In-frame indel",
}
CLASS_PRIORITY = ["Truncating", "Missense", "In-frame indel"]  # most to least severe

RIGHT_COLON = {"Cecum", "Ascending colon", "Hepatic flexure of colon", "Transverse colon"}
LEFT_COLON = {"Splenic flexure of colon", "Descending colon", "Sigmoid colon"}
RECTUM = {"Rectosigmoid junction", "Rectum, NOS"}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--mutations", nargs="+", required=True)
    p.add_argument("--expression", nargs="+", required=True, help="only the header (sample IDs) is read")
    p.add_argument("--clinical", nargs="+", required=True)
    p.add_argument("--genes", required=True, help="TSV with a gene column")
    p.add_argument("--min-freq", type=float, default=0.05, help="minimum fraction of tumours mutated")
    p.add_argument("--exome-mb", type=float, default=38.0, help="exome size used for TMB")
    p.add_argument("--hypermut-tmb", type=float, default=10.0, help="TMB (mut/Mb) at or above which a tumour is hypermutated")
    p.add_argument("--exclude-hypermutated", choices=["true", "false"], default="false")
    p.add_argument("--require-expression", choices=["true", "false"], default="true")
    return p.parse_args()


def project_of(path):
    m = re.search(r"TCGA-(COAD|READ)", path)
    return m.group(1) if m else "NA"


def expression_samples(paths):
    samples = set()
    for path in paths:
        with gzip.open(path, "rt") as f:
            samples.update(f.readline().rstrip("\n").split("\t")[1:])
    return samples


def classify(effect):
    """Most severe protein-altering class in a ';'-joined effect string, or None."""
    classes = {CLASS_OF_TERM[t] for t in str(effect).split(";") if t in CLASS_OF_TERM}
    for c in CLASS_PRIORITY:
        if c in classes:
            return c
    return None


def pick_one_aliquot_per_patient(samples, with_expression):
    """TCGA sample IDs look like TCGA-XX-XXXX-01A. Prefer aliquots with RNA-seq, then vial A."""
    table = pd.DataFrame({"sample": sorted(samples)})
    table["patient"] = table["sample"].str[:12]
    table["has_expr"] = table["sample"].isin(with_expression)
    table = table.sort_values(["patient", "has_expr", "sample"], ascending=[True, False, True])
    return set(table.drop_duplicates("patient")["sample"])


def collapse_stage(stage):
    m = re.match(r"Stage (IV|III|II|I)", str(stage))
    return m.group(1) if m else "NA"


def side_of(site):
    if site in RIGHT_COLON:
        return "Right colon"
    if site in LEFT_COLON:
        return "Left colon"
    if site in RECTUM:
        return "Rectum"
    return "NA"


def main():
    a = parse_args()
    counts = {}

    # 1. pool mutation calls
    muts = pd.concat(
        [pd.read_csv(p, sep="\t", usecols=["sample", "gene", "effect"]).assign(project=project_of(p))
         for p in a.mutations],
        ignore_index=True,
    )
    counts["samples_in_mutation_files"] = int(muts["sample"].nunique())

    # 2. sample selection
    with_expr = expression_samples(a.expression)
    primary = {s for s in muts["sample"].unique() if s[13:15] == "01"}
    counts["primary_tumour_samples"] = len(primary)
    keep = pick_one_aliquot_per_patient(primary, with_expr)
    counts["after_one_aliquot_per_patient"] = len(keep)
    if a.require_expression == "true":
        keep &= with_expr
    counts["after_requiring_expression"] = len(keep)
    muts = muts[muts["sample"].isin(keep)]
    samples = sorted(keep)

    # 3. protein-altering mutations only
    muts = muts.assign(mclass=muts["effect"].map(classify)).dropna(subset=["mclass"])
    counts["protein_altering_mutations"] = int(len(muts))

    # 4. TMB per tumour
    project = muts.drop_duplicates("sample").set_index("sample")["project"]
    n_nonsyn = muts.groupby("sample").size().reindex(samples, fill_value=0)
    ann = pd.DataFrame(index=pd.Index(samples, name="sample"))
    ann["project"] = project.reindex(samples).fillna("NA")
    ann["n_protein_altering"] = n_nonsyn
    ann["tmb"] = (n_nonsyn / a.exome_mb).round(3)
    ann["hypermutated"] = np.where(ann["tmb"] >= a.hypermut_tmb, "Hypermutated", "Non-hypermutated")

    clin = pd.concat([pd.read_csv(p, sep="\t", low_memory=False) for p in a.clinical], ignore_index=True)
    clin = clin.drop_duplicates("sample").set_index("sample")
    ann["stage"] = clin["ajcc_pathologic_stage.diagnoses"].reindex(samples).map(collapse_stage)
    ann["tumour_site"] = clin["tissue_or_organ_of_origin.diagnoses"].reindex(samples).fillna("NA")
    ann["side"] = ann["tumour_site"].map(side_of)
    ann["has_expression"] = ann.index.isin(with_expr)

    counts["hypermutated_tumours"] = int((ann["hypermutated"] == "Hypermutated").sum())
    if a.exclude_hypermutated == "true":
        ann = ann[ann["hypermutated"] != "Hypermutated"]
        samples = ann.index.tolist()
        muts = muts[muts["sample"].isin(samples)]
    counts["after_hypermutation_filter"] = len(samples)

    # 5. driver panel x frequency filter (frequency is among the remaining tumours)
    panel = pd.read_csv(a.genes, sep="\t")
    pm = muts[muts["gene"].isin(panel["gene"])]
    # one class per sample x gene; two or more distinct mutations in the same gene = "Multi-hit"
    per_cell = pm.groupby(["sample", "gene"]).agg(n=("mclass", "size"), mclass=("mclass", "first"))
    worst = pm.assign(rank=pm["mclass"].map(CLASS_PRIORITY.index)).groupby(["sample", "gene"])["rank"].min()
    per_cell["mclass"] = worst.map(lambda r: CLASS_PRIORITY[r])
    per_cell.loc[per_cell["n"] > 1, "mclass"] = "Multi-hit"

    types = per_cell["mclass"].unstack("gene").reindex(index=samples, columns=panel["gene"])
    binary = types.notna().astype(int)
    freq = binary.mean()

    gene_summary = panel.set_index("gene")
    gene_summary["n_mutated"] = binary.sum()
    gene_summary["freq"] = freq.round(4)
    gene_summary["kept"] = freq >= a.min_freq
    gene_summary = gene_summary.sort_values("freq", ascending=False)
    kept = gene_summary.index[gene_summary["kept"]].tolist()
    counts["panel_genes"] = int(len(panel))
    counts["panel_genes_kept"] = len(kept)
    counts["tumours_with_no_kept_driver_mutation"] = int((binary[kept].sum(axis=1) == 0).sum())
    counts["min_freq"] = a.min_freq

    binary[kept].to_csv("mutation_matrix.tsv", sep="\t")
    types[kept].fillna("").to_csv("mutation_types.tsv", sep="\t")
    ann.to_csv("sample_annotation.tsv", sep="\t")
    gene_summary.to_csv("gene_summary.tsv", sep="\t")
    with open("build_summary.json", "w") as f:
        json.dump(counts, f, indent=2)
    print(json.dumps(counts, indent=2))


if __name__ == "__main__":
    main()
