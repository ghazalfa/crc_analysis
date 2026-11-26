#!/usr/bin/env python3
"""MAKE_REPORT: fill templates/report_template.md with this run's numbers and tables -> report.md.

All wording lives in the template. {{name}} marks a value or table supplied here (every pipeline
parameter is available under its own name). Every figure is copied to figures/ next to report.md,
so the figures/<name>.png links in the template resolve.
"""
import argparse
import json
import re
import shutil
from importlib.metadata import version
from pathlib import Path

import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for n in ["build-summary", "cluster-info", "cluster-summary", "k-selection", "annotation-tests",
              "counts-summary", "de-summary"]:
        p.add_argument(f"--{n}", required=True)
    p.add_argument("--template", required=True)
    p.add_argument("--figures", nargs="+", required=True)
    p.add_argument("--settings", nargs="+", required=True, help="key=value pipeline parameters")
    return p.parse_args()


def load(path):
    return json.loads(Path(path).read_text())


def table(df):
    """DataFrame as a Markdown table ('|' in cells is escaped)."""
    body = df.astype(object).where(df.notna(), "n/a")
    cell = lambda x: str(x).replace("|", "\\|")
    lines = [list(body.columns), ["---"] * body.shape[1], *body.values.tolist()]
    return "\n".join("| " + " | ".join(cell(x) for x in row) + " |" for row in lines)


def fill(template, values):
    """Replace every {{name}}; an unknown name stops the run instead of reaching the report."""
    names = set(re.findall(r"\{\{\s*(\w+)\s*\}\}", template))
    unknown = sorted(names - set(values))
    if unknown:
        raise SystemExit(f"Unknown placeholder(s) in the report template: {', '.join(unknown)}\n"
                         f"Available: {', '.join(sorted(values))}")
    return re.sub(r"\{\{\s*(\w+)\s*\}\}", lambda m: str(values[m.group(1)]), template)


def check_figures(text, figures):
    """Stop the run if the template links to a figure that was not made."""
    made = {Path(f).name for f in figures}
    missing = sorted(set(re.findall(r"\(figures/([^)\s]+\.png)\)", text)) - made)
    if missing:
        raise SystemExit(f"The report template links to figure(s) that were not made: {', '.join(missing)}")


def main():
    a = parse_args()
    s = dict(x.split("=", 1) for x in a.settings)
    build, info, counts = load(a.build_summary), load(a.cluster_info), load(a.counts_summary)
    clusters = pd.read_csv(a.cluster_summary, sep="\t")
    clusters = clusters.drop(columns=["median_tmb", "pct_right_colon", "enriched_genes_q05"])
    ksel = pd.read_csv(a.k_selection, sep="\t")
    tests = pd.read_csv(a.annotation_tests, sep="\t")
    tests["p_value"] = tests["p_value"].map("{:.2g}".format)
    tests["variable"] = tests["variable"].replace({"project": "READ or COAD"})
    de = pd.read_csv(a.de_summary, sep="\t")

    best = ksel.loc[ksel["mean_silhouette"].idxmax()]
    chosen = ksel[ksel["chosen"]].iloc[0]
    k_text = (f"k = {info['k']} in this run.") 
             
    sample_flow = pd.DataFrame([
        ("Tumour samples in the mutation files", build["samples_in_mutation_files"]),
        ("Primary tumours", build["primary_tumour_samples"]),
        ("One sample per patient", build["after_one_aliquot_per_patient"]),
        ("After the RNA-seq requirement", build["after_requiring_expression"]),
        ("After the hypermutation filter", build["after_hypermutation_filter"]),
        ("Clustered (at least one kept driver gene mutated)", info["n_clustered"]),
        ("Set aside with no kept driver gene mutated (C0)", info["n_no_driver"]),
    ], columns=["Step", "Tumours"])
    versions = pd.DataFrame([(p, version(p)) for p in ["numpy", "pandas", "scipy", "scikit-learn", "statsmodels",
                                                       "pydeseq2", "matplotlib"]], columns=["Package", "Version"])
    cov = "" if s["de_covariates"] == "none" else s["de_covariates"].replace(",", " + ") + " + "

    values = {
        **s,
        **{key: build[key] for key in ["panel_genes", "panel_genes_kept", "hypermutated_tumours"]},
        "k_text": k_text,
        "genes_kept": f"{counts['genes_kept']:,}",
        "genes_in_file": f"{counts['genes_in_file']:,}",
        "count_filter": counts["filter"],
        "tumours": counts["tumours"],
        "batch_levels": counts["batch_levels"],
        "design": f"~ {cov}cluster",
        "table_sample_flow": table(sample_flow),
        "table_clusters": table(clusters),
        "table_tests": table(tests),
        "table_de": table(de),
        "table_software": table(versions),
    }
    text = fill(Path(a.template).read_text(), values)
    check_figures(text, a.figures)
    Path("report.md").write_text(text, encoding="utf-8")
    Path("figures").mkdir(exist_ok=True)
    for f in a.figures:
        shutil.copy(f, Path("figures") / Path(f).name)


if __name__ == "__main__":
    main()
