#!/usr/bin/env python3
"""MAKE_REPORT: fill report/report_template.md with this run's numbers and tables -> report.html.

All wording lives in the template. {{name}} marks a value or table supplied here (every pipeline
parameter is available under its own name). PNG figures are embedded, so report.html is a single
self-contained file; every figure is also copied to figures/.
"""
import argparse
import base64
import json
import re
import shutil
from importlib.metadata import version
from pathlib import Path

import markdown
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


STYLE = """
body{font:16px/1.55 system-ui,-apple-system,"Segoe UI",sans-serif;color:#0b0b0b;background:#fcfcfb;
max-width:1040px;margin:0 auto;padding:24px 16px 64px}h1{line-height:1.2}
h2{margin-top:56px;padding-top:16px;border-top:1px solid #c3c2b7}h3{margin-top:28px}
.scroll{overflow-x:auto}table{border-collapse:collapse;font-size:14px}
th,td{text-align:left;padding:6px 12px 6px 0;border-bottom:1px solid #e1e0d9;vertical-align:top}
img{max-width:100%;max-height:820px;width:auto;height:auto}
"""


def to_html(text, figures):
    """Markdown -> one HTML page with figures/<name>.png links replaced by embedded images."""
    body = markdown.markdown(text, extensions=["tables", "toc", "sane_lists"])
    body = body.replace("<table>", '<div class="scroll"><table>').replace("</table>", "</table></div>")
    by_name = {Path(f).name: f for f in figures}

    def embed(m):
        if m.group(1) not in by_name:
            raise SystemExit(f"The report template links to a figure that was not made: {m.group(1)}")
        data = base64.b64encode(Path(by_name[m.group(1)]).read_bytes()).decode()
        return f'src="data:image/png;base64,{data}"'

    body = re.sub(r'src="figures/([^"]+\.png)"', embed, body)
    return ('<!doctype html><html lang="en"><head><meta charset="utf-8">'
            '<meta name="viewport" content="width=device-width,initial-scale=1">'
            f"<title>CRC mutation clusters report</title><style>{STYLE}</style></head><body>\n{body}\n</body></html>\n")


def main():
    a = parse_args()
    s = dict(x.split("=", 1) for x in a.settings)
    build, info, counts = load(a.build_summary), load(a.cluster_info), load(a.counts_summary)
    clusters = pd.read_csv(a.cluster_summary, sep="\t")
    ksel = pd.read_csv(a.k_selection, sep="\t")
    tests = pd.read_csv(a.annotation_tests, sep="\t")
    tests["p_value"] = tests["p_value"].map("{:.2g}".format)
    de = pd.read_csv(a.de_summary, sep="\t")

    best = ksel.loc[ksel["mean_silhouette"].idxmax()]
    chosen = ksel[ksel["chosen"]].iloc[0]
    k_text = (f"k = {info['k']} was set by the user (mean silhouette {chosen['mean_silhouette']}). The highest "
              f"silhouette is at k = {best['k']} ({best['mean_silhouette']}), so silhouette alone would "
              f"{'agree' if best['k'] == info['k'] else 'not have picked this k'}."
              if s["k"] != "0" else f"k = {info['k']} was chosen automatically: {info['k_reason']}.")

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
    Path("report.html").write_text(to_html(text, a.figures), encoding="utf-8")
    Path("figures").mkdir(exist_ok=True)
    for f in a.figures:
        shutil.copy(f, Path("figures") / Path(f).name)


if __name__ == "__main__":
    main()
