# Resubtyping Colorectal Cancer

A Nextflow pipeline that groups TCGA colorectal tumours by their driver-gene mutations and then compares gene expression between the groups.

BIOF 501 final project, Ghazal Fallahpour Sichani

- [Background](#background)
- [Question and rationale](#question-and-rationale)
- [Pipeline](#pipeline)
- [Data](#data)
- [Dependencies](#dependencies)
- [How to run](#how-to-run)
- [Expected results](#expected-results)
- [Directory structure](#directory-structure)
- [Limitations](#limitations)
- [References](#references)

## Background

Colorectal cancer (CRC) is a leading cause of cancer-related death worldwide, and patients differ widely in outcome and drug response. Tumours are commonly subtyped with the consensus molecular subtypes (CMS), which are predicted from gene expression by a machine learning model [1].

| Subtype | Typical features |
|---|---|
| CMS1 | High microsatellite instability (MSI), right-sided, immune activation |
| CMS2 | Low MSI, high chromosomal instability (CIN), left-sided, WNT signalling |
| CMS3 | More frequent KRAS mutation, low CIN, mixed MSI status |
| CMS4 | Low MSI, high CIN, TGF-beta and integrin signalling |

The CMS work provides two classifiers, a random forest (RF) and a single sample predictor (SSP). In a comparison of 100 TCGA patients made for this project's presentation, the two often disagreed on the same patient, and changing the protocol also changed the calls. That comparison is motivation only; this pipeline does not reproduce it.

## Question and rationale

**Can clustering somatic mutations, then interpreting the clusters by differential expression, give more reliable CRC subtypes?**

Expression-based calls depend on normalisation and on the cohort used for centring. Whether a driver gene is mutated is a simpler yes/no measurement, so it is a reasonable alternative basis for grouping tumours. Expression is then used to interpret the groups rather than to define them.

| Design choice | Reason |
|---|---|
| Nextflow | The whole analysis is automated, resumable and repeatable on another machine. |
| Driver-gene panel [2] | In an unrestricted matrix, long genes are often mutated by chance. |
| Jaccard distance | Ignores genes that neither tumour has mutated, which suits a mostly-zero matrix. |
| Ward linkage | Gave groups large enough to compare; average linkage mostly split off single tumours. |
| Silhouette and bootstrap | Measure separation and reproducibility, which are different properties. |
| PyDESeq2 with a site covariate | Models counts directly and adjusts for TCGA collection site. |

## Pipeline

```mermaid
flowchart LR
    S[data/sources.csv] --> F[FETCH]
    F -->|mutations, clinical, expression sample IDs| B[BUILD_MUTATION_MATRIX]
    P[CRC driver panel] --> B
    B --> C[CLUSTER_MUTATIONS]
    B --> PC[PLOT_CLUSTERS]
    C --> PC
    P --> PC
    F -->|STAR counts, gene annotation| R[PREPARE_COUNTS]
    B --> R
    C --> R
    R --> D[DIFFERENTIAL_EXPRESSION]
    D --> PD[PLOT_DE]
    R --> PD
    B --> M[MAKE_REPORT]
    C --> M
    R --> M
    D --> M
    PC --> M
    PD --> M
    T[report_template.md] --> M
    M --> H[report.html]
```

| Stage | What it does | Output folder |
|---|---|---|
| `FETCH` (setup) | Downloads each file and checks its md5 | `data/raw` |
| `BUILD_MUTATION_MATRIX` | Keeps one primary tumour per patient with RNA-seq, then builds a tumour x driver-gene 0/1 matrix of protein-altering mutations | `01_mutation_matrix` |
| `CLUSTER_MUTATIONS` | Hierarchical clustering; silhouette for k = 2 to 8; 200 bootstrap resamples; tests of genes and clinical variables | `02_clustering` |
| `PLOT_CLUSTERS` | Oncoprint and cluster figures | `05_report/figures` |
| `PREPARE_COUNTS` | Converts log2(count + 1) back to integer counts and removes low-count genes | `03_counts` |
| `DIFFERENTIAL_EXPRESSION` | PyDESeq2, `~ batch + cluster`, Benjamini-Hochberg adjusted | `04_differential_expression` |
| `PLOT_DE` | Volcano plot, heatmap, PCA and gene counts | `05_report/figures` |
| `MAKE_REPORT` | Fills the template with this run's numbers and embeds the figures in one HTML file | `05_report` |

## Data

All data are public TCGA COAD (colon) and READ (rectum) files, downloaded automatically from the UCSC Xena GDC hub. The total download is about 81 MB. URLs and md5 checksums are in [data/sources.csv](data/sources.csv).

| Data | Contents | Format | Size |
|---|---|---|---|
| Somatic mutations | Whole-exome masked calls, one row per mutation per tumour | MAF-like table, `.tsv.gz` | 7.95 MB |
| Gene expression | STAR counts stored as log2(count + 1), genes x samples | Matrix, `.tsv.gz` | 70.05 MB |
| Clinical | Stage, tumour site and survival fields | Table, `.tsv.gz` | 0.17 MB |
| Gene annotation | GENCODE v36 Ensembl ID to symbol map | Text | 3.21 MB |

The driver panel of 109 genes ([data/crc_driver_genes.tsv](data/crc_driver_genes.tsv), from [2]) is included in the repository. These data suit the question because mutations and expression were measured in the same tumours.

| Selection step | Tumours |
|---|---:|
| In the mutation files | 596 |
| Primary tumours | 593 |
| One sample per patient | 582 |
| With RNA-seq | 577 |
| Clustered (at least one kept driver gene mutated) | **552** |

The 25 tumours with no kept driver mutation are set aside as `C0_no_driver`.

## Dependencies

**Core (you install these):** Nextflow 24.04 or newer (tested with 26.04.7; the install below includes Java) and conda.

**Pipeline (Nextflow installs these through conda, versions pinned in [envs/environment.yml](envs/environment.yml)):**

| Tool | Version | Used for |
|---|---|---|
| Python | 3.12.16 | All stages |
| NumPy, pandas | 2.5.4, 2.3.3 | Data handling, bootstrap resampling |
| SciPy | 1.18.1 | Hierarchical clustering, statistical tests |
| scikit-learn | 1.9.1 | Silhouette score |
| statsmodels | 0.15.0 | Multiple-testing correction |
| PyDESeq2, anndata | 0.5.4, 0.13.4 | Differential expression |
| Matplotlib | 3.11.2 | Figures |
| Markdown | 3.11.3 | HTML report |

## How to run

Requirements: macOS or Linux (amd64 or arm64), internet access and about 5 GB of free disk space. On Windows, use WSL2.

**1. Install conda**, if `conda --version` fails:

```bash
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash "Miniforge3-$(uname)-$(uname -m).sh"
```

Answer `yes` when asked to initialise conda, then open a new terminal.

**2. Install Nextflow** (this also installs Java):

```bash
conda create -n nextflow --override-channels -c conda-forge -c bioconda nextflow
conda activate nextflow
nextflow -version
```

This environment is only for Nextflow. Nextflow creates the analysis environment itself.

**3. Get the code:**

```bash
git clone <repository-url>
cd <repository-folder>
```

**4. Run:**

```bash
nextflow run main.nf
```

The first run downloads the data and builds the conda environment (about 1 GB). The analysis uses up to 6 cores and took about 3 minutes on a laptop, including the environment build but not the data download.

**5. Open the report**, `results/05_report/report.html`, in a browser.

### Options

Defaults are in [nextflow.config](nextflow.config). Override any of them on the command line:

| Option | Default | Meaning |
|---|---|---|
| `--k` | 2 | Number of clusters; 0 picks k by silhouette |
| `--min_freq` | 0.05 | Minimum mutation frequency for a panel gene |
| `--exclude_hypermutated` | false | Cluster only non-hypermutated tumours |
| `--de_alpha` | 0.05 | Adjusted p-value cutoff |
| `--de_min_lfc` | 1 | Absolute log2 fold-change cutoff |
| `--max_cpus` | 6 | Maximum cores |
| `--outdir` | results | Output folder |

Add `-resume` to reuse finished steps, for example `nextflow run main.nf --k 3 -resume`.

### Troubleshooting

- **`conda: command not found`:** open a new terminal, or run `conda init` first.
- **Checksum or download error:** delete the affected file from `data/raw` and rerun.

## Expected results

Reference outputs from the default settings are in [expected_results/](expected_results/). On the same machine, two complete runs gave identical tables and PNG figures. Across machines, small floating-point differences in p-values are possible.

| Stage | Output in `results/` | Expected value | Reference copy |
|---|---|---|---|
| `BUILD_MUTATION_MATRIX` | `01_mutation_matrix/build_summary.json` | 577 tumours, 36 of 109 panel genes kept, 84 hypermutated | [build_summary.json](expected_results/build_summary.json) |
| `CLUSTER_MUTATIONS` | `02_clustering/cluster_summary.tsv` | C1 = 464, C2 = 88, C0 = 25 | [cluster_summary.tsv](expected_results/cluster_summary.tsv) |
| `CLUSTER_MUTATIONS` | `02_clustering/cluster_stability.tsv` | C1 = 0.924, C2 = 0.786 | [cluster_stability.tsv](expected_results/cluster_stability.tsv) |
| `CLUSTER_MUTATIONS` | `02_clustering/k_selection.tsv` | Silhouette 0.2233 at k = 2 | [k_selection.tsv](expected_results/k_selection.tsv) |
| `PREPARE_COUNTS` | `03_counts/prepare_counts_summary.json` | 552 tumours, 21,901 of 60,616 genes | [prepare_counts_summary.json](expected_results/prepare_counts_summary.json) |
| `DIFFERENTIAL_EXPRESSION` | `04_differential_expression/de_summary.tsv` | 1,814 genes: 742 up and 1,072 down in C2 | [de_summary.tsv](expected_results/de_summary.tsv) |
| `MAKE_REPORT` | `05_report/report.html` | HTML report with 8 figures | [report.html](expected_results/report.html) |

### Two stable but weakly separated clusters

| Cluster | Tumours | Stability (mean Jaccard) | Hypermutated | Right colon | Enriched genes |
|---|---:|---:|---:|---:|---|
| C1 | 464 | 0.924 | 1.1% | 36.5% | APC, TP53, KRAS |
| C2 | 88 | 0.786 | 89.8% | 77.0% | RNF43, BMPR2, BRAF, ACVR2A, KMT2B |

Both clusters pass the 0.75 stability guideline [3]. The silhouette is highest at k = 2, but at 0.22 it indicates weak separation. **C2 is the hypermutated group.**

![Oncoprint of driver mutations by cluster](expected_results/figures/oncoprint.png)

![Silhouette by number of clusters](expected_results/figures/k_selection.png)

### Differential expression, C2 vs C1

1,814 of 21,887 tested genes differ (adjusted p < 0.05 and at least a two-fold change): 742 are higher and 1,072 lower in C2.

![Volcano plot](expected_results/figures/volcano_plots.png)

![Heatmap of top genes](expected_results/figures/top_de_heatmap.png)

### Supplementary analyses (not part of the Nextflow run)

These R scripts were run separately on the pipeline outputs, with R 4.6.1.

**CMS composition** ([cms/run_cms.R](cms/run_cms.R), CMSclassifier random forest, log2(CPM + 1) input). C2 is 57% CMS1; C1 is mainly CMS2 (31%) and CMS4 (28%). 22% of tumours are unclassified (top probability below 0.5).

![CMS composition by cluster](expected_results/supplementary/cms_proportions.png)

**Overall survival** ([survival/run_km.R](survival/run_km.R)). No significant difference was detected in 522 patients (log-rank p = 0.33; hazard ratio C2 vs C1 0.76, 95% CI 0.44 to 1.32; adjusted for age and stage 1.08, 0.60 to 1.94). C2 has only 15 deaths, so power is low.

![Overall survival by cluster](expected_results/supplementary/km_clusters.png)

To rerun them after the pipeline, from the project folder:

```bash
Rscript -e 'install.packages("survival", repos="https://cloud.r-project.org")'
Rscript survival/run_km.R

mkdir -p cms/rlib
Rscript -e 'install.packages(c("remotes","BiocManager"), repos="https://cloud.r-project.org"); .libPaths(c("cms/rlib", .libPaths())); BiocManager::install("org.Hs.eg.db", lib="cms/rlib", update=FALSE, ask=FALSE); install.packages("randomForest", lib="cms/rlib", repos="https://cloud.r-project.org"); remotes::install_github("dFenix7/CMSclassifier@ba02ec5c88194e34c2a2b6a0c81bc07c36c8e9fa", lib="cms/rlib", upgrade="never")'
Rscript cms/run_cms.R
```

The original `Sage-Bionetworks/CMSclassifier` repository was unavailable in November 2026. The copy above has the same commit and tree hashes as the Software Heritage archive of the original.

### Conclusion

Clustering driver-gene mutations gives a reproducible grouping, but it mostly recovers known biology (hypermutation) rather than new subtypes. Mutation data alone does not yet give a finer subtyping than CMS.

## Directory structure

```text
.
├── main.nf                      # workflow
├── nextflow.config              # parameters and conda settings
├── envs/environment.yml         # pinned software
├── src/
│   ├── fetch_file.py            # FETCH
│   ├── build_mutation_matrix.py # BUILD_MUTATION_MATRIX
│   ├── cluster_mutations.py     # CLUSTER_MUTATIONS
│   ├── plot_clusters.py         # PLOT_CLUSTERS
│   ├── prepare_counts.py        # PREPARE_COUNTS
│   ├── differential_expression.py # DIFFERENTIAL_EXPRESSION
│   ├── plot_de.py               # PLOT_DE
│   ├── make_report.py           # MAKE_REPORT
│   └── plot_style.py            # shared figure style
├── report/report_template.md    # report wording
├── data/
│   ├── sources.csv              # download URLs and md5
│   ├── crc_driver_genes.tsv     # driver panel
│   └── raw/                     # downloads (created by the run)
├── expected_results/            # reference outputs
├── cms/run_cms.R                # supplementary: CMS composition
├── survival/run_km.R            # supplementary: survival
├── results/                     # outputs (created by the run)
│   ├── 01_mutation_matrix/
│   ├── 02_clustering/
│   ├── 03_counts/
│   ├── 04_differential_expression/
│   ├── 05_report/
│   └── pipeline_info/           # runtime, timeline and DAG
└── work/                        # Nextflow cache and conda environments (created by the run)
```

## Limitations

- The clusters overlap (silhouette 0.22), and silhouette cannot test whether k = 1 is better.
- C2 membership and hypermutation almost coincide, so expression differences cannot be attributed to specific driver mutations.
- MSI status is not available in this pipeline; hypermutation is defined as at least 10 protein-altering mutations per Mb.
- Fold changes are not shrunken, and tumour purity is not modelled.
- All tumours come from TCGA; the clusters were not validated in an independent cohort.

## References

1. Guinney J, et al. The consensus molecular subtypes of colorectal cancer. *Nature Medicine* 21, 1350-1356 (2015). https://doi.org/10.1038/nm.3967
2. Nunes L, et al. Prognostic genome and transcriptome signatures in colorectal cancers. *Nature* 633, 137-146 (2024). https://doi.org/10.1038/s41586-024-07769-3
3. Hennig C. Cluster-wise assessment of cluster stability. *Computational Statistics & Data Analysis* 52, 258-271 (2007). https://doi.org/10.1016/j.csda.2006.11.025
