# Driver-mutation clusters in colorectal cancer

TCGA colon (COAD) and rectal (READ) tumours are grouped by which driver genes are mutated, and gene expression is then compared between the groups. Numbers in the tables are read from this run's outputs.

- [Main finding](#main-finding)
- [1. Mutation matrix](#1-mutation-matrix)
- [2. Hierarchical clustering](#2-hierarchical-clustering)
- [3. Count preparation](#3-count-preparation)
- [4. Differential expression](#4-differential-expression)
- [Limitations](#limitations)
- [References](#references)
- [Software](#software)

## Main finding

**Cluster C2 is the hypermutated group.** With k = 2 the tumours split into a large, mostly non-hypermutated cluster (C1) and a smaller hypermutated cluster (C2). Both clusters pass the 0.75 stability guideline (see the table in section 2). The expression analysis is a single comparison, C2 against C1.

The separation is weak in absolute terms (see the silhouette in section 2), and a hypermutated group in colorectal cancer is already well known. Recheck this section if k or the gene panel changes.

## 1. Mutation matrix

Masked somatic mutations are reduced to a tumour by gene table of 0 (not mutated) and 1 (at least one protein-altering mutation). Only genes on a colorectal driver panel (Nunes et al. 2024 [1]) are used, because in an unrestricted matrix long genes are mutated often by chance. The panel contains oncogenes and tumour suppressors

{{table_sample_flow}}

{{panel_genes_kept}} of {{panel_genes}} panel genes passed the frequency filter; {{hypermutated_tumours}} tumours are hypermutated.

### Hyperparameters

| Parameter | Value | Why this value |
|---|---|---|
| min_freq | {{min_freq}} | Rarely mutated genes add little information about groups. Convention; not tuned. |
| exome_mb | {{exome_mb}} | Approximate size of the sequenced coding exome, used to express TMB per Mb. |
| hypermut_tmb | {{hypermut_tmb}} | Commonly used TMB-high cutoff. It only annotates tumours; it does not affect clustering unless exclude_hypermutated is true. |
| exclude_hypermutated | {{exclude_hypermutated}} | false keeps hypermutated tumours, because they form the main cluster of interest. true clusters only the other tumours. |
| require_expression | {{require_expression}} | Keeps clustered tumours usable for differential expression. |

## 2. Hierarchical clustering

**Distance.** Jaccard distance is 1 minus (genes mutated in both tumours / genes mutated in either). It ignores genes that neither tumour has mutated, which suits a matrix that is mostly zeros. It carries no information for a tumour with no mutated genes, so those tumours form a separate group, C0, and are not clustered.

**Linkage.** Clustering starts with each tumour alone and repeatedly merges two groups. Ward linkage merges the pair that least increases within-group variance, giving compact groups; average linkage uses the mean distance between groups, and complete linkage the largest. Ward is derived for Euclidean distance, so with Jaccard it is a heuristic. During development, average linkage mostly split off single tumours.

**Silhouette.** For each tumour, `(b - a) / max(a, b)`, where a is the mean distance to its own cluster and b to the nearest other cluster. Values near 1 are well separated, near 0 are on a boundary, below 0 fit another cluster better. {{k_text}}

**Stability.** Tumours are resampled with replacement and reclustered {{n_boot}} times (repeated draws of a tumour are used once). For each cluster, the Jaccard overlap with its best match is averaged. Guideline (Hennig 2007): above 0.75 stable, 0.6 to 0.75 a pattern with uncertain membership, at or below 0.5 dissolved.

**Gene and clinical tests.** Fisher's exact test per gene (cluster against the other clusters, Benjamini-Hochberg adjusted) and chi-square tests for clinical variables. The gene p-values are descriptive only: the same mutations defined the clusters.

### Hyperparameters

| Parameter | Value | Why this value |
|---|---|---|
| linkage | {{linkage}} | Gave groups large enough to compare; average linkage gave singletons. |
| k_min | {{k_min}} | Smallest possible split. |
| k_max | {{k_max}} | Upper limit of the k values compared. |
| k | {{k}} | 0 selects by silhouette; any other value is a user choice. 2 has the highest silhouette of the k values compared and is stable in both clusters. Larger k kept the same hypermutated group at k = 3 and 7 but split the rest into clusters with low stability. |
| min_cluster_size | {{min_cluster_size}} | Minimum size used in automatic selection of k, so groups can support expression tests. Judgement call. |
| n_boot | {{n_boot}} | Resamples averaged for stability; more gives less noisy scores. |
| seed | {{seed}} | Arbitrary; makes the resampling repeatable. |

### Results

{{table_clusters}}

Percentages of side exclude tumours with unknown site.

{{table_tests}}

![Oncoprint](figures/oncoprint.png)

*Oncoprint. Columns are tumours in dendrogram order, rows are genes. Pale cells are unmutated or unknown.*

![Silhouette by k](figures/k_selection.png)

*Mean silhouette for each number of clusters; the circled value was used.*

![Gene frequency per cluster](figures/cluster_gene_frequency.png)

*Percent of tumours with each gene mutated, per cluster. Stars mark enrichment (descriptive).*

![Cluster characteristics](figures/cluster_characteristics.png)

*Cluster composition and bootstrap stability.*

## 3. Count preparation

Expression is stored as log2(count + 1) and is converted back to integer counts, as the DESeq2 model requires. Low-count genes are removed. Tissue source site, read from the TCGA barcode, is recorded as a batch variable.

{{genes_kept}} of {{genes_in_file}} genes were kept ({{count_filter}}); {{tumours}} tumours; {{batch_levels}} site levels.

### Hyperparameters

| Parameter | Value | Why this value |
|---|---|---|
| min_count | {{min_count}} | Conventional low-count filter. The required number of tumours is the smallest cluster size, so genes expressed in one cluster survive. |
| min_site | {{min_site}} | Very small sites cannot be estimated reliably, so they are pooled as other. Judgement call. |

## 4. Differential expression

PyDESeq2 (a Python implementation of DESeq2) normalises samples with median-of-ratios size factors, models counts with a negative binomial distribution, shrinks each gene's dispersion towards a trend, and tests cluster differences with a Wald test. One model, `{{design}}`, is fitted to all clustered tumours and the clusters are compared pairwise. P-values are Benjamini-Hochberg adjusted. Fold changes are not shrunken.

### Hyperparameters

| Parameter | Value | Why this value |
|---|---|---|
| de_covariates | {{de_covariates}} | Adjusts for technical differences between TCGA collection sites. |
| de_alpha | {{de_alpha}} | Conventional 5% false discovery rate. |
| de_min_lfc | {{de_min_lfc}} | Keeps changes of at least two-fold when 1. Applied after testing, so it is an effect-size filter. |
| max_cpus | {{max_cpus}} | Speed only. |

### Results

C2_vs_C1 compares the hypermutated cluster with the other clustered tumours; positive fold changes are higher in C2.

{{table_de}}

![Significant genes per comparison](figures/de_gene_counts.png)

*Significant genes per comparison.*

![Volcano plots](figures/volcano_plots.png)

*Volcano plots. Dashed lines are the cutoffs; labelled genes have the largest combined significance and fold change.*

![Heatmap of top genes](figures/top_de_heatmap.png)

*Top 25 genes by adjusted p-value per comparison, z-scored log2 normalised counts (not batch-corrected).*

![Expression PCA](figures/expression_pca.png)

*PCA of the 2,000 most variable genes, log2 normalised counts (not batch-corrected).*

## Limitations

- Silhouette and stability show overlapping groups; check them above before treating a cluster as a subtype.
- C1 is a broad group: no further split of it was stable, so differences within it (such as KRAS status) are not analysed.
- A high stability score is easier to reach for a cluster that holds most of the tumours, as C1 does.
- Silhouette cannot evaluate k = 1, so it does not show that clusters exist.
- C2 reflects mutation burden as well as specific drivers: tumours with many mutations are more likely to have panel genes mutated. MSI status is not available in this pipeline.
- Site adjustment assumes every cluster is represented across sites; sites with few tumours of a cluster contribute little.
- Bulk expression reflects tissue mixture as well as tumour cells. Tumour purity is not modelled.

## References

1. Nunes L, Li F, Wu M, et al. Prognostic genome and transcriptome signatures in colorectal cancers. *Nature* 633, 137-146 (2024). https://doi.org/10.1038/s41586-024-07769-3

## Software

{{table_software}}
