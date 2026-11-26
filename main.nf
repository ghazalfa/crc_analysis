#!/usr/bin/env nextflow

/*
 * crc-mutation-clusters
 *
 *   FETCH (setup, one task per file in data/sources.csv)
 *
 *   Module 1, mutation clusters:
 *     BUILD_MUTATION_MATRIX -> CLUSTER_MUTATIONS -> PLOT_CLUSTERS
 *
 *   Module 2, differential expression between the clusters:
 *     CLUSTER_MUTATIONS -> PREPARE_COUNTS -> DIFFERENTIAL_EXPRESSION -> PLOT_DE
 *
 *   MAKE_REPORT: all figures, algorithms and hyperparameters in results/05_report
 *
 *   Module 3 (next): pathway analysis of the DE results
 */

nextflow.enable.dsl = 2

process FETCH {
    tag "$name"
    storeDir "${params.cache_dir}"      // downloads are reused across runs

    input:
    tuple val(name), val(role), val(url), val(md5)

    output:
    tuple val(role), path(name)

    script:
    """
    python3 ${projectDir}/src/fetch_file.py --name "${name}" --url "${url}" --md5 "${md5}" --project-dir "${projectDir}"
    """
}

// fingerprint of src/, passed to every process as an input so that editing a script invalidates -resume
def scriptsHash() {
    def scripts = files("${projectDir}/src/*.py").sort { f -> f.name }
    return scripts.collect { f -> f.text }.join('').md5()
}

process BUILD_MUTATION_MATRIX {
    publishDir "${params.outdir}/01_mutation_matrix", mode: 'copy'

    input:
    path mutations
    path expression
    path clinical
    path genes
    val scripts

    output:
    path "mutation_matrix.tsv",   emit: matrix
    path "mutation_types.tsv",    emit: types
    path "sample_annotation.tsv", emit: annotation
    path "gene_summary.tsv"
    path "build_summary.json",    emit: summary

    script:
    """
    # scripts: ${scripts}
    python3 ${projectDir}/src/build_mutation_matrix.py \\
        --mutations ${mutations} \\
        --expression ${expression} \\
        --clinical ${clinical} \\
        --genes ${genes} \\
        --min-freq ${params.min_freq} \\
        --exome-mb ${params.exome_mb} \\
        --hypermut-tmb ${params.hypermut_tmb} \\
        --exclude-hypermutated ${params.exclude_hypermutated} \\
        --require-expression ${params.require_expression}
    """
}

process CLUSTER_MUTATIONS {
    publishDir "${params.outdir}/02_clustering", mode: 'copy'

    input:
    path matrix
    path annotation
    val scripts

    output:
    path "cluster_assignments.tsv",     emit: assignments
    path "linkage.tsv",                 emit: linkage
    path "linkage_samples.tsv",         emit: linkage_samples
    path "k_selection.tsv",             emit: k_selection
    path "cluster_gene_enrichment.tsv", emit: enrichment
    path "cluster_summary.tsv",         emit: summary
    path "cluster_summary.json",        emit: info
    path "cluster_assignments_all_k.tsv"
    path "cluster_stability.tsv"
    path "cluster_annotation_tests.tsv", emit: tests

    script:
    """
    # scripts: ${scripts}
    python3 ${projectDir}/src/cluster_mutations.py \\
        --matrix ${matrix} \\
        --annotation ${annotation} \\
        --linkage ${params.linkage} \\
        --k-min ${params.k_min} --k-max ${params.k_max} --k ${params.k} \\
        --min-cluster-size ${params.min_cluster_size} \\
        --n-boot ${params.n_boot} --seed ${params.seed}
    """
}

process PLOT_CLUSTERS {
    input:
    path types
    path annotation
    path assignments
    path linkage
    path linkage_samples
    path k_selection
    path enrichment
    path summary
    path info
    path genes
    val scripts

    output:
    path "*.png", emit: figures
    path "oncoprint.pdf", emit: pdf

    script:
    """
    # scripts: ${scripts}
    python3 ${projectDir}/src/plot_clusters.py \\
        --types ${types} --annotation ${annotation} --assignments ${assignments} \\
        --linkage ${linkage} --linkage-samples ${linkage_samples} \\
        --k-selection ${k_selection} --enrichment ${enrichment} --summary ${summary} \\
        --cluster-info ${info} --genes ${genes} --hypermut-tmb ${params.hypermut_tmb}
    """
}

process PREPARE_COUNTS {
    publishDir "${params.outdir}/03_counts", mode: 'copy', pattern: "{sample_metadata.tsv,gene_annotation.tsv,prepare_counts_summary.json}"

    input:
    path expression
    path probemap
    path assignments
    path annotation
    val scripts

    output:
    path "counts.tsv.gz",        emit: counts
    path "sample_metadata.tsv",  emit: metadata
    path "gene_annotation.tsv",  emit: genes
    path "prepare_counts_summary.json", emit: summary

    script:
    """
    # scripts: ${scripts}
    python3 ${projectDir}/src/prepare_counts.py \\
        --expression ${expression} --probemap ${probemap} \\
        --assignments ${assignments} --annotation ${annotation} \\
        --min-count ${params.min_count} --min-site ${params.min_site}
    """
}

process DIFFERENTIAL_EXPRESSION {
    label 'parallel'
    publishDir "${params.outdir}/04_differential_expression", mode: 'copy'

    input:
    path counts
    path metadata
    path genes
    val scripts

    output:
    path "de_*_vs_*.tsv",               emit: results
    path "de_summary.tsv",              emit: summary
    path "top_genes_log2norm.tsv",      emit: top
    path "expression_pca.tsv",          emit: pca
    path "expression_pca_variance.tsv", emit: pca_variance

    script:
    """
    # scripts: ${scripts}
    python3 ${projectDir}/src/differential_expression.py \\
        --counts ${counts} --metadata ${metadata} --genes ${genes} \\
        --alpha ${params.de_alpha} --min-lfc ${params.de_min_lfc} \\
        --covariates ${params.de_covariates} --cpus ${task.cpus}
    """
}

process PLOT_DE {
    input:
    path de_results
    path summary
    path top
    path metadata
    path pca
    path pca_variance
    val scripts

    output:
    path "*.png", emit: figures

    script:
    """
    # scripts: ${scripts}
    python3 ${projectDir}/src/plot_de.py \\
        --de ${de_results} --summary ${summary} --top-genes ${top} --metadata ${metadata} \\
        --pca ${pca} --pca-variance ${pca_variance} \\
        --alpha ${params.de_alpha} --min-lfc ${params.de_min_lfc}
    """
}

// report.md: the editable wording in templates/report_template.md, filled with this run's numbers
process MAKE_REPORT {
    publishDir "${params.outdir}/05_report", mode: 'copy'

    input:
    path template
    path figures
    path build_summary
    path info
    path cluster_summary
    path k_selection
    path tests
    path counts_summary
    path de_summary
    val settings
    val scripts

    output:
    path "report.md"
    path "figures"

    script:
    """
    # scripts: ${scripts}
    python3 ${projectDir}/src/make_report.py --template ${template} \\
        --figures ${figures} --build-summary ${build_summary} --cluster-info ${info} \\
        --cluster-summary ${cluster_summary} --k-selection ${k_selection} --annotation-tests ${tests} \\
        --counts-summary ${counts_summary} --de-summary ${de_summary} \\
        --settings ${settings}
    """
}

// all fetched files of one role (column 2 of sources.csv), as a single sorted list
def byRole(files, role) {
    return files.filter { r, f -> r == role }.map { r, f -> f }.collect(sort: true)
}

workflow {
    scripts = scriptsHash()
    genes = file(params.driver_genes, checkIfExists: true)

    files = channel.fromPath(params.sources, checkIfExists: true)
        .splitCsv(header: true)
        .map { row -> tuple(row.name, row.role, row.url, row.md5) }
    files = FETCH(files)

    BUILD_MUTATION_MATRIX(
        byRole(files, 'mutation'), byRole(files, 'expression'), byRole(files, 'clinical'), genes, scripts
    )

    CLUSTER_MUTATIONS(BUILD_MUTATION_MATRIX.out.matrix, BUILD_MUTATION_MATRIX.out.annotation, scripts)

    PLOT_CLUSTERS(
        BUILD_MUTATION_MATRIX.out.types,
        BUILD_MUTATION_MATRIX.out.annotation,
        CLUSTER_MUTATIONS.out.assignments,
        CLUSTER_MUTATIONS.out.linkage,
        CLUSTER_MUTATIONS.out.linkage_samples,
        CLUSTER_MUTATIONS.out.k_selection,
        CLUSTER_MUTATIONS.out.enrichment,
        CLUSTER_MUTATIONS.out.summary,
        CLUSTER_MUTATIONS.out.info,
        genes, scripts
    )

    PREPARE_COUNTS(
        byRole(files, 'expression'), byRole(files, 'annotation'),
        CLUSTER_MUTATIONS.out.assignments, BUILD_MUTATION_MATRIX.out.annotation, scripts
    )

    DIFFERENTIAL_EXPRESSION(PREPARE_COUNTS.out.counts, PREPARE_COUNTS.out.metadata, PREPARE_COUNTS.out.genes, scripts)

    PLOT_DE(
        DIFFERENTIAL_EXPRESSION.out.results,
        DIFFERENTIAL_EXPRESSION.out.summary,
        DIFFERENTIAL_EXPRESSION.out.top,
        PREPARE_COUNTS.out.metadata,
        DIFFERENTIAL_EXPRESSION.out.pca,
        DIFFERENTIAL_EXPRESSION.out.pca_variance, scripts
    )

    // parameters shown in the report
    settings = [
        'min_freq', 'exome_mb', 'hypermut_tmb', 'exclude_hypermutated', 'require_expression',
        'linkage', 'k_min', 'k_max', 'k', 'min_cluster_size', 'n_boot', 'seed',
        'min_count', 'min_site', 'de_covariates', 'de_alpha', 'de_min_lfc', 'max_cpus'
    ].collect { name -> "${name}=${params[name]}".toString() }.join(' ')

    MAKE_REPORT(
        file(params.report_template, checkIfExists: true),
        PLOT_CLUSTERS.out.figures.mix(PLOT_CLUSTERS.out.pdf, PLOT_DE.out.figures).collect(),
        BUILD_MUTATION_MATRIX.out.summary,
        CLUSTER_MUTATIONS.out.info,
        CLUSTER_MUTATIONS.out.summary,
        CLUSTER_MUTATIONS.out.k_selection,
        CLUSTER_MUTATIONS.out.tests,
        PREPARE_COUNTS.out.summary,
        DIFFERENTIAL_EXPRESSION.out.summary,
        settings, scripts
    )
}
