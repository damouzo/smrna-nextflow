/*
 * Phase 1 analysis modules: QC and composition, PCA and DE, summarize DE.
 */

/*
 * QC flags and mirTrace composition (qc_composition.R).
 * Samples are flagged, never excluded.
 */

process QC_COMPOSITION {

    tag "QC and composition"
    label 'process_single'

    publishDir "${params.outdir}", mode: 'copy'

    when:
    params.nfcore_dir

    input:
    path qc_table
    path samplesheet
    path contrasts
    path qc_script
    val  qc_min_mirna_frac
    val  qc_max_trna_frac
    val  qc_max_rrna_frac
    val  qc_min_mirna_reads
    val  seed

    output:
    path 'deliverables/qc/qc_flags.tsv',              emit: qc_flags
    path 'deliverables/qc/flagged_samples.txt',       emit: flagged_samples
    path 'deliverables/qc/plots/qc_overview.pdf',     emit: qc_overview_pdf
    path 'deliverables/qc/plots/qc_overview.png',     emit: qc_overview_png
    path 'deliverables/qc/composition_by_sample.tsv', emit: composition_by_sample
    path 'deliverables/qc/composition_delta.tsv',     emit: composition_delta
    path 'deliverables/qc/plots/composition_stacked_bar.pdf', emit: composition_stacked_bar_pdf
    path 'deliverables/qc/plots/composition_stacked_bar.png', emit: composition_stacked_bar_png
    path 'run_info/qa_report.tsv',                    emit: qa_report

    script:
    """
    Rscript ${qc_script} \\
        --qc_table ${qc_table} \\
        --samplesheet ${samplesheet} \\
        --contrasts ${contrasts} \\
        --qc_min_mirna_frac ${qc_min_mirna_frac} \\
        --qc_max_trna_frac ${qc_max_trna_frac} \\
        --qc_max_rrna_frac ${qc_max_rrna_frac} \\
        --qc_min_mirna_reads ${qc_min_mirna_reads} \\
        --seed ${seed}
    """
}

/*
 * DESeq2 per norm_group (deseq2_group.R).
 * One invocation per group (plus 'all', explore only) and every run.
 * main writes deliverables; sensitivities stay in raw/.
 */

process DESEQ2_GROUP {

    tag "DESeq2 $group"
    label 'process_medium'

    publishDir "${params.outdir}", mode: 'copy', pattern: 'deliverables/**'
    publishDir "${params.outdir}", mode: 'copy', pattern: 'raw/**'

    when:
    params.nfcore_dir

    input:
    tuple val(group), path(samplesheet), path(contrasts), path(counts_mirtop), path(counts_mature), path(flagged), path(deseq_script)
    val  count_source
    val  min_cpm
    val  min_total_count
    val  alpha
    val  lfc_flag
    val  shrink_type
    val  n_top_var
    val  batch_col
    val  seed

    output:
    path 'deliverables/**',   emit: deliverables, arity: '0..*'
    path 'raw/**',            emit: raw, arity: '0..*'
    path 'deseq_group.done',  emit: done

    script:
    def batch_arg = batch_col ? "--batch_col ${batch_col}" : ""
    """
    Rscript ${deseq_script} \\
        --counts_mirtop ${counts_mirtop} \\
        --counts_mature ${counts_mature} \\
        --samplesheet ${samplesheet} \\
        --contrasts ${contrasts} \\
        --flagged ${flagged} \\
        --group ${group} \\
        --count_source ${count_source} \\
        --min_cpm ${min_cpm} \\
        --min_total_count ${min_total_count} \\
        --alpha ${alpha} \\
        --lfc_flag ${lfc_flag} \\
        --shrink_type ${shrink_type} \\
        --n_top_var ${n_top_var} \\
        ${batch_arg} \\
        --seed ${seed}
    """
}

/*
 * Summarize DE runs and concordance (summarize_de.R).
 * Reads the published outdir; gated by the collected DESeq2 completion stamps.
 */

process SUMMARIZE_DE {

    tag "summarize DE"
    label 'process_single'

    publishDir "${params.outdir}", mode: 'copy', pattern: 'deliverables/de/**'

    input:
    path sum_script
    val done
    val outdir
    val alpha

    output:
    path 'deliverables/de/de_summary.tsv',  emit: de_summary
    path 'deliverables/de/concordance.tsv', emit: concordance

    script:
    """
    Rscript ${sum_script} \\
        --outdir ${outdir} \\
        --alpha ${alpha}
    """
}
