/*
 * Validate inputs (validate_inputs.R).
 * Reads the nf-core results dir plus samplesheet/contrasts and writes the audit,
 * the validated tables and the run_info metadata (sample_index, design_check,
 * software_versions, run_manifest) into run_info/.
 */

process VALIDATE_INPUTS {

    tag "validate inputs"
    label 'process_single'

    publishDir "${params.outdir}/run_info", mode: 'copy'

    when:
    params.nfcore_dir && params.samplesheet && params.contrasts

    input:
    path nfcore_dir
    path samplesheet
    path contrasts
    path validate_script
    val  batch_col
    val  expected_adapter
    val  outdir
    val  container
    val  nextflow_version

    output:
    path 'audit.md',                 emit: audit
    path 'samplesheet_validated.csv', emit: samplesheet
    path 'contrasts_validated.csv',  emit: contrasts
    path 'sample_index.tsv',         emit: sample_index
    path 'design_check/*.tsv',       emit: design_check
    path 'software_versions.tsv',    emit: software_versions
    path 'run_manifest.tsv',         emit: run_manifest

    script:
    """
    Rscript ${validate_script} \\
        --nfcore_dir ${nfcore_dir} \\
        --samplesheet ${samplesheet} \\
        --contrasts ${contrasts} \\
        --batch_col "${batch_col}" \\
        --expected_adapter ${expected_adapter} \\
        --source_dir ${params.nfcore_dir} \\
        --outdir ${outdir} \\
        --container ${container} \\
        --nextflow_version ${nextflow_version}
    """
}

/*
 * Build counts (build_counts.R).
 * Produces the miRNA count matrix, the QC table and the sensitivity dumps.
 */

process BUILD_COUNTS {

    tag "build counts"
    label 'process_medium'

    publishDir "${params.outdir}", mode: 'copy'

    when:
    params.nfcore_dir

    input:
    path nfcore_dir
    path samplesheet_validated
    path build_script
    val  isomir_exclude_variants
    val  isomir_max_shift_5p
    val  isomir_max_shift_3p
    val  isomir_min_len

    output:
    path 'deliverables/counts/counts_mirtop.tsv', emit: counts_mirtop
    path 'deliverables/qc/qc_table.tsv',          emit: qc_table
    path 'raw/counts/counts_mature.tsv',          emit: counts_mature
    path 'raw/counts/isomir_filtered.tsv.gz',     emit: isomir_filtered

    script:
    """
    Rscript ${build_script} \\
        --nfcore_dir ${nfcore_dir} \\
        --samplesheet ${samplesheet_validated} \\
        --isomir_exclude_variants "${isomir_exclude_variants}" \\
        --isomir_max_shift_5p ${isomir_max_shift_5p} \\
        --isomir_max_shift_3p ${isomir_max_shift_3p} \\
        --isomir_min_len ${isomir_min_len}
    """
}
