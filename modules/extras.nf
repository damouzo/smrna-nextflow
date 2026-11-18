/*
 * Phase 3 modules and the report: report, isomiRs, functional.
 */

/*
 * HTML report (report.Rmd).
 * Reads only the published outdir; fails if a phase-1 file is missing.
 */

process REPORT {

    tag "report"
    label 'process_single'

    publishDir "${params.outdir}", mode: 'copy', pattern: 'deliverables/report/**'

    input:
    path report_script
    val done
    val outdir
    val run_biotype
    val run_isomir
    val run_functional

    output:
    path 'deliverables/report/smrna_report.html', emit: report

    script:
    def rb = run_biotype ? 'TRUE' : 'FALSE'
    def ri = run_isomir ? 'TRUE' : 'FALSE'
    def rf = run_functional ? 'TRUE' : 'FALSE'
    """
    mkdir -p deliverables/report
    Rscript -e "rmarkdown::render('${report_script}', output_file='smrna_report.html', output_dir=file.path(getwd(),'deliverables/report'), intermediates_dir=tempdir(), params=list(outdir='${outdir}', run_biotype=${rb}, run_isomir=${ri}, run_functional=${rf}))"
    """
}

/*
 * isomiR profile (isomir_profile.R, phase 3, run_isomir).
 * Requires reads with a clean 3' end.
 */

process ISOMIR_PROFILE {

    tag "isomiR profile"
    label 'process_medium'

    publishDir "${params.outdir}", mode: 'copy', pattern: 'deliverables/isomir/**'

    when:
    params.run_isomir

    input:
    path isomir_filtered
    path samplesheet
    path contrasts
    path isomir_script
    val  isomir_min_reads
    val  min_cpm
    val  min_total_count
    val  alpha
    val  lfc_flag
    val  shrink_type
    val  seed

    output:
    path 'deliverables/isomir/**', emit: deliverables, arity: '0..*'

    script:
    """
    Rscript ${isomir_script} \\
        --isomir_filtered ${isomir_filtered} \\
        --samplesheet ${samplesheet} \\
        --contrasts ${contrasts} \\
        --isomir_min_reads ${isomir_min_reads} \\
        --min_cpm ${min_cpm} \\
        --min_total_count ${min_total_count} \\
        --alpha ${alpha} \\
        --lfc_flag ${lfc_flag} \\
        --shrink_type ${shrink_type} \\
        --seed ${seed}
    """
}
