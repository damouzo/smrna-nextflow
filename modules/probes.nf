/*
 * Read probes (read_probe.R).
 * Counts reads derived from given probe sequences (shRNA constructs, etc.) in
 * the mirtrace collapsed FASTA, and optionally dumps the top 'unknown'
 * sequences per sample. Generic: no knockdown-specific logic.
 */

process READ_PROBES {

    tag "read probes"
    label 'process_medium'

    publishDir "${params.outdir}", mode: 'copy', pattern: 'deliverables/probes/**'

    when:
    params.nfcore_dir

    input:
    path nfcore_dir
    path samplesheet
    path probe_file
    path probe_script
    val  probe_mismatches
    val  probe_min_len
    val  top_unknown_n
    val  top_hits_n
    val  seed

    output:
    path 'deliverables/probes/**', emit: deliverables, arity: '0..*'

    script:
    """
    Rscript ${probe_script} \\
        --nfcore_dir ${nfcore_dir} \\
        --samplesheet ${samplesheet} \\
        --probe_file ${probe_file} \\
        --probe_mismatches ${probe_mismatches} \\
        --probe_min_len ${probe_min_len} \\
        --top_unknown_n ${top_unknown_n} \\
        --top_hits_n ${top_hits_n} \\
        --seed ${seed}
    """
}
