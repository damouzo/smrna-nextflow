/*
 * Phase 2 — full biotype composition (biotype_full.R).
 * Aligns the mirtrace collapsed FASTA with bowtie, intersects with bedtools and
 * assigns each read to the highest-priority overlapped feature.
 */

process BIOTYPE_FULL {

    tag "biotype composition"
    label 'process_high'

    publishDir "${params.outdir}", mode: 'copy', pattern: 'deliverables/biotype/**'
    publishDir "${params.outdir}", mode: 'copy', pattern: 'raw/biotype/**'

    when:
    params.run_biotype

    input:
    path nfcore_dir
    path samplesheet
    path contrasts
    path gtf
    path ncrna_bed
    path genome_fasta
    path biotype_script
    val  bowtie_index
    val  trim3p
    val  bowtie_mismatches
    val  bowtie_max_hits
    val  min_overlap_frac
    val  biotype_priority
    val  de_classes
    val  min_cpm
    val  min_total_count
    val  alpha
    val  lfc_flag
    val  shrink_type
    val  batch_col
    val  seed

    output:
    path 'deliverables/biotype/**', emit: deliverables, arity: '0..*'
    path 'raw/biotype/**',          emit: raw, arity: '0..*'

    script:
    def bi = (bowtie_index && bowtie_index.toString().size() > 0) ? "--bowtie_index ${bowtie_index}" : ""
    def bc = (batch_col && batch_col.toString().size() > 0) ? "--batch_col ${batch_col}" : ""
    """
    Rscript ${biotype_script} \\
        --nfcore_dir ${nfcore_dir} \\
        --genome_fasta ${genome_fasta} \\
        ${bi} \\
        --gtf ${gtf} \\
        --ncrna_bed ${ncrna_bed} \\
        --samplesheet ${samplesheet} \\
        --contrasts ${contrasts} \\
        --trim3p ${trim3p} \\
        --bowtie_mismatches ${bowtie_mismatches} \\
        --bowtie_max_hits ${bowtie_max_hits} \\
        --min_overlap_frac ${min_overlap_frac} \\
        --biotype_priority "${biotype_priority}" \\
        --de_classes "${de_classes}" \\
        --min_cpm ${min_cpm} \\
        --min_total_count ${min_total_count} \\
        --alpha ${alpha} \\
        --lfc_flag ${lfc_flag} \\
        --shrink_type ${shrink_type} \\
        ${bc} \\
        --threads ${task.cpus} \\
        --seed ${seed}
    """
}
