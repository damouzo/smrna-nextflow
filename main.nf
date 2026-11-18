#!/usr/bin/env nextflow

/*
 * smrna-nextflow — downstream small RNA-seq analysis on nf-core/smrnaseq 2.4.x results.
 *
 * Phase 1: validate inputs -> counts -> QC/composition -> PCA/DE -> summary -> report
 * Phase 2: full biotype composition
 * Phase 3: isomiRs, functional enrichment
 */

nextflow.enable.dsl=2

include { VALIDATE_INPUTS; BUILD_COUNTS } from './modules/inputs.nf'
include { QC_COMPOSITION } from './modules/analysis.nf'
include { DESEQ2_GROUP } from './modules/analysis.nf'
include { SUMMARIZE_DE } from './modules/analysis.nf'
include { REPORT } from './modules/extras.nf'
include { ISOMIR_PROFILE } from './modules/extras.nf'
include { BIOTYPE_FULL } from './modules/reads.nf'
include { READ_PROBES } from './modules/probes.nf'

// Fail early with an actionable message when a top-level input is missing.
def requireInput(name, value) {
    if (value == null || value == '' ) {
        error("Missing required input '${name}'.\n" +
              "  Where: params.${name}\n" +
              "  Fix:   set it in your params.yaml and pass it with -params-file params.yaml")
    }
    def f = file(value)
    if (!f.exists()) {
        error("Input '${name}' not found: ${value}\n" +
              "  Where: params.${name}\n" +
              "  Fix:   check the path (relative paths resolve against the launch directory)")
    }
    return f
}

workflow {

    // ---- Run metadata for traceability -----------------------------------
    // Effective params and thresholds are dumped before any task runs, so a run
    // can be reproduced from run_info/ alone.
    def run_info_dir      = new File("${params.outdir}/run_info")
    def pipeline_info_dir = new File(run_info_dir, "pipeline_info")
    pipeline_info_dir.mkdirs()
    new File(pipeline_info_dir, "params.json").text =
        groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(params))

    def pkeys = params.keySet().findAll { k -> params[k] != null }.toList().sort()
    new File(run_info_dir, "params_used.tsv").text =
        "parameter\tvalue\n" + pkeys.collect { k -> "${k}\t${params[k]}" }.join("\n") + "\n"

    def thr_keys = ['count_source', 'isomir_exclude_variants', 'isomir_max_shift_5p',
                    'isomir_max_shift_3p', 'isomir_min_len', 'qc_min_mirna_frac',
                    'qc_max_trna_frac', 'qc_max_rrna_frac', 'qc_min_mirna_reads',
                    'sensitivity_exclude_flags', 'min_cpm', 'min_total_count', 'alpha',
                    'lfc_flag', 'shrink_type', 'n_top_var', 'isomir_min_reads',
                    'min_sig_mirnas', 'seed', 'trim3p', 'bowtie_mismatches',
                    'bowtie_max_hits', 'min_overlap_frac', 'probe_mismatches',
                    'probe_min_len', 'probe_top_unknown_n', 'probe_top_hits_n']
    new File(run_info_dir, "thresholds_used.tsv").text =
        "parameter\tvalue\n" + thr_keys.findAll { k -> params[k] != null }
            .collect { k -> "${k}\t${params[k]}" }.join("\n") + "\n"

    // Ship the output-tree READMEs with the results.
    new File("${params.outdir}/deliverables").mkdirs()
    new File("${params.outdir}/raw").mkdirs()
    ['deliverables', 'raw', 'run_info'].each { tree ->
        def readme = new File("${projectDir}/assets/readme_${tree}.md")
        if (readme.exists()) new File("${params.outdir}/${tree}/README.md").text = readme.text
    }

    def nxf_version = nextflow.version.toString()

    // ---- Inputs -----------------------------------------------------------
    ch_nfcore_dir  = requireInput('nfcore_dir',  params.nfcore_dir)
    ch_samplesheet = requireInput('samplesheet', params.samplesheet)
    ch_contrasts   = requireInput('contrasts',   params.contrasts)

    // ---- Validate inputs --------------------------------------------------
    // Writes audit.md, samplesheet_validated.csv and contrasts_validated.csv to run_info/.
    VALIDATE_INPUTS(ch_nfcore_dir, ch_samplesheet, ch_contrasts,
                    file("${projectDir}/bin/validate_inputs.R"), params.batch_col ?: '',
                    params.expected_adapter,
                    params.outdir, params.container ?: 'ghcr.io/damouzo/smrna-nextflow:latest',
                    nxf_version)

    // ---- Build counts (bin/build_counts.R) --------------------------------
    BUILD_COUNTS(
        ch_nfcore_dir,
        VALIDATE_INPUTS.out.samplesheet,
        file("${projectDir}/bin/build_counts.R"),
        params.isomir_exclude_variants,
        params.isomir_max_shift_5p,
        params.isomir_max_shift_3p,
        params.isomir_min_len
    )

    // ---- QC and composition (bin/qc_composition.R) ------------------------
    QC_COMPOSITION(
        BUILD_COUNTS.out.qc_table,
        VALIDATE_INPUTS.out.samplesheet,
        VALIDATE_INPUTS.out.contrasts,
        file("${projectDir}/bin/qc_composition.R"),
        params.qc_min_mirna_frac,
        params.qc_max_trna_frac,
        params.qc_max_rrna_frac,
        params.qc_min_mirna_reads,
        params.seed
    )

    // ---- PCA and DE (bin/deseq2_group.R) ----------------------------------
    // One process per norm_group, plus 'all' (explore only, no DE).
    ch_groups = Channel.of(tuple('all'))
        .mix(
            VALIDATE_INPUTS.out.samplesheet
                .splitCsv(header: true)
                .map { row -> tuple(row.norm_group) }
                .unique()
        )

    ch_deseq_files = VALIDATE_INPUTS.out.samplesheet
        .combine(VALIDATE_INPUTS.out.contrasts)
        .combine(BUILD_COUNTS.out.counts_mirtop)
        .combine(BUILD_COUNTS.out.counts_mature)
        .combine(QC_COMPOSITION.out.flagged_samples)

    ch_deseq_input = ch_groups
        .combine(ch_deseq_files)
        .map { g, meta, con, cm, cmat, flagged -> tuple(g, meta, con, cm, cmat, flagged, file("${projectDir}/bin/deseq2_group.R")) }

    DESEQ2_GROUP(
        ch_deseq_input,
        params.count_source,
        params.min_cpm,
        params.min_total_count,
        params.alpha,
        params.lfc_flag,
        params.shrink_type,
        params.n_top_var,
        params.batch_col ?: '',
        params.seed
    )

    // ---- Summarize DE (bin/summarize_de.R) --------------------------------
    SUMMARIZE_DE(
        file("${projectDir}/bin/summarize_de.R"),
        DESEQ2_GROUP.out.done.collect(),
        params.outdir,
        params.alpha
    )

    // ---- Full biotype composition (phase 2, bin/biotype_full.R) -----------
    if (params.run_biotype) {
        BIOTYPE_FULL(
            ch_nfcore_dir,
            VALIDATE_INPUTS.out.samplesheet,
            VALIDATE_INPUTS.out.contrasts,
            file(params.gtf),
            file(params.ncrna_bed),
            file(params.genome_fasta),
            file("${projectDir}/bin/biotype_full.R"),
            params.bowtie_index ?: '',
            params.trim3p,
            params.bowtie_mismatches,
            params.bowtie_max_hits,
            params.min_overlap_frac,
            params.biotype_priority,
            params.de_classes,
            params.min_cpm,
            params.min_total_count,
            params.alpha,
            params.lfc_flag,
            params.shrink_type,
            params.batch_col ?: '',
            params.seed
        )
    }

    // ---- isomiR profile (phase 3, bin/isomir_profile.R) -------------------
    if (params.run_isomir) {
        ISOMIR_PROFILE(
            BUILD_COUNTS.out.isomir_filtered,
            VALIDATE_INPUTS.out.samplesheet,
            VALIDATE_INPUTS.out.contrasts,
            file("${projectDir}/bin/isomir_profile.R"),
            params.isomir_min_reads,
            params.min_cpm,
            params.min_total_count,
            params.alpha,
            params.lfc_flag,
            params.shrink_type,
            params.seed
        )
    }

    // ---- Functional enrichment (phase 3, bin/functional.R) ---------------
    // ...

    // ---- Read probes (bin/read_probe.R) ----------------------------------
    // Counts reads derived from given probe sequences (shRNA constructs, etc.)
    // in the mirtrace collapsed FASTA. Runs when a probe file is given or when
    // the top-N unknown table is requested.
    if (params.probe_file || params.probe_top_unknown_n > 0) {
        READ_PROBES(
            ch_nfcore_dir,
            VALIDATE_INPUTS.out.samplesheet,
            params.probe_file ? file(params.probe_file) : file("${projectDir}/assets/probes.tsv"),
            file("${projectDir}/bin/read_probe.R"),
            params.probe_mismatches,
            params.probe_min_len,
            params.probe_top_unknown_n,
            params.probe_top_hits_n,
            params.seed
        )
    }

    // ---- Report (bin/report.Rmd) — runs last ------------------------------
    // Gated by phase 1 and, when enabled, by phase 2/3 and the probes module.
    ch_report_gate = Channel.empty()
    ch_report_gate = ch_report_gate.mix(SUMMARIZE_DE.out.de_summary.collect())
    if (params.run_biotype) ch_report_gate = ch_report_gate.mix(BIOTYPE_FULL.out.deliverables.collect())
    if (params.run_isomir)  ch_report_gate = ch_report_gate.mix(ISOMIR_PROFILE.out.deliverables.collect())
    if (params.probe_file || params.probe_top_unknown_n > 0) ch_report_gate = ch_report_gate.mix(READ_PROBES.out.deliverables.collect())
    REPORT(
        file("${projectDir}/bin/report.Rmd"),
        ch_report_gate.collect(),
        params.outdir,
        params.run_biotype,
        params.run_isomir,
        params.run_functional
    )

}
