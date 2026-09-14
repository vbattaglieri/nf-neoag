#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * nf-neoag — neoantigen prediction from SNVs (step 10) and Indels (step 11)
 * Nextflow port of the original TMB/neoantigens SOP.
 * Original scripts & method: see CITATIONS.md (credit to the original author).
 */

include { NEO_SNV      } from './modules/local/neo_snv'
include { NEO_INDEL    } from './modules/local/neo_indel'
include { MERGE_SNV_FS } from './modules/local/merge_snv_fs'

// ------------------------------------------------------------------
// Parameter checks
// ------------------------------------------------------------------
if (!params.input)  { error "Missing --input samplesheet (CSV)." }
if (!params.netmhc) { error "Missing --netmhc: path to the netMHC-4.0 install directory (license-restricted, bring-your-own)." }

def netmhc_ch = Channel.value( file(params.netmhc, checkIfExists: true) )

// Resolve a reference file for a given genome + key from conf/genomes.config
def ref_for = { genome, key ->
    def g = params.genomes?.get(genome)
    if (!g)      { error "Unknown genome '${genome}'. Configured: ${params.genomes?.keySet()}" }
    if (!g[key]) { error "Genome '${genome}' has no '${key}' set. Provide it in conf/genomes.config or via --genomes.${genome}.${key} (see docs/references.md)." }
    return file(g[key], checkIfExists: true)
}

workflow {

    // ------------------------------------------------------------------
    // Parse samplesheet: sample,genome,snv_table,indel_table,hla
    // ------------------------------------------------------------------
    rows = Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true, strip: true)
        .map { row ->
            if (!row.sample) { error "A row is missing the 'sample' column: ${row}" }
            if (!row.genome) { error "Sample '${row.sample}' is missing the 'genome' column." }
            if (!row.hla)    { error "Sample '${row.sample}' is missing the 'hla' column." }
            def meta = [ id: row.sample, genome: row.genome, hla: row.hla ]
            [ meta, row.snv_table?.trim(), row.indel_table?.trim() ]
        }

    // ------------------------------------------------------------------
    // Step 10 — neo from SNVs
    // ------------------------------------------------------------------
    snv_in = rows
        .filter { meta, snv, indel -> snv }
        .map    { meta, snv, indel ->
            [ meta, file(snv, checkIfExists: true), ref_for(meta.genome, 'cdna_tab') ]
        }

    NEO_SNV( snv_in, netmhc_ch )

    // ------------------------------------------------------------------
    // Step 11 — neo from Indels
    // ------------------------------------------------------------------
    indel_in = rows
        .filter { meta, snv, indel -> indel }
        .map    { meta, snv, indel ->
            [ meta, file(indel, checkIfExists: true), ref_for(meta.genome, 'yar') ]
        }

    NEO_INDEL( indel_in, netmhc_ch )

    // ------------------------------------------------------------------
    // Optional merge (indel + SNV consensus) — requires both for a sample
    // ------------------------------------------------------------------
    if (params.merge) {
        merge_in = NEO_INDEL.out.consensus
            .join( NEO_SNV.out.consensus )   // by meta -> [meta, indel_cons, snv_cons]
        MERGE_SNV_FS( merge_in )
    }

    // ------------------------------------------------------------------
    // Collate versions
    // ------------------------------------------------------------------
    NEO_SNV.out.versions
        .mix( NEO_INDEL.out.versions )
        .collectFile(name: 'software_versions.yml', storeDir: params.outdir)
}

workflow.onComplete {
    log.info ( workflow.success
        ? "\n[nf-neoag] Completed. Results in: ${params.outdir}\n"
        : "\n[nf-neoag] Failed. See .nextflow.log\n" )
}
