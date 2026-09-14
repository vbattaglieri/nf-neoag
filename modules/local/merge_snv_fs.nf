process MERGE_SNV_FS {
    tag "$meta.id"
    label 'process_low'

    publishDir "${params.outdir}/merged", mode: 'copy'

    input:
    tuple val(meta), path(indel_consensus), path(snv_consensus)

    output:
    tuple val(meta), path("${meta.id}.neoantigens.consensus.all")         , emit: all
    tuple val(meta), path("${meta.id}.neoantigens.consensus.all.expanded"), emit: expanded

    script:
    """
    neo-merge \\
        ${indel_consensus} \\
        ${snv_consensus} \\
        ${meta.id}
    """

    stub:
    """
    touch ${meta.id}.neoantigens.consensus.all
    touch ${meta.id}.neoantigens.consensus.all.expanded
    """
}
