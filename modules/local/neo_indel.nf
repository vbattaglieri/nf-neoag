process NEO_INDEL {
    tag "$meta.id"
    label 'process_medium'

    publishDir "${params.outdir}/neo_indel", mode: 'copy'

    input:
    tuple val(meta), path(indel_table), path(yar)   // yar = refFlat_mRNA.<genome>.yar
    path netmhc_dir                                  // directory containing the netMHC executable

    output:
    tuple val(meta), path("${meta.id}.neoantigens.results.consensus"), emit: consensus
    path "versions.yml"                                              , emit: versions

    script:
    def hla = meta.hla
    def kmer_args = "--min-kmer ${params.min_kmer} --max-kmer ${params.max_kmer}"
    """
    export NETMHC="\$(readlink -e ${netmhc_dir})/${params.netmhc_bin}"

    neo-indel \\
        ${indel_table} \\
        ${yar} \\
        ${meta.id} \\
        ${hla} \\
        ${kmer_args} \\
        --dist ${params.indel_family_dist}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        neoag: \$(python -c "import neoag; print(neoag.__version__)")
        python: \$(python --version 2>&1 | sed 's/Python //')
        netmhc: "${params.netmhc_bin}"
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.neoantigens.results.consensus
    touch versions.yml
    """
}
