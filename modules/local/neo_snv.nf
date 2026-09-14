process NEO_SNV {
    tag "$meta.id"
    label 'process_medium'

    publishDir "${params.outdir}/neo_snv", mode: 'copy'

    input:
    tuple val(meta), path(snv_table), path(cdna_tab)   // cdna_tab = refFlat_mRNA.<genome>.noalt.parsed.tab
    path netmhc_dir                                     // directory containing the netMHC executable

    output:
    tuple val(meta), path("${meta.id}.neoantigens.results.consensus"), emit: consensus
    path "versions.yml"                                              , emit: versions

    script:
    def hla = meta.hla
    def kmer_args = "--min-kmer ${params.min_kmer} --max-kmer ${params.max_kmer}"
    """
    export NETMHC="\$(readlink -e ${netmhc_dir})/${params.netmhc_bin}"

    neo-snv \\
        ${snv_table} \\
        ${cdna_tab} \\
        ${meta.id} \\
        ${hla} \\
        ${kmer_args}

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
