process filter{
    publishDir "${params.outdir}/filter", mode:'copy'

    input:
    file(fasta) from sequences
    file(metadata) from filter_metadata
    file(exclude) from dropped_strains

    output:
    file "filtered.fasta" into input_sequences

    shell:
    """
    augur filter \
    --sequences ${fasta} \
    --metadata ${metadata} \
    --exclude ${exclude} \
    --output filtered.fasta \
    --group-by ${params.group_by} \
    --sequences-per-group ${params.per_group} \
    ${params.filter_addn}
    """
}