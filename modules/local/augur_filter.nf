process AUGUR_FILTER {

    container = 'staphb/augur:16.0.3'

    input:
    file(multifasta)
    file(metadata)
    file(exclude)

    output:
    file "filtered.fasta", emit: filtered_fasta

    shell:
    """
    augur filter \
    --sequences ${multifasta} \
    --metadata ${metadata} \
    --exclude ${exclude} \
    --output filtered.fasta \
    --group-by ${params.group_by} \
    --sequences-per-group ${params.per_group} \
    ${params.filter_addn}
    """
}
