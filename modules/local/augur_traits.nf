process AUGUR_TRAITS {

    container = 'staphb/augur:16.0.3'

    input:
    file(refined_tree)
    file(metadata)

    output:
    file "traits.json" , emit: traits

    shell:
    """
    augur traits \
        --tree ${refined_tree} \
        --metadata ${metadata} \
        --output traits.json \
        --columns ${params.columns} \
        --confidence \
        ${params.traits_addn}
    """
}