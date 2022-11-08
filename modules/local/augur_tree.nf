process tree{

    container = 'staphb/augur:16.0.3'

    input:
    file(alignment)

    output:
    file "raw_tree.newick" , emit: raw_tree

    shell:
    """
    augur tree \
        --alignment ${alignment} \
        --output raw_tree.newick \
        ${params.tree_addn}
    """
}