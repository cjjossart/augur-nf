process AUGUR_TRANSLATE {

    container = 'staphb/augur:16.0.3'

    input:
    file(refined_tree)
    file(ancestral_nt)
    file(reference)

    output:
    file "aa_muts.json" , emit: ancestral_aa

    shell:
    """
    augur translate \
        --tree ${refine_tree} \
        --ancestral-sequences ${ancestral_nt} \
        --reference-sequence ${reference} \
        --output aa_muts.json \
        ${params.translate_addn}
    """
}