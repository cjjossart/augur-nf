process AUGUR_ALIGN {

    container = 'staphb/augur:16.0.3'

    input:
    file(filtered_fasta)
    file(reference)
    file(bed)

    output:
    file "aligned.fasta" , emit: alignment

    script:
    if(params.mask!=""){
        """
        augur align \
        --sequences ${filtered_fasta} \
        --reference-sequence ${reference} \
        --output pre_mask.fasta \
        --fill-gaps \
        --nthreads ${params.threads} \
        ${params.align_addn}
        
        augur mask \
        --sequences pre_mask.fasta \
        --mask ${bed} \
        --output aligned.fasta \
        ${params.mask_addn}
        """
    } else{
        """
        augur align \
        --sequences ${filtered_fasta} \
        --reference-sequence ${reference} \
        --output aligned.fasta \
        --fill-gaps \
        --nthreads ${params.threads} \
        ${params.align_addn}
        """
    }
}