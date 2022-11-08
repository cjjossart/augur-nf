process AUGUR_EXPORT {
    
    container = 'staphb/augur:16.0.3'

    input:
    file(refined_tree)
    file(metadata)
    file(branch_lengths)
    file(ancestral_nt)
    file(ancestral_aa)
    file(colors)
    file(lat_long)
    file(config)
    file(traits)

    output:
    file "auspice.json" , emit: auspice_json

    script:
    if(params.traits){
    """
    augur export v2 \
    --tree ${refined_tree} \
    --metadata ${metadata} \
    --node-data ${branch_lengths} \
                ${ancestral_nt} \
                ${ancestral_aa} \
                ${traits} \
    --colors ${colors} \
    --lat-longs ${lat_long} \
    --auspice-config ${config} \
    --output auspice.json \
    ${params.export_addn}
    """
    } else{
    """
    augur export v2 \
        --tree ${refined_tree} \
        --metadata ${metadata} \
        --node-data ${branch_lengths} \
                    ${ancestral_nt} \
                    ${ancestral_aa} \
        --colors ${colors} \
        --lat-longs ${lat_long} \
        --auspice-config ${config} \
        --output auspice.json \
        --color-by-metadata ${params.metadata_color}
        ${params.export_addn}
    """
    }
}
