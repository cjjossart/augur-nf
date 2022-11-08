#!/usr/bin/env nextflow

//Description: Nextflow implementation of Nextstrain's Augur pipeline
//Available: https://github.com/nextstrain/augur
//Authors of this Nextflow: Abigail Shockey
//Email: abigail.shockey@slh.wisc.edu

include { AUGUR_FILTER      } from '../modules/local/augur_filter'
include { AUGUR_ALIGN       } from '../modules/local/augur_align'
include { AUGUR_TREE        } from '../modules/local/augur_tree'
include { AUGUR_REFINE      } from '../modules/local/augur_refine'
include { AUGUR_ANCESTRAL   } from '../modules/local/augur_ancestral'
include { AUGUR_TRANSLATE   } from '../modules/local/augur_translate'
include { AUGUR_TRAITS      } from '../modules/local/augur_traits'
include { AUGUR_EXPORT      } from '../modules/local/augur_export'

workflow AUGUR {
  take:
  multifasta
  reference
  metadata
  colors
  lat_long
  
  main:

  ch_versions = Channel.empty()
  ch_multifasta = file(params.multifasta)
  ch_reference = file(params.reference)
  ch_metadata = file(metadata)
  ch_colors = file(colors)
  ch_lat_long = file(lat_long)
  ch_exclude = file(params.exclude)
  ch_bed = file(params.bed)

  AUGUR_FILTER( ch_multifasta, ch_metadata, ch_exclude )

  AUGUR_ALIGN( AUGUR_FILTER.out.filtered_fasta, ch_reference, ch_bed )

  AUGUR_TREE( AUGUR_ALIGN.out.alignment )

  AUGUR_REFINE( AUGUR_TREE.out.raw_tree, AUGUR_ALIGN.out.alignment, ch_metadata )

  AUGUR_ANCESTRAL( AUGUR_REFINE.out.refined_tree, AUGUR_ALIGN.out.alignment ) 

  AUGUR_TRANSLATE( AUGUR_REFINE.out.refined_tree, AUGUR_ANCESTRAL.out.ancestral_aa, ch_reference )

  AUGUR_TRAITS( AUGUR_REFINE.out.traits, ch_metadata )

  AUGUR_EXPORT( AUGUR_REFINE.out.refined_tree, ch_metadata, AUGUR_REFINE.out.branch_lengths, AUGUR_ANCESTRAL.out.ancestral_nt, AUGUR_TRANSLATE.out.ancestral_aa, ch_colors, ch_lat_long, ch_config )


  emit:

  versions = ch_versions
  auspice_json = AUGUR_EXPORT.out.



}

// Input channels

Channel
    .fromPath("${params.reference}")
    .ifEmpty { exit 1, "Cannot find reference sequence in: ${params.reference}" }
    .into { reference_alignment; reference_translate }

Channel
    .fromPath( "${params.metadata}")
    .ifEmpty { exit 1, "Cannot find metadata file in: ${params.metadata}" }
    .into { filter_metadata; refine_tree_metadata; traits_metadata; metadata_export }

Channel
    .fromPath( "${params.colors}")
    .ifEmpty { exit 1, "Cannot find colors file in: ${params.colors}" }
    .set { colors_export }

Channel
    .fromPath( "${params.lat_long}")
    .ifEmpty { exit 1, "Cannot find latitude and longitude file in: ${params.lat_long}" }
    .set { lat_long_export }

Channel
    .fromPath( "${params.config}")
    .ifEmpty { exit 1, "Cannot find Auspice config file in: ${params.auspice_config}" }
    .set { config_export }

if (params.mask != "") {
    Channel
        .fromPath( "${params.mask}")
        .ifEmpty { exit 1, "Cannot find input BED file in: ${params.mask}" }
        .set { bed_file }
}

else {
    Channel
        .from( "${params.mask}")
        .set { bed_file }
}

// Filter sequences if sequences to drop are included
if (params.filter) {
    Channel
        .fromPath( "${params.sequences}")
        .ifEmpty { exit 1, "Cannot find input sequences in: ${params.sequences}" }
        .set { sequences }
    Channel
        .fromPath( "${params.filter}")
        .ifEmpty { exit 1, "Cannot find filter file in: ${params.filter}" }
        .set { dropped_strains }

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
}

else {
    Channel
        .fromPath( "${params.sequences}")
        .ifEmpty { exit 1, "Cannot find input sequences in: ${params.sequences}" }
        .set { input_sequences }
}

process align{
  publishDir "${params.outdir}/align", mode:'copy'

  input:
  file(fasta) from input_sequences
  file(ref) from reference_alignment
  file bed from bed_file

  output:
  file "aligned.fasta" into raw_tree_alignment, refine_tree_alignment, ancestral_alignment

  script:
  if(params.mask!=""){
    """
    augur align \
      --sequences ${fasta} \
      --reference-sequence ${ref} \
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
      --sequences ${fasta} \
      --reference-sequence ${ref} \
      --output aligned.fasta \
      --fill-gaps \
      --nthreads ${params.threads} \
      ${params.align_addn}
    """
  }
}

process tree{
  publishDir "${params.outdir}/tree", mode:'copy'

  input:
  file(msa) from raw_tree_alignment

  output:
  file "raw_tree.newick" into raw_tree

  shell:
  """
  augur tree \
    --alignment ${msa} \
    --output raw_tree.newick \
    ${params.tree_addn}
  """
}

process refine{
  publishDir "${params.outdir}/refine", mode:'copy'

  input:
  file(tree) from raw_tree
  file(msa) from refine_tree_alignment
  file(metadata) from refine_tree_metadata

  output:
  file "refined_tree.newick" into refined_tree_traits, refined_tree_ancestral, refined_tree_translate, refined_tree_export
  file "branch_lengths.json" into branch_lengths_export

  shell:
  """
  augur refine \
    --tree ${tree} \
    --alignment ${msa} \
    --metadata ${metadata} \
    --output-tree refined_tree.newick \
    --output-node-data branch_lengths.json \
    --timetree \
    --coalescent opt \
    --date-confidence \
    --date-inference marginal \
    --clock-filter-iqd 4 \
    ${params.refine_addn}
  """
}

process ancestral{
  publishDir "${params.outdir}/ancestral_nt", mode:'copy'

  input:
  file(tree) from refined_tree_ancestral
  file(msa) from ancestral_alignment

  output:
  file "nt_muts.json" into ancestral_nt_translate, ancestral_nt_export

  shell:
  """
  augur ancestral \
    --tree ${tree} \
    --alignment ${msa} \
    --output-node-data nt_muts.json \
    ${params.ancestral_addn}
  """
}

process translate{
  publishDir "${params.outdir}/ancestral_aa", mode:'copy'

  input:
  file(tree) from refined_tree_translate
  file(nt_muts) from ancestral_nt_translate
  file(ref) from reference_translate

  output:
  file "aa_muts.json" into ancestral_aa_export

  shell:
  """
  augur translate \
    --tree ${tree} \
    --ancestral-sequences ${nt_muts} \
    --reference-sequence ${ref} \
    --output aa_muts.json \
    ${params.translate_addn}
  """
}

if (params.traits) {
    process traits {
        publishDir "${params.outdir}/traits", mode:'copy'

        input:
        file(tree) from refined_tree_traits
        file(metadata) from traits_metadata

        output:
        file "traits.json" into traits_export
        
        shell:
        """
        augur traits \
            --tree ${tree} \
            --metadata ${metadata} \
            --output traits.json \
            --columns ${params.columns} \
            --confidence \
            ${params.traits_addn}
        """
    }
}

else {
    Channel
        .from( "${params.traits}")
        .set { traits_export }
}

process export {
  publishDir "${params.outdir}", mode:'copy'

  input:
  file(tree) from refined_tree_export
  file(metadata) from metadata_export
  file(branch_lengths) from branch_lengths_export
  file(nt_muts) from ancestral_nt_export
  file(aa_muts) from ancestral_aa_export
  file(colors) from colors_export
  file(lat_long) from lat_long_export
  file(config) from config_export
  file(traits) from traits_export

  output:
  file "auspice.json"

  script:
  if(params.traits){
  """
  augur export v2 \
  --tree ${tree} \
  --metadata ${metadata} \
  --node-data ${branch_lengths} \
              ${nt_muts} \
              ${aa_muts} \
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
    --tree ${tree} \
    --metadata ${metadata} \
    --node-data ${branch_lengths} \
                ${nt_muts} \
                ${aa_muts} \
    --colors ${colors} \
    --lat-longs ${lat_long} \
    --auspice-config ${config} \
    --output auspice.json \
    ${params.export_addn}
  """
  }
}
