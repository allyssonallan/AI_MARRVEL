#!/usr/bin/env nextflow
nextflow.enable.dsl=2
process JOIN_TIER_PHRANK {
    tag "${scores.simpleName}"

    input:
        path scores
        path phrank
        path ref_annot_dir
        path ref_var_tier_dir
        path ref_merge_expand_dir

    output:
        path "${scores.simpleName}_scores.txt.gz", emit: compressed_scores
        path "${scores.simpleName}_Tier.v2.tsv", emit: tier

    script:
        """
        mv $scores scores.csv
        VarTierDiseaseDBFalse.R ${params.ref_ver}
        generate_new_matrix_2.py ${params.run_id} ${params.ref_ver}
        mv scores.txt.gz  ${scores.simpleName}_scores.txt.gz
        mv Tier.v2.tsv ${scores.simpleName}_Tier.v2.tsv
        """
}