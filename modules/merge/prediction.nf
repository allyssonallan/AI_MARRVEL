#!/usr/bin/env nextflow
nextflow.enable.dsl=2
// This process merges scores from different chromosomes into a single matrix file
// It combines the tier information and compressed scores into a unified output
process PREDICTION {
    publishDir "${params.outdir}/${params.run_id}/prediction/", mode: "copy"

    input:
    path merged_matrix  
    path merged_compressed_scores  

    path ref_predict_new_dir
    path ref_model_inputs_dir

    output:
    path "conf_4Model/*.csv"
    path "conf_4Model/integrated/*.csv"

    script:
    """
    mkdir final_matrix_expanded
    mkdir conf_4Model

    run_final.py ${params.run_id}
    merge_rm.py ${params.run_id}
    extraModel_main.py -id ${params.run_id}
    """
}