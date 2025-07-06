#!/usr/bin/env nextflow
nextflow.enable.dsl=2
process MERGE_SCORES_BY_CHROMOSOME {
    publishDir "${params.outdir}/${params.run_id}/merged", mode: "copy"

    input:
        path phrank
        path tier
        path compressed_scores
        path ref_annot_dir
        path ref_mod5_diffusion_dir
        path ref_merge_expand_dir

    output:
        path "${params.run_id}.matrix.txt", emit: merged_matrix
        path "scores.txt.gz", emit: merged_compressed_scores

    script:
        """
        # Merge matrices
        cat ${tier} > Tier.v2.tsv

        # Merge compressed scores
        zcat -f ${compressed_scores} | gzip > scores.txt.gz

        post_processing.py ${params.run_id} ${params.ref_ver}
        """
}