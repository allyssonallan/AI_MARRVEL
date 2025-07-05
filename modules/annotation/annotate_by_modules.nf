#!/usr/bin/env nextflow
nextflow.enable.dsl=2
// This process annotates a VEP file using the Phrank modules for disease and gene
// annotation. It uses HGMD and OMIM similarity files to compute patient HPO similarity scores
// and outputs a CSV file with the scores.
process ANNOTATE_BY_MODULES {
    tag "${vep.simpleName}"

    input:
    path vep
    path hgmd_sim, stageAs: "hgmd_sim.tsv"
    path omim_sim, stageAs: "omim_sim.tsv"
    path ref_annot_dir

    output:
    path "${vep.baseName}_scores.csv", emit: scores

    script:
    """
    feature.py \\
        -patientHPOsimiOMIM $omim_sim \\
        -patientHPOsimiHGMD $hgmd_sim \\
        -varFile ${vep} \\
        -inFileType vepAnnotTab \\
        -patientFileType one \\
        -genomeRef ${params.ref_ver} \\
        -diseaseInh AD \\
        -modules curate,conserve
    
        mv scores.csv ${vep.baseName}_scores.csv
    """
}