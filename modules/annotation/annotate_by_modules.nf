#!/usr/bin/env nextflow
nextflow.enable.dsl=2
process ANNOTATE_BY_MODULES {
    tag "${vep.simpleName}"

    input:
        path vep
        path hgmd_sim
        path omim_sim
        path ref_annot_dir

    output:
        path "${vep.baseName}_scores.csv", emit: scores

    script:
        """
        feature.py \
            -patientHPOsimiOMIM ${omim_sim} \
            -patientHPOsimiHGMD ${hgmd_sim} \
            -varFile ${vep} \
            -inFileType vepAnnotTab \
            -patientFileType one \
            -genomeRef ${params.ref_ver} \
            -diseaseInh AD \
            -modules curate,conserve
        
        mv scores.csv ${vep.baseName}_scores.csv
        """
}