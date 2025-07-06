#!/usr/bin/env nextflow
nextflow.enable.dsl=2
process HPO_SIM {
    input:
        path hpo
        path omim_hgmd_phen
        path omim_obo
        path omim_genemap2
        path omim_pheno

    output:
        path "${params.run_id}-cz", emit: hgmd_sim
        path "${params.run_id}-dx", emit: omim_sim

    script:
        """
        if [[ -z \$(egrep 'HP:[0-9]{7}' $hpo) ]] ; then
            echo "HP:0000001" > $hpo
        fi

        phenoSim.R \
            ${hpo} \
            ${omim_hgmd_phen} \
            ${omim_obo} \
            ${omim_genemap2} \
            ${omim_pheno} \
            ${params.run_id}-cz \
            ${params.run_id}-dx
        """
}