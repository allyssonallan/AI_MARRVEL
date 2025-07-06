#!/usr/bin/env nextflow
nextflow.enable.dsl=2
process ENSEMBL_TO_GENESYM {
    input:
        path ensmbl
        path ref_to_sym
        path ref_sorted_sym

    output:
        path "${params.run_id}-gene.txt", emit: gene

    script:
        """
        cat ${ensmbl} | sort -k5,5 | join -1 5 -2 1 - ${ref_to_sym}  | sed 's/ /\t/g' | cut -f2- > genesym.txt
        cat genesym.txt | cut -f5 | sort -u | join -t\t' -1 1 -2 2 - ${ref_sorted_sym} | cut -f2 | sort -u > ${params.run_id}-gene.txt
        """
}
