#!/usr/bin/env nextflow
nextflow.enable.dsl=2
process GENESYM_TO_PHRANK {
    publishDir "${params.outdir}/${params.run_id}/phrank/", mode: 'copy'

    input:
        path gene
        path hpo
        path dagfile
        path disease_annotation
        path gene_annotation
        path disease_gene

    output:
        path "${params.run_id}.phrank.txt", emit: phrank

    script:
        """
        run_phrank.py \
            ${gene} \
            ${hpo} \
            ${dagfile} \
            ${disease_annotation} \
            ${gene_annotation} \
            ${disease_gene} > ${params.run_id}.phrank.txt
        """
}
