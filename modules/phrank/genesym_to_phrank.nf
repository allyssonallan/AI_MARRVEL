#!/usr/bin/env nextflow
nextflow.enable.dsl=2
// This process converts gene symbols to Phrank format for disease-gene association analysis.
// It takes gene symbols, HPO terms, a DAG file, disease annotations, gene annotations, and disease-gene associations as inputs.
// The output is a Phrank formatted file with the specified run ID.
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
    run_phrank.py \\
        $gene $hpo $dagfile $disease_annotation $gene_annotation $disease_gene > ${params.run_id}.phrank.txt
    """
}