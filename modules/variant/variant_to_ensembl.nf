#!/usr/bin/env nextflow
nextflow.enable.dsl=2
// This process converts a VCF file to a simple variant format with chromosome, position, reference, and alternate alleles.
process VARIANTS_TO_ENSEMBL {
    input:
    path var
    path ref

    output:
    path "${params.run_id}-ensmbl.txt"

    script:
    """
    location_to_gene.py $var $ref | \\
     sed 's/:/\\t/g' | sed 's/X\\t/23\\t/g' | sed 's/Y\\t/24\\t/g' | \\
     sed 's/MT\\t/25\\t/g' > ${params.run_id}-ensmbl.txt
    """
}