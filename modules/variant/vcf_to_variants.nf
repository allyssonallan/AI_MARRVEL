#!/usr/bin/env nextflow
nextflow.enable.dsl=2
// This process converts a VCF file to a simple variant format with chromosome, position, reference, and alternate alleles.
process VCF_TO_VARIANTS {
    input:
    path vcf

    output:
    path "${params.run_id}-var-filt.txt", emit: var

    script:
    """
    zcat $vcf | awk 'substr(\$0, 1, 1) != "#"' | cut -f1,2,4,5 | sed 's/\t/:/g' > ${params.run_id}-var.txt
    cat ${params.run_id}-var.txt | sed 's/chr//g' | sort -u > ${params.run_id}-var-filt.txt
    """
}