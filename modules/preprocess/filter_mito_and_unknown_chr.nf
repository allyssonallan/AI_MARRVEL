#!/usr/bin/env nextflow
nextflow.enable.dsl=2
// This process filters the VCF file to remove mitochondrial and unknown chromosomes.
// It retains only the standard chromosomes (1-22, X, Y) and outputs a gzipped VCF file with an index.
process FILTER_MITO_AND_UNKNOWN_CHR {
    publishDir "${params.outdir}/${params.run_id}/vcf/", mode: 'copy'
    input:
    path vcf
    path tbi

    output:
    path "${params.run_id}.filt.rmMT.vcf.gz", emit: vcf
    path "${params.run_id}.filt.rmMT.vcf.gz.tbi", emit: tbi


    script:
    """
    bcftools view -r 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y $vcf -o ${params.run_id}.filt.rmMT.vcf

    bgzip ${params.run_id}.filt.rmMT.vcf
    tabix ${params.run_id}.filt.rmMT.vcf.gz
    """
}