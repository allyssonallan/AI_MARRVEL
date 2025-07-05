#!/usr/bin/env nextflow
nextflow.enable.dsl=2
// This process filters a VCF file based on a reference BED file.
process FILTER_BED {
    input:
    path vcf
    path tbi
    path ref_filter_bed

    output:
    path "${params.run_id}.recode.vcf.gz", emit: vcf
    path "${params.run_id}.recode.vcf.gz.tbi", emit: tbi

    script:
    """
    if [ -f ${ref_filter_bed} ]; then
        awk '{gsub(/^chr/, ""); print}' ${ref_filter_bed} > bed
        bcftools filter --regions-file bed ${vcf} -Oz -o "${params.run_id}.recode.vcf.gz"
        tabix -p vcf "${params.run_id}.recode.vcf.gz"
    else
        cp ${vcf} "${params.run_id}.recode.vcf.gz"
        cp ${tbi} "${params.run_id}.recode.vcf.gz.tbi"
    fi
    """
}