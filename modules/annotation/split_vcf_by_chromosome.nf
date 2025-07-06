#!/usr/bin/env nextflow
nextflow.enable.dsl=2
process SPLIT_VCF_BY_CHROMOSOME {
    input:
        path vcf 

    output:
        path "chr*.vcf.gz", emit: chr_vcfs

    script:
        """
        # Get the list of chromosomes from the VCF file

        bgzip ${vcf}
        bcftools index ${vcf}.gz

        bcftools query -f '%CHROM\n' ${vcf}.gz | sort | uniq > chrom_list.txt
        # Split the VCF file by chromosome
        while read chrom; do
            bcftools view -r \${chrom} ${vcf}.gz -Oz -o chr\${chrom}.vcf.gz
        done < chrom_list.txt
        """
}