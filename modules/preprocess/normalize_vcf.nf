#!/usr/bin/env nextflow
nextflow.enable.dsl=2
// This process normalizes a VCF file by compressing it with bgzip and indexing it with tabix.
// It checks if the input VCF is already in BGZF format and skips compression if it is.
// If the input VCF is a symbolic link, it resolves the link to get the actual file path.
// The output is a compressed VCF file and its corresponding index file.
// The output files are named "input.vcf.gz" and "input.vcf.gz.tbi".
process NORMALIZE_VCF {
    input:
        path vcf

    output:
        path "input.vcf.gz", emit: vcf
        path "input.vcf.gz.tbi", emit: tbi

    script:
    """
    INPUT_VCF_TYPE="\$(file -b $vcf)"
    if echo "\${INPUT_VCF_TYPE}" | grep -q 'symbolic link to'; then
        SYM_LINK="\$(readlink -f $vcf)"
        INPUT_VCF_TYPE="\$(file -b \${SYM_LINK})"
    fi
    if echo "\${INPUT_VCF_TYPE}" | grep -q 'BGZF'; then
        echo "The file is in BGZF format, ready for tabix."
        cp $vcf input.vcf.gz
        tabix -p vcf input.vcf.gz
    else
        bgzip -c $vcf > input.vcf.gz
        tabix -p vcf input.vcf.gz
    fi
    """
    stub:
    """
    touch input.vcf.gz input.vcf.gz.tbi
    """
}
