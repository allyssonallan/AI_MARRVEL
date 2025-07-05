#!/usr/bin/env nextflow
nextflow.enable.dsl=2
include { CONVERT_GVCF } from '../modules/preprocess/convert_gvcf.nf'
include { FILTER_UNPASSED } from '../modules/preprocess/filter_unpassed.nf'
include { FILTER_MITO_AND_UNKNOWN_CHR } from '../modules/preprocess/filter_mito_and_unknown_chr.nf'
include { FILTER_PROBAND } from '../modules/preprocess/filter_proband.nf'
include { FILTER_BED } from '../modules/preprocess/filter_bed.nf'

workflow VCF_PRE_PROCESS {
    take:
    vcf
    tbi
    fasta
    fasta_index
    fasta_dict

    main:
    CONVERT_GVCF(
        vcf,
        tbi,
        fasta,
        fasta_index,
        fasta_dict,
        params.chrmap
    )
    FILTER_UNPASSED(
        CONVERT_GVCF.out.vcf,
        CONVERT_GVCF.out.tbi,
        params.chrmap
    )
    FILTER_MITO_AND_UNKNOWN_CHR(
        FILTER_UNPASSED.out.vcf,
        FILTER_UNPASSED.out.tbi,
    )
    FILTER_BED(
        FILTER_MITO_AND_UNKNOWN_CHR.out.vcf,
        FILTER_MITO_AND_UNKNOWN_CHR.out.tbi,
        moduleDir.resolve(params.ref_filter_bed),
    )
    FILTER_PROBAND(
        FILTER_BED.out.vcf,
        FILTER_BED.out.tbi,
        params.ref_gnomad_genome,
        params.ref_gnomad_genome_idx,
        params.ref_gnomad_exome,
        params.ref_gnomad_exome_idx
    )

    emit:
    vcf = FILTER_PROBAND.out.vcf
}