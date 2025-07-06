#!/usr/bin/env nextflow
nextflow.enable.dsl=2
include { VCF_TO_VARIANTS } from '../modules/variant/vcf_to_variants.nf'
include { VARIANTS_TO_ENSEMBL } from '../modules/variant/variant_to_ensembl.nf'
include { ENSEMBL_TO_GENESYM } from '../modules/variant/ensembl_to_genesym.nf'
include { GENESYM_TO_PHRANK } from '../modules/phrank/genesym_to_phrank.nf'

/*
========================================================================================
    PHRANK SCORING SUBWORKFLOW
========================================================================================
*/
workflow PHRANK_SCORING {
    take:
        vcf

    main:
        VCF_TO_VARIANTS(vcf)
        VARIANTS_TO_ENSEMBL(VCF_TO_VARIANTS.variants, params.ref_loc)
        ENSEMBL_TO_GENESYM(VARIANTS_TO_ENSEMBL.ensembl, params.ref_to_sym, params.ref_sorted_sym)
        GENESYM_TO_PHRANK(
            ENSEMBL_TO_GENESYM.genesym,
            params.input_hpo,
            params.phrank_dagfile,
            params.phrank_disease_annotation,
            params.phrank_gene_annotation,
            params.phrank_disease_gene
        )

    emit:
        phrank = GENESYM_TO_PHRANK.out
}