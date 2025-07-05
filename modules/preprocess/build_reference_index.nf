#!/usr/bin/env nextflow
nextflow.enable.dsl=2
// This process builds a reference index for the specified genome version.
//// It downloads the genome sequence from UCSC, processes it to remove 'chr' prefixes,
// and creates a FASTA file along with its index and dictionary.
// The output files are named according to the genome version specified in the parameters.
process BUILD_REFERENCE_INDEX {
    container "broadinstitute/gatk"
    storeDir "${params.outdir}/general/reference_index/"

    output:
    path "final_${params.ref_ver}.fa", emit: fasta
    path "final_${params.ref_ver}.fa.fai", emit: fasta_index
    path "final_${params.ref_ver}.dict", emit: fasta_dict

    script:
    """
    wget --quiet http://hgdownload.soe.ucsc.edu/goldenPath/${params.ref_ver}/bigZips/${params.ref_ver}.fa.gz
    gunzip ${params.ref_ver}.fa.gz
    sed 's/>chr/>/g' ${params.ref_ver}.fa > num_prefix_${params.ref_ver}.fa
    samtools faidx num_prefix_${params.ref_ver}.fa 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M > final_${params.ref_ver}.fa
    samtools faidx final_${params.ref_ver}.fa
    gatk CreateSequenceDictionary -R final_${params.ref_ver}.fa
    """
}