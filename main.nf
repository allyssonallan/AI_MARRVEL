#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * AI-MARRVEL: Main Workflow with Subworkflow (DSL2)
 * - All logic is inside workflow/process blocks
 * - Uses subworkflow for modularity and clarity
 * - No DSL1 constructs
 */

// Parameter definitions
// compute every derived variable once, in plain Groovy, using only ternaries
params.input_vcf = null
params.input_hpo = null
params.ref_dir = null
params.ref_ver = null
params.usage_file = 'USAGE.md'
params.version = false
params.help = false

// (moved variable assignments into workflow block below)


// Utility function for parameter validation
def showUsage() {
    if (params.help) {
        def helpFile = file(params.usage_file)  // Specify your Markdown file path here
        if (helpFile.exists()) {
            println helpFile.text
        } else {
            println """
            Sorry something went wrong, usage file not found!
            Please vist our website for more info : https://ai-marrvel.readthedocs.io/en/latest/
            """
        }
        exit 0
    }
}

def showVersion() {
    if (!params.version) {
        return
    }

    println "1.0.0"
    exit 0
}

def validateInputParams() {
    def checkPathParamMap = [
        "input_vcf": params.input_vcf,
        "input_hpo": params.input_hpo,
        "ref_dir"  : params.ref_dir,
        "ref_ver"  : params.ref_ver,
    ]

    checkPathParamMap.each { paramName, paramValue ->
        if (paramValue) {
            // Check if the file exists
            if(!(paramName == "ref_ver")) {
                def fileObj = file(paramValue, checkIfExists: true)
                //  println("Check file: '--${paramName}' value '${paramValue}' ")

                // Check the file type based on the parameter name
                if (paramName == "input_vcf" && !(paramValue.endsWith(".vcf") || paramValue.endsWith(".vcf.gz"))) {
                    println("Error: '--${paramName}' value '${paramValue}' should be a VCF file (.vcf) or (.vcf.gz)")
                    println("To see usage and available parameters run `nextflow run main.nf --help`")
                    exit 1
                } else if (paramName == "input_hpo" && !(paramValue.endsWith(".hpo") || paramValue.endsWith(".txt"))) {
                    println("Error: '--${paramName}' value '${paramValue}' should be an HPO file (.hpo) or (.txt)")
                    println("To see usage and available parameters run `nextflow run main.nf --help`")
                    exit 1
                } else if (paramName == "ref_dir" && !fileObj.isDirectory()) {
                    println("Error: '--${paramName}' value '${paramValue}' should be an directory.")
                    println("To see usage and available parameters run `nextflow run main.nf --help`")
                    exit 1
                }
            }

            if (paramName == "ref_ver" && !(paramValue.equals("hg19") || paramValue.equals("hg38")) ) { 
                println("Error: '--${paramName}' value ${paramValue} should be either set to 'hg19' or 'hg38'.")
                println("To see usage and available parameters run `nextflow run main.nf --help`")
                exit 1
            }

        } else {
            println("Input parameter '${paramName}' not specified or is null!")
            println("To see usage and available parameters run `nextflow run main.nf --help`")
            exit 1
        }
    }
}

// Process: Preprocessing
include { BUILD_REFERENCE_INDEX } from './modules/preprocess/build_reference_index.nf'
include { CONVERT_GVCF } from './modules/preprocess/convert_gvcf.nf'
include { FILTER_BED } from './modules/preprocess/filter_bed.nf'
include { FILTER_MITO_AND_UNKNOWN_CHR } from './modules/preprocess/filter_mito_and_unknown_chr.nf'
include { FILTER_PROBAND } from './modules/preprocess/filter_proband.nf'
include { FILTER_UNPASSED } from './modules/preprocess/filter_unpassed.nf'
include { NORMALIZE_VCF } from './modules/preprocess/normalize_vcf.nf'
// Process: HPO and Phrank
include { HPO_SIM } from './modules/phrank/hpo_sim.nf'
include { GENESYM_TO_PHRANK } from './modules/phrank/genesym_to_phrank.nf'
// Process: Variant
include { ENSEMBL_TO_GENESYM } from './modules/variant/ensembl_to_genesym.nf'
include { VARIANTS_TO_ENSEMBL } from './modules/variant/variant_to_ensembl.nf'
include { VCF_TO_VARIANTS } from './modules/variant/vcf_to_variants.nf'
// Process: Annotation
include { ANNOTATE_BY_MODULES } from './modules/annotation/annotate_by_modules.nf'
include { ANNOTATE_BY_VEP } from './modules/annotation/annotate_by_vep.nf'
include { JOIN_TIER_PHRANK } from './modules/annotation/join_tier_phrank.nf'
include { SPLIT_VCF_BY_CHROMOSOME } from './modules/annotation/split_vcf_by_chromosome.nf'
// Process: Merge
include { MERGE_SCORES_BY_CHROMOSOME } from './modules/merge/merge_scores_by_chromosome.nf'
include { PREDICTION } from './modules/merge/prediction.nf'
// Subworkflows
include { VCF_PRE_PROCESS } from './subworkflows/vcf_preprocess.nf'
include { PHRANK_SCORING } from './subworkflows/phrank_scoring.nf'

// Main workflow
workflow {
    // Move variable assignments here
    def refAssembly   = params.ref_ver == 'hg19' ? 'grch37' : 'grch38'
    def chrmapPath    = "${params.ref_dir}/bcf_annotate/chrmap.txt"
    def refLoc        = "${params.ref_dir}/phrank/${params.ref_ver}/${refAssembly}_symbol_to_location.txt"
    def refToSym      = "${params.ref_dir}/phrank/${params.ref_ver}/ensembl_to_symbol.txt"
    def sortedSym     = "${params.ref_dir}/phrank/${params.ref_ver}/gene_to_symbol_sorted.txt"
    def exonicBed     = "${params.ref_dir}/filter_exonic/${params.ref_ver}.bed"
    def filterBed     = params.bed_filter ? params.bed_filter : (params.exome_filter ? exonicBed : '/dev/null'

    )
    def vepCache      = "${params.ref_dir}/vep/${params.ref_ver}/"
    def vepPlugins    = "${vepCache}Plugins/"
    def dbnsfpName    = params.ref_ver == 'hg19' ? 'dbNSFP4.3a_grch37.gz' : 'dbNSFP4.1a_grch38.gz'
    def gnomadName    = params.ref_ver == 'hg19' ? 'gnomad.genomes.r2.1.sites.grch37_noVEP.vcf.gz' : 'gnomad.genomes.GRCh38.v3.1.2.sites.vcf.gz'
    def caddName      = params.ref_ver == 'hg19' ? 'hg19_whole_genome_SNVs.tsv.gz' : 'hg38_whole_genome_SNV.tsv.gz'

    showUsage()
    showVersion()
    validateInputParams()
    NORMALIZE_VCF(params.input_vcf)
    BUILD_REFERENCE_INDEX()
    VCF_PRE_PROCESS(
        NORMALIZE_VCF.out.vcf,
        NORMALIZE_VCF.out.tbi,
        BUILD_REFERENCE_INDEX.out.fasta,
        BUILD_REFERENCE_INDEX.out.fasta_index,
        BUILD_REFERENCE_INDEX.out.fasta_dict
    )
    SPLIT_VCF_BY_CHROMOSOME(VCF_PRE_PROCESS.out.vcf)
    ANNOTATE_BY_VEP(
        SPLIT_VCF_BY_CHROMOSOME.out.chr_vcfs.flatten(),
        params.vep_dir_cache,
        params.vep_dir_plugins,
        params.vep_custom_gnomad,
        params.vep_custom_clinvar,
        params.vep_custom_hgmd,
        params.vep_plugin_revel,
        params.vep_plugin_spliceai_snv,
        params.vep_plugin_spliceai_indel,
        params.vep_plugin_cadd,
        params.vep_plugin_dbnsfp,
        file(params.vep_idx)
    )
    HPO_SIM(
        params.input_hpo,
        params.omim_hgmd_phen,
        params.omim_obo,
        params.omim_genemap2,
        params.omim_pheno
    )
    ANNOTATE_BY_MODULES(
        ANNOTATE_BY_VEP.out.vep_output,
        HPO_SIM.out.hgmd_sim,
        HPO_SIM.out.omim_sim,
        file(params.ref_annot_dir),
    )
    PHRANK_SCORING(
        NORMALIZE_VCF.out.vcf
    )
    JOIN_TIER_PHRANK(
        ANNOTATE_BY_MODULES.out.scores,
        PHRANK_SCORING.out,
        file(params.ref_annot_dir),
        file(params.ref_var_tier_dir),
        file(params.ref_merge_expand_dir),
    )
    MERGE_SCORES_BY_CHROMOSOME(
        PHRANK_SCORING.out,
        JOIN_TIER_PHRANK.out.tier.collect(),
        JOIN_TIER_PHRANK.out.compressed_scores.collect(),
        file(params.ref_annot_dir),
        file(params.ref_mod5_diffusion_dir),
        file(params.ref_merge_expand_dir),
    )
    // Run Prediction on the final merged output
    PREDICTION(
        MERGE_SCORES_BY_CHROMOSOME.out.merged_matrix,
        MERGE_SCORES_BY_CHROMOSOME.out.merged_compressed_scores,
        file(params.ref_predict_new_dir),
        file(params.ref_model_inputs_dir),
    )
}
