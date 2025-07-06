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
workflow MAIN {
    // Reference assembly and file paths
    refAssembly             = params.ref_ver == 'hg19' ? 'grch37' : 'grch38'
    chrmap                  = "${params.ref_dir}/bcf_annotate/chrmap.txt"
    ref_loc                 = "${params.ref_dir}/phrank/${params.ref_ver}/${refAssembly}_symbol_to_location.txt"
    ref_to_sym              = "${params.ref_dir}/phrank/${params.ref_ver}/${params.ref_ver}/ensembl_to_symbol.txt"
    ref_sorted_sym          = "${params.ref_dir}/phrank/${params.ref_ver}/gene_to_symbol_sorted.txt"

    // Filter BED settings
    ref_exonic_filter_bed   = "${params.ref_dir}/filter_exonic/${params.ref_ver}.bed"
    ref_filter_bed          = params.bed_filter ? params.bed_filter : (params.exome_filter ? ref_exonic_filter_bed : '/dev/null')

    // Phrank files
    phrank_dagfile          = "${params.ref_dir}/phrank/${params.ref_ver}/child_to_parent.txt"
    phrank_disease_annotation = "${params.ref_dir}/phrank/${params.ref_ver}/disease_to_pheno.txt"
    phrank_gene_annotation  = "${params.ref_dir}/phrank/${params.ref_ver}/gene_to_phenotype.txt"
    phrank_disease_gene     = "${params.ref_dir}/phrank/${params.ref_ver}/disease_to_gene.txt"

    // OMIM/HPO annotation
    omim_hgmd_phen          = "${params.ref_dir}/omim_annotate/${params.ref_ver}/HGMD_phen.tsv"
    omim_obo                = "${params.ref_dir}/omim_annotate/hp.obo"
    omim_genemap2           = "${params.ref_dir}/omim_annotate/${params.ref_ver}/genemap2_v2022.rds"
    omim_pheno              = "${params.ref_dir}/omim_annotate/${params.ref_ver}/HPO_OMIM.tsv"

    // GNOMAD VCF paths
    ref_gnomad_genome       = "${params.ref_dir}/filter_vep/${params.ref_ver}/gnomad.${params.ref_ver}.blacklist.genomes.vcf.gz"
    ref_gnomad_genome_idx   = "${ref_gnomad_genome}.tbi"
    ref_gnomad_exome        = "${params.ref_dir}/filter_vep/${params.ref_ver}/gnomad.${params.ref_ver}.blacklist.exomes.vcf.gz"
    ref_gnomad_exome_idx    = "${ref_gnomad_exome}.tbi"

    // VEP plugin file names
    vep_dbnsfp_name         = params.ref_ver == 'hg19' ? 'dbNSFP4.3a_grch37.gz' : 'dbNSFP4.1a_grch38.gz'
    vep_gnomad_name         = params.ref_ver == 'hg19' ? 'gnomad.genomes.r2.1.sites.grch37_noVEP.vcf.gz' : 'gnomad.genomes.GRCh38.v3.1.2.sites.vcf.gz'
    vep_cadd_name           = params.ref_ver == 'hg19' ? 'hg19_whole_genome_SNVs.tsv.gz' : 'hg38_whole_genome_SNV.tsv.gz'

    // VEP cache and plugin paths
    vep_dir_cache           = "${params.ref_dir}/vep/${params.ref_ver}/"
    vep_dir_plugins         = "${params.ref_dir}/vep/${params.ref_ver}/Plugins/"
    vep_custom_gnomad       = "${vep_dir_cache}${vep_gnomad_name}"
    vep_custom_clinvar      = "${vep_dir_cache}clinvar_20220730.vcf.gz"
    vep_custom_hgmd         = "${vep_dir_cache}HGMD_Pro_2022.2_${params.ref_ver}.vcf.gz"
    vep_plugin_revel        = "${vep_dir_plugins}new_tabbed_revel_${refAssembly}.tsv.gz"
    vep_plugin_spliceai_snv = "${vep_dir_plugins}spliceai_scores.masked.snv.${params.ref_ver}.vcf.gz"
    vep_plugin_spliceai_indel = "${vep_dir_plugins}spliceai_scores.masked.indel.${params.ref_ver}.vcf.gz"
    vep_plugin_cadd         = "${vep_dir_plugins}${vep_cadd_name}"
    vep_plugin_dbnsfp       = "${vep_dir_plugins}${vep_dbnsfp_name}"
    vep_idx                 = "${params.ref_dir}/vep/${params.ref_ver}/*.tbi"

    // Pipeline directories
    ref_annot_dir           = "${params.ref_dir}/annotate"
    ref_var_tier_dir        = "${params.ref_dir}/var_tier"
    ref_merge_expand_dir    = "${params.ref_dir}/merge_expand"
    ref_mod5_diffusion_dir  = "${params.ref_dir}/mod5_diffusion"
    ref_predict_new_dir     = "${params.ref_dir}/predict_new"
    ref_model_inputs_dir    = "${params.ref_dir}/model_inputs"

    workflow.onComplete {
        println "[INFO] refAssembly: ${refAssembly}"
        println "[INFO] chrmapPath: ${chrmap}"
        println "[INFO] refLoc: ${ref_loc}"
        println "[INFO] refToSym: ${ref_to_sym}"
        println "[INFO] sortedSym: ${ref_sorted_sym}"
        println "[INFO] exonicBed: ${ref_exonic_filter_bed}"
        println "[INFO] filterBed: ${ref_filter_bed}"
        println "[INFO] vepCache: ${vep_dir_cache}"
        println "[INFO] vepPlugins: ${vep_dir_plugins}"
        println "[INFO] dbnsfpName: ${vep_plugin_dbnsfp}"
        println "[INFO] gnomadName: ${vep_custom_gnomad}"
        println "[INFO] caddName: ${vep_plugin_cadd}"
        println "[INFO] omimHgmdPhen: ${omim_hgmd_phen}"
        println "[INFO] omimOBO: ${omim_obo}"
        println "[INFO] omimGenemap2: ${omim_genemap2}"
        println "[INFO] omimPheno: ${omim_pheno}"
        println "[INFO] vepCustomClinvar: ${vep_custom_clinvar}"
        println "[INFO] vepCustomHgmd: ${vep_custom_hgmd}"
        println "[INFO] vepPluginRevel: ${vep_plugin_revel}"
        println "[INFO] vepPluginSpliceaiSNV: ${vep_plugin_spliceai_snv}"
        println "[INFO] vepPluginSpliceaiIndel: ${vep_plugin_spliceai_indel}"
        println "[INFO] vepIdx: ${vep_idx}"
    }

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
        vep_dir_cache,
        vep_dir_plugins,
        vep_custom_gnomad,
        vep_custom_clinvar,
        vep_custom_hgmd,
        vep_plugin_revel,
        vep_plugin_spliceai_snv,
        vep_plugin_spliceai_indel,
        vep_plugin_cadd,
        vep_plugin_dbnsfp,
        file(vep_idx)
    )
    HPO_SIM(
        params.input_hpo,
        omim_hgmd_phen,
        omim_obo,
        omim_genemap2,
        omim_pheno
    )
    ANNOTATE_BY_MODULES(
        ANNOTATE_BY_VEP.out.vep_output,
        HPO_SIM.out.hgmd_sim,
        HPO_SIM.out.omim_sim,
        file(ref_annot_dir),
    )
    PHRANK_SCORING(
        NORMALIZE_VCF.out.vcf
    )
    JOIN_TIER_PHRANK(
        ANNOTATE_BY_MODULES.out.scores,
        PHRANK_SCORING.out,
        file(ref_annot_dir),
        file(ref_var_tier_dir),
        file(ref_merge_expand_dir),
    )
    MERGE_SCORES_BY_CHROMOSOME(
        PHRANK_SCORING.out,
        JOIN_TIER_PHRANK.out.tier.collect(),
        JOIN_TIER_PHRANK.out.compressed_scores.collect(),
        file(ref_annot_dir),
        file(ref_mod5_diffusion_dir),
        file(ref_merge_expand_dir),
    )
    // Run Prediction on the final merged output
    PREDICTION(
        MERGE_SCORES_BY_CHROMOSOME.out.merged_matrix,
        MERGE_SCORES_BY_CHROMOSOME.out.merged_compressed_scores,
        file(ref_predict_new_dir),
        file(ref_model_inputs_dir),
    )
}
