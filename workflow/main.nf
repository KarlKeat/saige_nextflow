#!/usr/bin/env nextflow
workflow {
    phenos_ch = Channel
        .fromPath(params.tsv_input)
        .splitCsv(header: true, sep: '\t') // Parses rows as Map objects using the first line as keys
        .map { row -> 
            // Transform columns (e.g., turn strings into file paths)
            tuple(row.pheno_name, row.ancestry, row.pheno_covar_file, row.covar_list, row.qcovar_list)
        }

    raw_imputed_ch = Channel.fromPath(params.pmbb_imputed_chunks)
        .map { file -> file.parent / file.baseName }

    plink2_input_ch = phenos_ch.combine(raw_imputed_ch)
    PLINK2_PREPROCESS(plink2_input_ch)
    
    chrom_ch = Channel.of(1..22)
    merge_input_ch = PLINK2_PREPROCESS.out.plink_preprocessed_chunk
        .groupTuple(by: [0, 1, 2, 3, 4])
        .map { pheno_name, ancestry, pheno_covar_file, covar_list, qcovar_list, file_lists ->
            [pheno_name, ancestry, pheno_covar_file, covar_list, qcovar_list, file_lists.flatten()]
        }
        .combine(chrom_ch)  // collect all chunks sharing the same pheno_name + ancestry

    PLINK2_MERGE(merge_input_ch)
    plink2_merged_chroms = PLINK2_MERGE.out.plink_preprocessed_chromosomes.groupTuple(by: [0, 1, 2, 3, 4])
        .map { pheno_name, ancestry, pheno_covar_file, covar_list, qcovar_list, file_lists ->
            [pheno_name, ancestry, pheno_covar_file, covar_list, qcovar_list, file_lists.flatten()]
        } 
    
    MAKE_STEP1_INPUT(plink2_merged_chroms)
    SAIGE_STEP1(MAKE_STEP1_INPUT.out.step1_in)

    step2_in_ch = plink2_merged_chroms.join(SAIGE_STEP1.out.step1_out, by: [0, 1, 2, 3, 4])
    SAIGE_STEP2(step2_in_ch)
}

process PLINK2_PREPROCESS {
    //publishDir "${params.outputDir}/$pheno_name/preprocessed_by_chunk/"
    errorStrategy { task.exitStatus == 13 ? 'ignore' : 'terminate' }
    input: 
        tuple val(pheno_name), val(ancestry), val(pheno_covar_file), val(covar_list), val(qcovar_list), val(pmbb_imputed_prefix)

    output:
        tuple val(pheno_name), val(ancestry), val(pheno_covar_file), val(covar_list), val(qcovar_list), path("*.{psam,pvar,pgen}"), emit: plink_preprocessed_chunk

    script:
        """ 
        ${projectDir}/scripts/preprocess_imputed.sh $pmbb_imputed_prefix $ancestry $pheno_covar_file
        """
}

process PLINK2_MERGE {
    //publishDir "${params.outputDir}/$pheno_name/preprocessed_by_chrom/"
    input: 
        tuple val(pheno_name), val(ancestry), val(pheno_covar_file), val(covar_list), val(qcovar_list), path("chunks/*"), val(chrom)

    output:
        tuple val(pheno_name), val(ancestry), val(pheno_covar_file), val(covar_list), val(qcovar_list), path("*.{psam,pvar,pgen}"), emit: plink_preprocessed_chromosomes

    script:
        """ 
        ${projectDir}/scripts/merge_preprocessed_chunks.sh $pheno_name $ancestry $chrom chunks/
        """
}

process MAKE_STEP1_INPUT {
    //publishDir "${params.outputDir}/$pheno_name/step1_in/"
    input: 
        tuple val(pheno_name), val(ancestry), val(pheno_covar_file), val(covar_list), val(qcovar_list), path("chroms/*")

    output:
        tuple val(pheno_name), val(ancestry), val(pheno_covar_file), val(covar_list), val(qcovar_list), path(ld_pruned_merged_output), emit: step1_in

    script:
        ld_pruned_merged_output = pheno_name + "_" + ancestry + "_pruned.*"
        """ 
        ${projectDir}/scripts/make_step1_input.sh $pheno_name $ancestry chroms/
        """
}

process SAIGE_STEP1 {
    publishDir "${params.outputDir}/$pheno_name/step1_out/", mode: 'copy'
    input: 
        tuple val(pheno_name), val(ancestry), val(pheno_covar_file), val(covar_list), val(qcovar_list), path("ld_pruned_imputed/*")

    output:
        tuple val(pheno_name), val(ancestry), val(pheno_covar_file), val(covar_list), val(qcovar_list), path(grm_file), path(varianceRatio_file), emit: step1_out

    script:
        grm_file = "saige_step1_"+ pheno_name + "_" + ancestry + ".rda"
        varianceRatio_file = "saige_step1_"+ pheno_name + "_" + ancestry + ".varianceRatio.txt"
        """ 
        ${projectDir}/scripts/saige_step1.sh $pheno_name $ancestry $pheno_covar_file $covar_list $qcovar_list ld_pruned_imputed/${pheno_name}_${ancestry}_pruned ${params.trait_type} ${params.inv_normalize} ${task.cpus}
        """
}

process SAIGE_STEP2 {
    publishDir "${params.outputDir}/$pheno_name/step2_out/", mode: 'copy'
    input: 
        tuple val(pheno_name), val(ancestry), val(pheno_covar_file), val(covar_list), val(qcovar_list), path("chroms/*"), path(grm_file), path(varianceRatio_file)

    output:
        path "*.txt", emit: step2_out

    script:
        """ 
        ${projectDir}/scripts/saige_step2.sh $pheno_name $ancestry $pheno_covar_file $grm_file $varianceRatio_file chroms/
        """
}

