#!/usr/bin/env nextflow

in_bfile   = Channel.fromFilePairs("${params.bfile}{.bed,.bim,.fam}", size:3)
in_cov     = Channel.fromPath(params.cov)
out_dir    = file("$PWD/GWAS_results")

in_bfile
    .into { grm_bfile; combined_bfile; male_bfile; female_bfile }

out_dir.mkdir()

process transform {
    tag { "Transform" }
    publishDir "$PWD/pheno_files", mode: 'copy', overwrite: true

    input:
    file(cov) from in_cov

    output:
    file("Combined.pheno") into combined_pheno
    file("MALES.pheno") into males_pheno
    file("FEMALES.pheno") into females_pheno
    
    """
    Rscript ${projectDir}/bin/Trans.R ${cov}
    """
}

process makeGRM {
    publishDir "$PWD/grm_files", mode: 'copy', overwrite: true
    tag { sample }

    input:
    set sample, file("*") from grm_bfile

    output:
    set sample, file("${sample}.grm.id"), file("${sample}.grm.bin"), file("${sample}.grm.N.bin") into grm_combined, grm_male, grm_female
  
    """
    gcta64 --bfile ${sample} --make-grm --out ${sample} --thread-num 10
    """
}

combined_bfile
    .flatten()
    .toList()
    .join(grm_combined, by:0)
    .set { input_gwas_combined }

male_bfile
    .flatten()
    .toList()
    .join(grm_male, by:0)
    .set { input_gwas_males }

female_bfile
    .flatten()
    .toList()
    .join(grm_female, by:0)
    .set { input_gwas_females }

process combinedGWAS {
    publishDir "$PWD/out_combined", mode: 'copy', overwrite: true
    tag { sample }
    
    input:
    set sample, file(bed), file(bim), file(fam), file(id), file(bin), file(Nbin) from input_gwas_combined
    file(pheno) from combined_pheno

    output:
    set sample, file("*") into gwas_combined

    """
    gcta64 --mlma --bfile ${sample} --grm ${sample} --pheno ${pheno} --out ${sample}_combined_bmi --thread-num 10
    gcta64 --mlma --bfile ${sample} --grm ${sample} --pheno ${pheno} --mpheno 2 --out ${sample}_combined_varbmi --thread-num 10
    """
}

process malesGWAS {
    publishDir "$PWD/out_males", mode: 'copy', overwrite: true
    tag { sample }
    
    input:
    set sample, file(bed), file(bim), file(fam), file(id), file(bin), file(Nbin) from input_gwas_males
    file(pheno) from males_pheno

    output:
    set sample, file("*") into gwas_males

    """
    gcta64 --mlma --bfile ${sample} --grm ${sample} --pheno ${pheno} --out ${sample}_males_bmi --thread-num 10
    gcta64 --mlma --bfile ${sample} --grm ${sample} --pheno ${pheno} --mpheno 2 --out ${sample}_males_varbmi --thread-num 10
    """
}

process femaleGWAS {
    publishDir "$PWD/out_females", mode: 'copy', overwrite: true
    tag { sample }
    
    input:
    set sample, file(bed), file(bim), file(fam), file(id), file(bin), file(Nbin) from input_gwas_females
    file(pheno) from females_pheno

    output:
    set sample, file("*") into gwas_females

    """
    gcta64 --mlma --bfile ${sample} --grm ${sample} --pheno ${pheno} --out ${sample}_females_bmi --thread-num 10
    gcta64 --mlma --bfile ${sample} --grm ${sample} --pheno ${pheno} --mpheno 2 --out ${sample}_females_varbmi --thread-num 10
    """
}
