#!/usr/bin/env nextflow

params.bfile_prefix = "/data/gpfs/projects/punim1484/ad/adni/gwas/ADNI_Omni25_AD_PACC"
params.outdir = "tmp"
params.refpanel="/data/gpfs/projects/punim1484/ukb/genomics/project/2024_NF_GwasPipeline/HRC.r1-1.GRCh37.wgs.mac5.sites.tab"
params.token_file = 'token_topmed.txt'
params.token = file(params.token_file).text.trim()
params.jobname = "bwgoudey, WTCCC_US"
params.build = "hg19"
params.HRC_check_exe = "/home/bwgoudey/tools/HRC-1000G-check-bim-v4.2.13_NoReadKey/HRC-1000G-check-bim-NoReadKey.pl"
params.publish_dir = "./tmp"

workflow {

    /*println "Token: ${params.token}"
    println "Job Name: ${params.jobname}"
    println "Build: ${params.build}"
    println "Build: ${params.bfile_prefix}"
    */

    chrs= Channel.from( 1..22 )

    plink_data = Channel
    .fromFilePairs("${params.bfile_prefix}.{bed,fam,bim}", size:3)
    .ifEmpty {error "No matching plink files"}

    plink_data.view { "Plink data: $it" }

    chr_bfile_tuple = chrs.combine(plink_data)
    chr_bfile = splitChr(chr_bfile_tuple)
    freq_file = computeVariantFreq(chr_bfile)
    
    snp_fix_results = runWillRaynorScript(chr_bfile.join(freq_file, by: 0), params.refpanel)
    bfile_update = plinkUpdateSNPs(chr_bfile.join(snp_fix_results, by: 0))

    vcf_file = bfileToVcf(bfile_update)
    vcf_files = vcf_file.map { chr, file -> file }.collect()
    sendToTopMed(params.token, vcf_file, params.jobname, params.build)
    
    
}



/*
 * Split a fasta file into multiple files
 */
process splitChr {
    module 'FlexiBLAS/3.2.0'
    module 'PLINK/2.00a3.6'
    publishDir params.publish_dir, mode: 'symlink'

    input:
    tuple val(chr), val(bfile_prefix), path(bfiles)
 
    output:
    tuple val(chr), path("${bfile_prefix}_chr${chr}.{bed,fam,bim}")//, emit: chr_bfiles
 
    """
    plink2 --bfile $bfile_prefix \
          --chr $chr \
          --make-bed \
          --out ${bfile_prefix}_chr${chr} 
    """
}


/*
 * Output variamt frequencies from a given PLINK file
 */
process computeVariantFreq {
    module 'PLINK/1.9b_6.21-x86_64'
    publishDir params.publish_dir, mode: 'symlink'
    // Accepts the output from chr_bfiles
    input:
    tuple val(chr), path(bfiles)

    // Define the outputs of this process
    output:
    tuple val(chr), path("${bfiles[0].baseName}.frq")// emit: variant_freq_result

    """
    plink --bfile ${bfiles[0].baseName} --freq --out ${bfiles[0].baseName}
    """
}
/*
 * Run Will Rayner toolbox for SNP checking
 * This is really a SNP filtering stage
 * TODO: this could be broken up into multiple steps
 */
process runWillRaynorScript {
    input: 
      tuple val(chr), path(bfile), path(freq_file)
      path ref_panel

    output: 
      tuple val(chr), path("{Force-Allele1,Strand-Flip,Exclude,ID,LOG,Chromosome,Position}-${bfile[0].baseName}-HRC.txt")// emit: snp_fix_results

    """
    perl ${params.HRC_check_exe} -b ${bfile[0].baseName}.bim \
                -f $freq_file \
                -r $ref_panel \
                -h \
                -t=0.3
    """

}

/*
 * Run Will Rayner toolbox for SNP checking
 * This is really a SNP filtering stage
 * TODO: this could be broken up into multiple steps
 */
process plinkUpdateSNPs {
    module 'PLINK/1.9b_6.21-x86_64'
    input: 
      tuple val(chr), path(bfile), path(snp_fix_results)
      
    output: 
      tuple val(chr), path("${bfile[0].baseName}-updated.{bim,bed,fam}"), path("Force-Allele1-${bfile[0].baseName}-HRC.txt")

    """
    plink --bfile ${bfile[0].baseName} --exclude Exclude-${bfile[0].baseName}-HRC.txt --make-bed --out TEMP1
    plink --bfile TEMP1 --update-map Chromosome-${bfile[0].baseName}-HRC.txt --update-chr --make-bed --out TEMP2
    plink --bfile TEMP2 --update-map Position-${bfile[0].baseName}-HRC.txt --make-bed --out TEMP3
    plink --bfile TEMP3 --flip Strand-Flip-${bfile[0].baseName}-HRC.txt --make-bed --out TEMP4
    plink --bfile TEMP4 --reference-allele Force-Allele1-${bfile[0].baseName}-HRC.txt --make-bed --out ${bfile[0].baseName}-updated
    """

}

/*                                                                                                                                                                                       
 * Run Will Rayner toolbox for SNP checking                                                                                                                                              
 * This is really a SNP filtering stage                                                                                                                                                  
 * TODO: this could be broken up into multiple steps                                                                                                                                     
 */                                                                                                                                                                                      
process bfileToVcf {
    module "FlexiBLAS/3.2.0"
    publishDir params.publish_dir, mode: 'symlink'
    //module 'PLINK/2.00a3.6'
    module 'PLINK/1.9b_6.21-x86_64'
    module 'OpenSSL/1.1'
    module 'bzip2/1.0.8'
    module 'BCFtools/1.15.1'
    
    input:
      tuple val(chr), path(bfile), path(force_allele_file)

    output:                                                                                                                                                                              
      tuple val(chr), path("${bfile[0].baseName}.vcf.gz")// emit: vcf_file

    """
    echo $force_allele_file
    plink --bfile ${bfile[0].baseName} \
           --recode vcf \
           --out ${bfile[0].baseName} \
            --a2-allele ${force_allele_file}

    bcftools sort ${bfile[0].baseName}.vcf \
        -Oz -o ${bfile[0].baseName}.vcf.gz
    """
}


process sendToTopMed {
    publishDir params.publish_dir, mode: 'symlink'
    input: 
      val token
      tuple val(chr), path(vcf_file)
      val jobname
      val inital_build
      //tuple val(token), path(vcfs), val(job_name), val(inital_build)

    output:
        path "tmp_${chr}.txt"

    //script:
    //  def curlFiles = vcf_files.collect { file -> "-F \"files=@$file\"" }.join(' ')
    
    """
    echo 'curl https://imputation.biodatacatalyst.nhlbi.nih.gov/api/v2/jobs/submit/imputationserver \
      -X "POST" \
      -H "X-Auth-Token: ${token}" \
      -F "mode=qconly" \
      -F "job-name=${jobname}-${chr}" \
      -F "files=@${vcf_file}" \
      -F "refpanel=apps@topmed-r3" \
      -F "build=${inital_build}" \
      -F "phasing=eagle" \
      -F "population=all" \
      -F "meta=yes"' > tmp_${chr}.txt
    """
}


