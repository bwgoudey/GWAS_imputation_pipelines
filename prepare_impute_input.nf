#!/usr/bin/env nextflow

//params.bfile_prefix = "/data/gpfs/projects/punim1484/ad/adni/gwas/ADNI_Omni25_AD_PACC"
params.bfile_prefixes = ["/path/to/set1","/path/to/set2"] // Example
params.outdir = "tmp"
params.refpanel="/data/gpfs/projects/punim1484/ukb/genomics/project/2024_NF_GwasPipeline/HRC.r1-1.GRCh37.wgs.mac5.sites.tab"
params.token_file = 'token_topmed.txt'
params.token = file(params.token_file).text.trim()
params.jobname = "bwgoudey, WTCCC_US"
params.build = "hg19"
params.HRC_check_exe = "/home/bwgoudey/tools/HRC-1000G-check-bim-v4.2.13_NoReadKey/HRC-1000G-check-bim-NoReadKey.pl"
params.publish_dir = "./tmp"
params.chrStart = 1  // Default start chromosome
params.chrEnd = 22   // Default end chromosome
params.vcf_mode = "single" // "single" or "multiple", default is "single"


workflow {

    /*println "Token: ${params.token}"
    println "Job Name: ${params.jobname}"
    println "Build: ${params.build}"
    println "Build: ${params.bfile_prefix}"
    */

    chrs = Channel.from(params.chrStart..params.chrEnd)
    
    plink_data_sets= Channel                                                        
    .fromFilePairs("${params.bfile_prefixes}.{bed,fam,bim}", size:3)              
    .ifEmpty {error "No matching plink files"} 
    
    //plink_data_sets.toList().view { "Plink data: $it" }
    
    chr_bfile_tuple = chrs.combine(plink_data_sets)
    
    //If using all chromsomes in a given file
    if(params.chrStart == 0 && params.chrEnd == 0) {
        chr_bfile = chr_bfile_tuple
    } else {
        chr_bfile = splitChr(chr_bfile_tuple)
    }
    ///*
    freq_file = computeVariantFreq(chr_bfile)
    
    snp_fix_results = runWillRaynorScript(chr_bfile.join(freq_file, by: 0), params.refpanel)
    /*chr_bfile.toList().view { "Chr_bfile: $it" }
    snp_fix_results.toList().view { "Snp_fix_results data2: $it" }
    chr_bfile.join(snp_fix_results, by: [0,1]).toList().view { "Joined: $it" }
    */
    
    bfile_update = plinkUpdateSNPs(chr_bfile.join(snp_fix_results, by: [0,1], failOnMismatch: true, failOnDuplicate: true))
    chr_clean = removeRemappedChr(bfile_update)
    vcf_file = bfileToVcf(chr_clean)
    vcf_files = vcf_file.map { chr, file -> tuple(-1,file) }.groupTuple(by: 0)//vcf_file.map { chr, file -> file }.collect()
    //println "Build: ${params.build}"          
    //println "Mode: ${params.vcf_mode}"          
    if (params.vcf_mode == "single") {
        sendToTopMed(params.token, vcf_file, params.jobname, params.build, params.vcf_mode)
    } else if (params.vcf_mode == "multiple") {
        sendToTopMed(params.token, vcf_files, params.jobname, params.build, params.vcf_mode)
    }
    //*/
    
    
}



/*
 * Split a fasta file into multiple files
 */
process splitChr {
    module 'FlexiBLAS/3.2.0'
    module 'PLINK/2.00a3.6'
    //publishDir params.publish_dir, mode: 'symlink'

    input:
    tuple val(chr), val(bfile_prefix), path(bfiles)
 

    output:
    tuple val(chr), val(bfile_prefix), path("${bfile_prefix}_chr${chr}.{bed,fam,bim}")//, emit: chr_bfiles
 
    """
    plink2 --bfile ${bfile_prefix} \
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
    //publishDir params.publish_dir, mode: 'symlink'
    // Accepts the output from chr_bfiles
    input:
    tuple val(chr), val(bfile_prefix), path(bfiles)

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
      tuple val(chr), val(bfile_prefix), path(bfile), path(freq_file)
      path ref_panel

    output: 
      tuple val(chr), val(bfile_prefix), path("{Force-Allele1,Strand-Flip,Exclude,Chromosome,Position}-${bfile[0].baseName}-HRC.txt")// emit: snp_fix_results

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
      tuple val(chr), val(bfile_prefix), path(bfile), path(snp_fix_results)
      
    output: 
      tuple val(chr), path("${bfile[0].baseName}-updated.{bim,bed,fam}"), path("Force-Allele1-${bfile[0].baseName}-HRC.txt")

    /*"""
    plink --bfile ${bfile[0].baseName} --exclude Exclude-${bfile[0].baseName}-HRC.txt --make-bed --out TEMP1
    plink --bfile TEMP1 --update-map Chromosome-${bfile[0].baseName}-HRC.txt --update-chr --make-bed --out TEMP2
    plink --bfile TEMP2 --update-map Position-${bfile[0].baseName}-HRC.txt --make-bed --out TEMP3
    plink --bfile TEMP3 --flip Strand-Flip-${bfile[0].baseName}-HRC.txt --make-bed --out TEMP4
    plink --bfile TEMP4 --reference-allele Force-Allele1-${bfile[0].baseName}-HRC.txt --make-bed --out ${bfile[0].baseName}-updated
    """
    */


    // Define paths from snp_fix_results tuple
    // Assuming tuple structure:
    //{Force-Allele1,Strand-Flip,Exclude,ID,LOG,Chromosome,Position}
    // 0 - Chromosome file
    // 1 - Exclude file
    // 2 - Force Allele file
    // 3 - Position file
    // 444 Strand-Flip file
    """
    plink --bfile ${bfile[0].baseName} --exclude ${snp_fix_results[1]} --make-bed --out TEMP1
    plink --bfile TEMP1 --update-map ${snp_fix_results[0]} --update-chr --make-bed --out TEMP2
    plink --bfile TEMP2 --update-map ${snp_fix_results[3]} --make-bed --out TEMP3
    plink --bfile TEMP3 --flip ${snp_fix_results[4]} --make-bed --out TEMP4
    plink --bfile TEMP4 --reference-allele ${snp_fix_results[2]} --make-bed --out ${bfile[0].baseName}-updated
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
      tuple val(chr), path(vcf_files)
      val jobname
      val initial_build
      val vcf_mode
      //tuple val(token), path(vcfs), val(job_name), val(inital_build)

    output:
        path("${jobname}-${chr}*.txt")

    script:
    // Preparing a string for CURL file upload based on the mode
    def curlFiles = ""
    def jobName_full = "${jobname}-${chr}"
    def outfile = "curl_file_"
    if (vcf_mode == "single") {
        curlFiles = "-F \"files=@${vcf_files}\"" // Assuming vcf_files is a single file in this mode
        jobName_full= "${jobName_full}_${vcf_files.baseName}"
    } else if (vcf_mode == "multiple") {
        curlFiles = vcf_files.collect { file -> "-F \"files=@$file\"" }.join(' ')
        jobName_full= "${jobName_full}_all_chr"
    }

    """
    echo 'curl https://imputation.biodatacatalyst.nhlbi.nih.gov/api/v2/jobs/submit/imputationserver \
    -X "POST" \\
    -H "X-Auth-Token: ${token}" \\
    -F "mode=qconly" \\
    -F "job-name=${jobName_full}" \\
    $curlFiles \\
    -F "refpanel=apps@topmed-r3" \\
    -F "build=${initial_build}" \\
    -F "phasing=eagle" \\
    -F "population=all" \\
    -F "meta=yes"' > ${jobName_full}.txt
    """
}



/*
 * Output variamt frequencies from a given PLINK file
 */
process removeRemappedChr {
    module 'PLINK/1.9b_6.21-x86_64'
    // Accepts the output from chr_bfiles
    input:
    tuple val(chr), path(bfiles), path(force_allele_file)

    // Define the outputs of this process
    output:
    tuple val(chr), path("${bfiles[0].baseName}_filt.{bed,fam,bim}"), path(force_allele_file)

    """
    plink --bfile ${bfiles[0].baseName} --chr ${chr} --out ${bfiles[0].baseName}_filt --make-bed
    """
}
