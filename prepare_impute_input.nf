#!/usr/bin/env nextflow

params.bfile_prefix = "/data/gpfs/projects/punim1484/ad/adni/gwas/ADNI_Omni25_AD_PACC"
params.outdir = "tmp"
params.refpanel="/data/gpfs/projects/punim1484/ukb/genomics/project/2024_NF_GwasPipeline/HRC.r1-1.GRCh37.wgs.mac5.sites.tab"
params.token = file('token_topmed.txt').text.trim()
params.jobname = "bwgoudey, AD_US"
params.build = "hg19"

workflow {

    /*println "Token: ${params.token}"
    println "Job Name: ${params.jobname}"
    println "Build: ${params.build}"
    println "Build: ${params.bfile_prefix}"
    */

    chrs= Channel.from( 1..2 )

    plink_data = Channel
    .fromFilePairs("${params.bfile_prefix}.{bed,fam,bim}", size:3)
    .ifEmpty {error "No matching plink files"}

    
    plink_data.view { "Plink data: $it" }
      
    chr_bfiles = splitChr(plink_data.combine(chrs))
    freq_file = computeVariantFreq(chr_bfiles)
    snp_fix_results = runWillRaynorScript(chr_bfiles, freq_file, params.refpanel)
    bfile_update=plinkUpdateSNPs(chr_bfiles,snp_fix_results)
    vcf_file=bfileToVcf(bfile_update)
    
    vcf_files = vcf_file.collect()
    
    sendToTopMed(params.token, vcf_files, params.jobname, params.build)
    
}



/*
 * Split a fasta file into multiple files
 */
process splitChr {
    module 'FlexiBLAS/3.2.0'
    module 'PLINK/2.00a3.6'
    publishDir "./tmp", mode: 'symlink'

    input:
    tuple val(input_bfile), path(files), val(chr)
 
    output:
    path "${input_bfile}_chr${chr}.{bed,fam,bim}", emit: chr_bfiles
 
    """
    plink2 --bfile $input_bfile \
          --chr $chr \
          --make-bed \
          --out ${input_bfile}_chr${chr} 
    """
}


/*
 * Output variamt frequencies from a given PLINK file
 */
process computeVariantFreq {
    module 'PLINK/1.9b_6.21-x86_64'
    //publishDir "./tmp", mode: 'symlink'
    // Accepts the output from chr_bfiles
    input:
    path bfiles 

    // Define the outputs of this process
    output:
    path "${bfiles[0].baseName}.frq", emit: variant_freq_result

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
      path bfile
      path freq_file
      path ref_panel
    output: 
      path "{Force-Allele1,Strand-Flip,Exclude,ID,LOG,Chromosome,Position}-${bfile[0].baseName}-HRC.txt", emit: snp_fix_results

    """
    perl /home/bwgoudey/tools/HRC-1000G-check-bim-v4.2.7/HRC-1000G-check-bim.pl -b ${bfile[0].baseName}.bim \
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
      path bfile
      path snp_fix_results
      
    output: 
      path "${bfile[0].baseName}-updated.{bim,bed,fam}"

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
    publishDir "./tmp", mode: 'symlink'
    module 'PLINK/2.00a3.6'
    input:                                                                                                                                                                               
      path bfile                                                                                                                                                                         
                                                                                                                                                                                         
    output:                                                                                                                                                                              
      path "${bfile[0].baseName}.vcf", emit: vcf_file

    """
    plink2 --bfile ${bfile[0].baseName} --recode vcf --out ${bfile[0].baseName}
    """
}


process sendToTopMed {
    publishDir "./tmp", mode: 'symlink'
    input: 
      val token
      path vcf_files
      val jobname
      val inital_build
      //tuple val(token), path(vcfs), val(job_name), val(inital_build)

    output:
        path "tmp.txt"

    script:
      def curlFiles = vcf_files.collect { file -> "-F \"files=@$file\"" }.join(' ')
    
    """
    echo 'curl https://imputation.biodatacatalyst.nhlbi.nih.gov/api/v2/jobs/submit/imputationserver \
      -X "POST" \
      -H "X-Auth-Token: ${token}" \
      -F "mode=qconly" \
      -F "job-name=${jobname}" \
      $curlFiles \
      -F "refpanel=apps@topmed-r3" \
      -F "build=${inital_build}" \
      -F "phasing=eagle" \
      -F "population=all" \
      -F "meta=yes"' > tmp.txt
    """
}


