#!/usr/bin/env nextflow

params.bfile_prefix = "/data/gpfs/projects/punim1484/ad/adni/gwas/ADNI_Omni25_AD_PACC"
params.outdir = "tmp"
 

workflow {
    chrs= Channel.from( 1..2 )

    plink_data = Channel
    .fromFilePairs("${params.bfile_prefix}.{bed,fam,bim}", size:3)
    .ifEmpty {error "No matching plink files"}
    

    chrs.subscribe { println "value: $it" }
    plink_data.view()
    chr_bfiles = splitChr(plink_data.combine(chrs))
    computeVariantFreq(chr_bfiles)
}



/*
 * Split a fasta file into multiple files
 */
process splitChr {
    module 'PLINK/2.00a3.6'

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
    module 'PLINK/2.00a3.6'
    input: 
    path bfile
    
    output: 
    path "${bfile}_freq", emit: chr_bfiles

    """
    plink2 --bfile $bfile  --freq --out ${bfile}_freq
    """
}

/*
 * Run Will Rayner toolbox for SNP checking
 * This is really a SNP filtering stage
 * TODO: this could be broken up into multiple steps
 *
process computeVariantFreq {
    input: 
    
    output: 

    """
    perl HRC-1000G-check-bim.pl -b <bim file> 
                -f <freq-file> 
                -r <reference> 
                -h
    """
}
*/



