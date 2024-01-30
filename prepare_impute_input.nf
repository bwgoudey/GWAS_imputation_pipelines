#!/usr/bin/env nextflow

params.bfile = "./gwas/ADNI_Omni25_AD_PACC"
params.outdir = "tmp"
 

workflow {
    channel.from
    chrs= Channel.from( 1..22 )
    chrs.subscribe { println "value: $it" }
    
    splitChr(params.bfile, chrs)
    computeVariantFreq(chr_bfiles)
}

/*
 * Split a fasta file into multiple files
 */
process splitChr {
 
    input:
    path input_bfile
    val chr
 
    output:
    path "${input_bfile}_chr${chr}", emit: chr_bfiles
 
    """
    echo plink --bfile $input_bfile 
          --chr $chr 
          --make-bed
          --out ${input_bfile}_chr${chr} 
    """
}


/*
 * Output variamt frequencies from a given PLINK file
 */
process computeVariantFreq {
    input: 
    path bfile
    
    output: 
    path "${bfile}_freq", emit: chr_bfiles

    """
    plink --bfile $bfile  --freq --out ${bfile}_freq
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



