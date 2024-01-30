#!/usr/bin/env nextflow

params.bfile_prefix = "/data/gpfs/projects/punim1484/ad/adni/gwas/ADNI_Omni25_AD_PACC"
params.outdir = "tmp"
 

workflow {
    chrs= Channel.from( 1..2 )

    plink_data = Channel
    .fromFilePairs("${params.bfile_prefix}.{bed,fam,bim}", size:3)
    .ifEmpty {error "No matching plink files"}
    
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

    // Accepts the output from chr_bfiles
    input:
    path bfiles 

    // Define the outputs of this process
    output:
    path "${bfiles[0].baseName}.afreq", emit: variant_freq_result

    """
    plink2 --bfile ${bfiles[0].baseName} --freq --out ${bfiles[0].baseName}
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



