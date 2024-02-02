!/usr/bin/env nextflow

params.bfile_prefix = "/data/gpfs/projects/punim1484/ad/adni/gwas/ADNI_Omni25_AD_PACC"
params.outdir = "tmp"
params.refpanel=/HRC.r1-1.GRCh37.wgs.mac5.sites.tab

workflow {
    chrs= Channel.from( 1..2 )

    plink_data = Channel
    .fromFilePairs("${params.bfile_prefix}.{bed,fam,bim}", size:3)
    .ifEmpty {error "No matching plink files"}
    
    chr_bfiles = splitChr(plink_data.combine(chrs))
    freq_file = computeVariantFreq(chr_bfiles)
    runWillRaynorScript(chr_bfiles, freq_file, params.refpanel)
}



/*
 * Split a fasta file into multiple files
 */
process splitChr {
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
    publishDir "./tmp", mode: 'symlink'
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
      path "${bfiles[0].baseName}-updated.{bed,fam,bim}", emit: update_bfiles

    """
    perl ~tools/HRC-1000G-check-bim.pl -b $bfile
                -f $bfile
                -r $freq_file
                -r $ref_panel
                -h 
                -t=0.3
 '   """

}
