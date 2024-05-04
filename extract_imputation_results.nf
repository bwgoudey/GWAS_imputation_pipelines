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

params.csvfile = '/data/gpfs/projects/punim1484/wtc/imputation/wtc_impute_zips.csv' // Path to CSV file
params.imputation_outdir = "/data/gpfs/projects/punim1484/wtc/imputation"
params.ukb_path = '../../ukb/genomics/impute/' 
params.info_thresh = 0.5
params.password = "alongpassword"
params.ukb_proj_id="ukb56110"


workflow {

    zipped_files = Channel
        .fromPath(params.csvfile)
        .splitCsv(header: false, sep: ',', strip: true)
        .map { row -> tuple(row[0], row[1]) }

    vcf_files = UnzipImputed(zipped_files, params.bfile_prefix)//, params.password)
    indexed_vcf_files = IndexVcf(vcf_files)//, params.password)
    filt_vcf  = FilterVcf(indexed_vcf_files, params.info_thresh, params.ukb_path, params.ukb_proj_id)
    plink_files = ConvertToPlink(filt_vcf)
    //plink_files_aligned = AlignToUKB(plink_files)

}

process IndexVcf {
    module 'OpenSSL/1.1'
    module 'bzip2/1.0.8'
    module 'BCFtools/1.15.1'

    input:
        path vcf
    
    output:
        path("${renamed_vcf_base}.dose.vcf.{gz,gz.csi}")
    
    script:
        """
        bcftools index ${vcf}
        """
}

process UnzipImputed {
    input:
        tuple path(zip_file), val(password) 
        path bfile_prefix
    
    output:
        path "${bfile_prefix.baseName}${vcf_base}.dose.vcf.gz" 
    
    script:
        def vcf_base = zip_file.baseName.replace('_', '')
        def renamed_vcf_base= "${bfile_prefix.baseName}${vcf_base}.dose.vcf.gz"
        """
        unzip -P "$password" -j $zip_file "${vcf_base}.dose.vcf.gz" 
        mv ${vcf_base}.dose.vcf.gz ${renamed_vcf_base}
        """
}

process FilterVcf {
    input:
        path vcf_file
        val threshold
        path ukb_path
        val ukb_proj_id
    
    output:
        path "${vcf_file[0].baseName}.filtered.vcf.gz" 
    
    script:
        """
        bcftools filter -e 'INFO/info<=${threshold}' ${vcf_file[0].basename}.vcf.gz -Oz -o ${vcf_file[0].baseName}.filtered.vcf.gz
        """
}


process ConvertToPlink {
    input:
        path vcf_file 
    
    output:
        path "${vcf_file.baseName}" 
    
    script:
        """
        plink --vcf $vcf_file --make-bed --out ${vcf_file.baseName}
        """
}
