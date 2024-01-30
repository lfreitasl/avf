process VCFTOOLS_MISSING {
    tag "missingness_filters"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vcftools:0.1.16--he513fc3_4' :
        'biocontainers/vcftools:0.1.16--he513fc3_4' }"

    input:
    // We're gonna make some filters in the file to keep only biallelic variants (wthout indels). Optionally, you can use a parameter to thin the m>    // 1 snp each 1000 basepairs

    path variant_file
    val  maxmissing

    output:
    path "*.vcf"                                      , emit: vcf
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def input_file = ("$variant_file".endsWith(".vcf")) ? "--vcf ${variant_file}" :
        ("$variant_file".endsWith(".vcf.gz")) ? "--gzvcf ${variant_file}" :
        ("$variant_file".endsWith(".bcf")) ? "--bcf ${variant_file}" : ''
    def out_name = ("$variant_file".endsWith("thinned.recode.vcf")) ? "all_samples_filtered_thinned_$maxmissing_Comp" : "all_samples_filtered_$maxmissing"
    
    """
	vcftools \\
	$input_file \\
	 --out $out_name \\
        --max-missing $maxmissing \\
        --recode \\
        --recode-INFO-all \\


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcftools: \$(echo \$(vcftools --version 2>&1) | sed 's/^.*VCFtools (//;s/).*//')
    END_VERSIONS
    """
}	
