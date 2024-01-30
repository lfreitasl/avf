process VCFTOOLS_BASE {
    tag "base_filters"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vcftools:0.1.16--he513fc3_4' :
        'biocontainers/vcftools:0.1.16--he513fc3_4' }"

    input:
    // We're gonna make some filters in the file to keep only biallelic variants (wthout indels). Optionally, you can use a parameter to thin the matrix
    // 1 snp each 1000 basepairs
    path variant_file
    val thin

    output:
    path "*.vcf"				      , emit: vcf
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def input_file = ("$variant_file".endsWith(".vcf")) ? "--vcf ${variant_file}" :
        ("$variant_file".endsWith(".vcf.gz")) ? "--gzvcf ${variant_file}" :
        ("$variant_file".endsWith(".bcf")) ? "--bcf ${variant_file}" : ''
    def out_name = thin ? "all_samples_filtered_thinned" : "all_samples_filtered"
    def thining = thin ? "--thin 1000" : ''

    """
    vcftools \\
        $input_file \\
        --out $out_name \\
        --min-alleles 2 \\
        --max-alleles 2 \\
        --remove-indels \\
        $thining \\
        --recode \\
        --recode-INFO-all \\
	

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcftools: \$(echo \$(vcftools --version 2>&1) | sed 's/^.*VCFtools (//;s/).*//')
    END_VERSIONS
    """
}
