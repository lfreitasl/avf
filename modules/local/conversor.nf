process PGDSPIDER {
    tag "converting_files"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://lfreitasl/pgdspider:latest':
        'docker.io/lfreitasl/pgdspider:latest' }"

    input:
    // We're gonna make some filters in the file to keep only biallelic variants (wthout indels). Optionally, you can use a parameter to thin the m>    // 1 snp each 1000 basepairs

    path variant_file
    val thin
    val maxmissing

    output:
    path "strformat_*"                                      , emit: str
    
    when:
    task.ext.when == null || task.ext.when

    script:
    
    def out_name = thin ? "strformat_filtered_thinned_$maxmissing_Comp" : "strformat_filtered_$maxmissing"

    """
        java \\
        -Xmx1024m \\
        -Xms512m \\
        -jar /bin/PGDSpider2-cli.jar \\
        -inputfile $variant_file \\
        -inputformat VCF \\
	-outputfile $out_name \\
        -outputformat STRUCTURE \\
        -spid /template_VCF_STRUCTURE.spid
    """

}
