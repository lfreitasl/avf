process BCFTOOLS_MPILEUP {
    tag "$meta"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.18--h8b25389_0':
        'biocontainers/bcftools:1.18--h8b25389_0' }"

    input:
    tuple val(meta), path(bam)
    val tmp // bam paths to .collect()
    path fasta


    output:
//    tuple val(meta), path("*vcf.gz")     , emit: vcf
 //   tuple val(meta), path("*vcf.gz.tbi") , emit: tbi
  //  tuple val(meta), path("*stats.txt")  , emit: stats
   // path  "versions.yml"                 , emit: versions
    path "*.list"			 , emit: listas
    
    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta}"

    """
    echo "${meta}" > sample_name.list

    
    for b in \$(echo "${tmp}" | tr ',' '\n' | tr -d '[' | tr -d ']'); do
        echo "\$b" >> bam.list
    done
    """
}

