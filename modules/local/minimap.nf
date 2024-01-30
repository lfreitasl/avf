process ALIGNMENT {
    tag "$meta"
    label 'process_medium'

	conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:365b17b986c1a60c1b82c6066a9345f38317b763-0' :
        'biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:365b17b986c1a60c1b82c6066a9345f38317b763-0' }"

    input:
    tuple val(meta), path(genomes)
    path reference_fasta

    output:
    tuple val(meta), path("*.bam"),  emit: bam
    path "*.bam",		     emit: tmp
    val meta,			     emit: names

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta}"

    """
     minimap2 -ax asm5 --cs $reference_fasta $genomes |
	samtools sort |
	samtools view -@ ${task.cpus} -b -h -o ${prefix}.bam 

    """
}
