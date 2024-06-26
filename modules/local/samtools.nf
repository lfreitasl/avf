process SAMTOOLS_MPILEUP {
    tag "$meta.id"
    label 'process_single'

conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.18--h50ea8bc_1' :
        'biocontainers/samtools:1.18--h50ea8bc_1' }"

    input:
    tuple val(meta), path(input), path(intervals)
    path  fasta

    output:
    tuple val(meta), path("*.mpileup.gz"), emit: mpileup
    path  "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def intervals = intervals ? "-l ${intervals}" : ""
    
    """
    samtools mpileup -Ou -b $bamlist -f reference/EC20.fna | bcftools call -mv --ploidy 1 -Ob -o test_all_together.bcf
    
	cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
