process BCFTOOLS_MPILEUP {
    tag "all_samples"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.18--h8b25389_0':
        'biocontainers/bcftools:1.18--h8b25389_0' }"

    input:
   // tuple val(meta), path(bam)
    val tmp
    path fasta
    val names


    output:
    path "*vcf.gz"       , emit: vcf
    path "*vcf.gz.tbi"   , emit: tbi
    path "*stats.txt"    , emit: stats
    path  "versions.yml" , emit: versions
    path  "*.list"       , emit: listas

    when:
    task.ext.when == null || task.ext.when

    script:
    // def prefix = task.ext.prefix ?: "${meta}"

    """
    for b in \$(echo "${names}" | tr ',' '\n' | tr -d '[' | tr -d ']'); do
        echo "\$b" >> samplenames.list
    done       

    for b in \$(echo "${tmp}" | tr ',' '\n' | tr -d '[' | tr -d ']'); do
        echo "\$b" >> bam.list
    done

    bcftools \\
        mpileup \\
        --fasta-ref $fasta \\
        -b bam.list \\
        | bcftools call -mv --ploidy 1 -O v \\
        | bcftools reheader -s samplenames.list \\
        | bcftools view -O z -o all_samps.vcf.gz
    

    tabix -p vcf -f all_samps.vcf.gz

    bcftools stats all_samps.vcf.gz > all_samps.bcftools_stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
