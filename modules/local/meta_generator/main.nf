process META_VCF {
    tag "$meta"
    label 'process_single'
   // errorStrategy 'ignore'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-7cbdd89e4a230442dcff2455c5693faa81d0cda4:f36d68574b255379a8cabacc75a9ec9c783aff6e-0':
        'biocontainers/mulled-v2-7cbdd89e4a230442dcff2455c5693faa81d0cda4:f36d68574b255379a8cabacc75a9ec9c783aff6e-0' }"

    input:
    path  meta
    path  vcf
    path  str

    output:
    tuple path('vcfs_info.csv'), path('*sorted_meta.csv'), path(str), emit: vcf_meta
    path "versions.yml"                                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/hybrider/bin/
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${vcf.baseName}"

    """
    report_vcfs_meta.R \\
        $meta \\
        $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        r-dartr: \$(Rscript -e "library(dartR); cat(as.character(packageVersion('dartR')))")
        r-vcfr: \$(Rscript -e "library(vcfR); cat(as.character(packageVersion('vcfR')))")
    END_VERSIONS
    """
}
