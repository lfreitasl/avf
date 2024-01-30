include { VCFTOOLS_BASE } from '../../modules/nf-core/vcftools/main'
include { VCFTOOLS_MISSING } from '../../modules/local/vcftools_missing.nf'

workflow VCFTOOLS_FILTER {
    take:
    vcf
    thining
    maxmiss

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    VCFTOOLS_BASE (
        vcf,
        thining
    )

    VCFTOOLS_MISSING (
        VCFTOOLS_BASE.out.vcf,
        maxmiss
    ) 

    emit:
    noncomp = VCFTOOLS_MISSING.out.vcf
    comp    = VCFTOOLS_BASE.out.vcf
}
