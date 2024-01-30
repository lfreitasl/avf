include { BCFTOOLS_MPILEUP } from '../../modules/nf-core/bcftools/mpileup/main'

workflow BCFTOOLS_VARIANT_CALL {

	take:
	tmp
	reference
	sampnames

	main:
	ch_versions = Channel.empty()

	BCFTOOLS_MPILEUP(tmp, reference, sampnames)
	
	// ch_versions = ch_versions.mix(BCFTOOLS_MPILEUP.out.versions.first())


	emit:
	vcf = BCFTOOLS_MPILEUP.out.vcf
	versions = ch_versions
}
