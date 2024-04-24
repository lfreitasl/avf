include { META_VCF          } from '../../modules/local/meta_generator/main'
include { PLOT_CLUSTERING   } from '../../modules/local/plot_clustering/main'
include { STRUCTURE         } from '../../modules/local/structure/main'

workflow STR_RUN_PLOT {

	take:
	meta
	vcf
	str

	main:
	ch_rep_per_k  	= Channel.empty()
	ch_kvalue     	= Channel.empty()
	ch_vcf_meta     = Channel.empty()
	ch_str_in    	= Channel.empty()
	ch_str		    = Channel.empty()
	ch_ffiles		= Channel.empty()
	ch_qfiles		= Channel.empty()
	ch_versions     = Channel.empty()

	META_VCF(meta, vcf, str)

	ch_vcf_meta     = ch_str.mix(META_VCF.out.vcf_meta.ifEmpty([]))
	ch_versions     = ch_versions.mix(META_VCF.out.versions.first().ifEmpty([]))
    // Transforming csv into a channel
	ch_vcf_meta
	.map { meta, sampmeta, vcfs ->
		return [ meta.splitCsv( header:true, sep:',' ), sampmeta, vcfs ]
	}
	.map { meta, sampmeta, vcfs ->
		return [create_csv_channel(meta[0]), sampmeta, vcfs] 
	}
	.set { ch_vcf_meta }

	ch_kvalue       = ch_kvalue.mix(Channel.from(1..params.k_value))
	ch_rep_per_k    = ch_rep_per_k.mix(Channel.from(1..params.rep_per_k))
	ch_str_in       = ch_str_in.mix(ch_vcf_meta.combine(ch_kvalue).combine(ch_rep_per_k))

	STRUCTURE(
		ch_str_in,
		params.noadmix,
		params.freqscorr,
		params.inferalpha,
		params.alpha,
		params.inferlambda,
		params.lambda,
		params.ploidy,
		params.burnin,
		params.mcmc
	)

	ch_ffiles       = ch_ffiles.mix(STRUCTURE.out.ffiles.groupTuple().map{meta,sampmeta,ffiles->return [meta,sampmeta[0],ffiles]}ifEmpty([]))
	ch_qfiles       = ch_qfiles.mix(STRUCTURE.out.qfiles.groupTuple().map{meta,sampmeta,ffiles->return [meta,sampmeta[0],ffiles]}ifEmpty([]))
	ch_versions     = ch_versions.mix(STRUCTURE.out.versions.first().ifEmpty([]))

	PLOT_CLUSTERING(
		ch_ffiles,
		params.plot_admix,
		params.plot_str,
		params.popinfo,
		params.writecsv
	)
	
	// ch_versions = ch_versions.mix(BCFTOOLS_MPILEUP.out.versions.first())


	emit:
	// vcf = BCFTOOLS_MPILEUP.out.vcf
	versions = ch_versions
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_csv_channel(LinkedHashMap row) {
    // create meta map
	def meta = [:]
	meta.id     = row.filenames
	meta.n_inds = row.n_inds
	meta.n_loc  = row.n_locs

    // add path(s) of the fastq file(s) to the meta map
	def csv_meta = meta 
	return csv_meta
}
