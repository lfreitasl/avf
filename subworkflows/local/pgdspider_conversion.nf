include { PGDSPIDER } from '../../modules/local/conversor.nf'

workflow VCF_TO_STR {

    take:
    vcfpath
    thin
    maxmissing

    main:

	PGDSPIDER(
	vcfpath,
	thin,
	maxmissing
	)


    emit:
	str = PGDSPIDER.out.str

}
