include { ALIGNMENT } from '../../modules/local/minimap.nf'

workflow ALIGNMENT_MINIMAP2 {
    take:
    fasta
    reference_fasta
    

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    ALIGNMENT (
        fasta,
        reference_fasta
    )
    

    emit:
    alignments = ALIGNMENT.out.bam
    bamnames = ALIGNMENT.out.tmp
    sampnames = ALIGNMENT.out.names
}
