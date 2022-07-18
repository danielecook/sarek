include { CAT_CAT as CAT_MPILEUP         } from '../../../../modules/nf-core/modules/cat/cat/main'
include { SAMTOOLS_MPILEUP               } from '../../../../modules/nf-core/modules/samtools/mpileup/main'

workflow RUN_MPILEUP {
    take:
        cram                    // channel: [mandatory] [meta, cram, interval]
        fasta                   // channel: [mandatory]

    main:

    ch_versions = Channel.empty()

    SAMTOOLS_MPILEUP(cram, fasta)
    mpileup = SAMTOOLS_MPILEUP.out.mpileup.branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }

    //Merge mpileup only when intervals and natural order sort them
    CAT_MPILEUP(mpileup.intervals
        .map{ meta, pileup ->
            new_meta = meta.tumor_id ? [patient:meta.patient, normal_id:meta.normal_id, tumor_id:meta.tumor_id, sex:meta.sex, id:meta.tumor_id + "_vs_" + meta.normal_id, num_intervals:meta.num_intervals] // not annotated, so no variantcaller necessary
                                        : [patient:meta.patient, sample:meta.sample, status:meta.status, sex:meta.sex, id:meta.sample, num_intervals:meta.num_intervals]
            [groupKey(new_meta, meta.num_intervals), pileup]
            }
        .groupTuple(sort:true))


    ch_versions = ch_versions.mix(SAMTOOLS_MPILEUP.out.versions)
    ch_versions = ch_versions.mix(CAT_MPILEUP.out.versions)

    emit:
    versions = ch_versions
    mpileup = Channel.empty().mix(CAT_MPILEUP.out.file_out, mpileup.no_intervals)
}