//
// DeepSomatic tumor-normal variant calling
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { GATK4_MERGEVCFS as MERGE_DEEPSOMATIC_VCF } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { DEEPSOMATIC                         } from '../../../modules/nf-core/deepsomatic/main'

workflow BAM_VARIANT_CALLING_SOMATIC_DEEPSOMATIC {
    take:
    cram          // channel: [mandatory] [ meta, normal_cram, normal_crai, tumor_cram, tumor_crai ]
    dict          // channel: [optional]  [ meta, dict ]
    fasta         // channel: [mandatory] [ fasta ]
    fasta_fai     // channel: [mandatory] [ fasta_fai ]
    intervals     // channel: [mandatory] [ intervals, num_intervals ] or [ [], 0 ] if no intervals

    main:
    versions = Channel.empty()

    cram_intervals = cram.combine(intervals)
        // Move num_intervals to meta map
        .map{ meta, normal_cram, normal_crai, tumor_cram, tumor_crai, intervals, num_intervals -> [ meta + [ num_intervals:num_intervals ], normal_cram, normal_crai, tumor_cram, tumor_crai, intervals ] }

    DEEPSOMATIC(cram_intervals, fasta.map{ fasta -> [ [ id:fasta.baseName ], fasta ] }, fasta_fai.map{ fasta_fai -> [ [ id:fasta_fai.baseName ], fasta_fai ] }, [ [ id:'null' ], [] ] )

    // Figuring out if there is one or more vcf(s) from the same sample
    vcf_out = DEEPSOMATIC.out.vcf.branch{
        // Use meta.num_intervals to asses number of intervals
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // Only when using intervals
    vcfs_to_merge = vcf_out.intervals.map{ meta, vcf -> [ groupKey(meta, meta.num_intervals), vcf ]}.groupTuple()
    MERGE_DEEPSOMATIC_VCF(vcfs_to_merge, dict)

    // Mix intervals and no_intervals channels together
    vcf = Channel.empty().mix(MERGE_DEEPSOMATIC_VCF.out.vcf, vcf_out.no_intervals)
        // add variantcaller to meta map and remove no longer necessary field: num_intervals
        .map{ meta, vcf -> [ meta - meta.subMap('num_intervals') + [ variantcaller:'deepsomatic' ], vcf ] }

    versions = versions.mix(DEEPSOMATIC.out.versions)

    emit:
    vcf

    versions
}
