process INTERVALS {
    tag "BedFileToIntervalList"

    input:
    path reference_genome
    path bedfile

    output:
    path "target_region.interval_list", emit: target_list

    script:
    def db = file(params.reference_genome).getName() + ".fna"
    """
    gatk BedToIntervalList \
				 I=${bedfile} \
				 O="target_region.interval_list" \
				 SD=${db}
    """

    stub:
    """
    touch "target_region.interval_list"
    """
}
