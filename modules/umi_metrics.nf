process UMI_METRICS {
    tag "UMI metrics - $sample_id"

    // Creates a series of metrics files:
    // **.family_sizes.txt**: metrics on the frequency of different types of families of different sizes
    // **.duplex_family_sizes.txt**: metrics on the frequency of duplex tag families by the number of observations from each strand
    // **.duplex_yield_metrics.txt**: summary QC metrics produced using 5%, 10%, 15%...100% of the data
    // **.umi_counts.txt**: metrics on the frequency of observations of UMIs within reads and tag families
    // **.duplex_qc.pdf**: a series of plots generated from the preceding metrics files for visualization
    // **.duplex_umi_counts.txt**: (optional) metrics on the frequency of observations of duplex UMIs within reads and tag families. This file is only produced _if_ the `--duplex-umi-counts` option is used as it requires significantly more memory to track all pairs of UMIs seen when a large number of UMI sequences are present.
    
    // Within the metrics files the prefixes CS, SS and DS are used to mean:

    // CS: tag families where membership is defined solely on matching genome coordinates and strand
    // SS: single-stranded tag families where membership is defined by genome coordinates, strand and UMI; ie. 50/A and 50/B are considered different tag families.
    // DS: double-stranded tag families where membership is collapsed across single-stranded tag families from the same double-stranded source molecule; i.e. 50/A and 50/B become one family

    // cpus = 8
    // memory = { 32.GB * task.attempt }
    // time = { 6.hour * task.attempt }

    publishDir "${params.outputDir}/log/umi_metrics/logs/", pattern: "${sample_id}.*.txt", mode:'copy'
    publishDir "${params.outputDir}/log/umi_metrics/plots/", pattern: "${sample_id}.*.pdf", mode:'copy'

    input:
    tuple val(sample_id), path(bam_file)

    output:
    path "${sample_id}.family_sizes.txt"
    path "${sample_id}.duplex_family_sizes.txt"
    path "${sample_id}.duplex_yield_metrics.txt"
    path "${sample_id}.umi_counts.txt"
    path "${sample_id}.duplex_qc.pdf"

    script:
    """
    fgbio CollectDuplexSeqMetrics   \
        --input=${bam_file}         \
        --output="${sample_id}"
    """
    
    stub:
    """
    touch "${sample_id}.family_sizes.txt"
    touch "${sample_id}.duplex_family_sizes.txt"
    touch "${sample_id}.duplex_yield_metrics.txt"
    touch "${sample_id}.umi_counts.txt"
    touch "${sample_id}.duplex_qc.pdf"
    """
}