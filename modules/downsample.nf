process DOWNSAMPLE {
    tag "$sample_id"
    // Downsample FastQ file to number of million reads.
    // Twist kit target region: 
    // 150x = 2.5M reads
    // 9600x = 167M reads
    // Approximately 40% is on target - so cap set at 400M reads.
    
    // cpus = 8
    // memory = { 32.GB * task.attempt }
    // time = { 6.hour * task.attempt }

    input:
    tuple val(sample_id), path(reads, stageAs: 'raw/*')
    val n_reads

    output:
    tuple val(sample_id), path("${sample_id}_R{1,2}.fq.gz"), emit: sampled_reads

    script:
    def total_reads = n_reads*1000000
    """
    reformat.sh \
        in1=${reads[0]} \
        in2=${reads[1]} \
        out1="${sample_id}_R1.fq.gz" \
        out2="${sample_id}_R2.fq.gz" \
        sampleseed=1 \
        samplereadstarget=${total_reads} \
    """
    
    stub:
    """
    touch "${sample_id}_R1.fq.gz"
    touch "${sample_id}_R2.fq.gz"
    """
}