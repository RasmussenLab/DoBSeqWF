process GENOMICSDB {
    tag "MERGE_GVCF"
    // Merge gVCF files into a single GenomicsDB
    
    // cpus = 8
    // memory = { 32.GB * task.attempt }
    // time = { 6.hour * task.attempt }

    publishDir "${params.outputDir}/log/", pattern: "genomicsdb.log", mode:'copy'
    publishDir "${params.outputDir}/", pattern: "./genomicsdb", mode:'copy'

    input:
    path vcfs
    path index
    path bedfile

    output:
    path "genomicsdb", emit: gendb
    path "genomicsdb.log"

    script:
    input_files_command = vcfs.collect(){"--variant ${it}"}.join(' ')
    """
    gatk --java-options "-Xmx8g -XX:-UsePerfData"   \
        GenomicsDBImport                            \
        $input_files_command                        \
        --genomicsdb-workspace-path genomicsdb      \
        --intervals ${bedfile}                      \
        --tmp-dir .                                 \
        > >(tee -a "genomicsdb.log")                \
        2> >(tee -a "genomicsdb.log" >&2)
    """

    stub:
    """
    touch "genomicsdb" "genomicsdb.log"
    """
}

