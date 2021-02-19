nextflow.enable.dsl=2

params.bam_list = "all_bam_samples.tsv"
params.maf = "consensus.maf"
params.fasta_base = "/juno/work/ci/resources/genomes/GRCh37/fasta/b37" // fasta
params.output_dir = "output"


num_threads = [1,2,4,8,12,24,32]
// number of samples per batch
limits = [1,2,4,8,16,32,48]
reps = [1,2,3]
process GBCMS_one_sample {
    tag "${threads}"
    // echo true
    maxForks 1
    cpus "${threads}"
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    tuple path(bam), path(bai), val(sampleId), path(maf), path(ref_fasta), file(fasta_indexes: "*")
    each threads

    output:
    path("${output_file}")

    script:
    output_file = "output.${threads}.maf"
    bam_arg = "${sampleId}:${bam}"
    """
    GetBaseCountsMultiSample \
    --omaf \
    --maq 20 \
    --baq 20 \
    --filter_improper_pair 0 \
    --output "${output_file}" \
    --thread "${threads}" \
    --maf "${maf}" \
    --bam "${bam_arg}" \
    --fasta "${ref_fasta}"
    """
}


process GBCMS_multi_sample {
    tag "${threads}:${limit}:${bams.size()}"
    // maxForks 3
    cpus "${threads}"
    publishDir "${params.output_dir}", mode: 'copy'
    // echo true

    input:
    tuple file(bams: '*'), file(bais: '*'), path(maf), path(ref_fasta), file(fasta_indexes: "*")
    each threads
    each limit
    each rep

    output:
    path("${output_file}")

    script:
    bam_arg = ''
    bams.eachWithIndex{ path, i ->
        if( i <= limit){
            bam_id = "${path}".replaceFirst(/.bam$/, "")
            bam_arg = bam_arg + ' --bam ' + "${bam_id}:${path}"
        }
    }
    output_file = "output.${threads}.${limit}.${bams.size()}.maf"
    """
    GetBaseCountsMultiSample \
    --omaf \
    --maq 20 \
    --baq 20 \
    --filter_improper_pair 0 \
    --output "${output_file}" \
    --thread "${threads}" \
    --maf "${maf}" \
    ${bam_arg} \
    --fasta "${ref_fasta}"
    """
}



workflow {
    log.info("----------------")
    log.info("workflow params:")
    log.info("${params}")
    log.info("----------------")


    ref_fasta = Channel.fromFilePairs("${params.fasta_base}{.dict,.fasta.amb,.fasta.ann,.fasta.bwt,.fasta.dict,.fasta.fai,.fasta.index,.fasta.pac,.fasta.sa}") | groupTuple | map { key, list ->
        return list
        } | flatten | collect | toList | set{ ref_fasta_indexes }

    ref_fasta = Channel.fromPath("${params.fasta_base}.fasta")
    maf = Channel.fromPath("${params.maf}")


    // Channel to one run sample bam
    bam_list = Channel.fromPath("${params.bam_list}")
    bam_list | splitCsv(sep:'\t') | map { row ->
        def bam = file(row[0])
        def bai = file(row[1])
        def sampleId = row[2]

        return [ bam, bai, sampleId]
    } | set { bam_bai_sample }
    // bam_bai_sample | combine(maf) | combine(ref_fasta) | combine(ref_fasta_indexes) | set { input_ch }
    // uncomment this section to run the single-sample process;
    // GBCMS_one_sample(input_ch, num_threads)


    // channels to run all bams
    Channel.fromPath("${params.bam_list}") | splitCsv(sep:'\t') | map { row ->
        def bam = file(row[0])
        def bai = file(row[1])
        def sampleId = row[2]

        return [ bam ]
    } | collect | toList | set { all_bams }

    Channel.fromPath("${params.bam_list}") | splitCsv(sep:'\t') | map { row ->
        def bam = file(row[0])
        def bai = file(row[1])
        def sampleId = row[2]

        return [ bai ]
    } | collect | toList | set { all_bais }

    all_bams | combine(all_bais) | combine(maf) | combine(ref_fasta) | combine(ref_fasta_indexes) | set { all_inputs }
    GBCMS_multi_sample(all_inputs, num_threads, limits, reps)
}
