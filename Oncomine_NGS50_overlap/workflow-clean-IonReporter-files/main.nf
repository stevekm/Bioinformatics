params.filtered_tsvs_list = "input_tsv_list.txt"
params.inputDir = "input"
params.output_dir = "output"

// params.inputDir = "IonReporter_files"
//
// Channel.fromPath( file(params.filtered_tsvs_list) )
//     .splitCsv() // read lines from file
Channel.fromPath("${params.inputDir}/**-oncomine.tsv")
    .map { item ->
        // parse header comment
        def lines = new File("${item}").readLines().grep(~/^#.*/)
        lines = lines.collect { "${it}".replaceFirst(/^##/, "") }
        def fileparams = lines.collect { "${it}".tokenize('=') }

        return [ item, fileparams ]
    }.set { all_tsvs }
    // .subscribe { println "${it}" }



process add_params {
    // publishDir "${params.output_dir}/cols", overwrite: true, mode: 'copy'
    executor "sge"

    input:
    set file(filtered_tsv), val(fileparams) from all_tsvs

    output:
    file("${output_file}") into updated_tsvs

    script:
    // output_file = "${filtered_tsv}".replaceFirst(/.tsv$/, ".headercols.tsv")
    output_file = "output.tsv"
    param_string = fileparams.collect { "paste-col.py --header '${it[0]}' -v '${it[1]}' --doublequote" }.join(" | ")
    """
    cat "${filtered_tsv}" | \
    grep -v '^#' | \
    ${param_string} | \
    paste-col.py --header 'Source' -v '${filtered_tsv}' --doublequote \
    > "${output_file}"
    """
}

process collect_tables {
    publishDir "${params.output_dir}", overwrite: true, mode: 'copy'

    input:
    file('t*') from updated_tsvs.collect()

    output:
    file("${output_file}")
    file("${output_zip}")

    script:
    output_file = 'Oncomine_all_filtered_annotations.tsv'
    output_zip = 'Oncomine_all_filtered_annotations.zip'
    """
    concat-tables.py * > "${output_file}"
    zip "${output_zip}" "${output_file}"
    """
}
