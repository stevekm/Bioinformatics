params.tcgaDir = "tcga"
def tcgaDir = "${params.tcgaDir}"
params.outputDir = "output"
def outputDir = "${params.outputDir}"
params.inputDir = "gdc"

Channel.fromPath("${params.inputDir}/**.maf").map{ item ->
    def basename = new File("${item}").getName()
    def parts = basename.tokenize(".")
    def study = parts[0]
    def project = parts[1]
    def caller = parts[2]
    return [item, study, project, caller]
}.set { input_mafs }
// .subscribe{ println "${it}" }

process filter_variants {
    tag "${prefix}"
    publishDir "${outputDir}/${prefix}", overwrite: true, mode: 'copy'
    label "R"
    echo true

    input:
    set file(maf), val(study), val(project), val(caller) from input_mafs

    output:
    set file("${output_maf}"), val(study), val(project), val(caller) into filtered_mafs

    script:
    prefix = "${study}.${project}.${caller}"
    output_maf = "${maf}".replaceFirst(/.maf$/, ".filtered.maf")
    """
    filter-variants.R "${maf}" "${output_maf}"
    """
}

process count_variants {
    tag "${prefix}"
    publishDir "${outputDir}/${prefix}", overwrite: true, mode: 'copy'
    label "R"

    input:
    set file(maf), val(study), val(project), val(caller) from filtered_mafs

    output:
    set file("${output_counts}"), val(study), val(project), val(caller) into (variant_counts, variant_counts2)
    set file("${output_aachange}"), val(study), val(project), val(caller) into aa_counts

    script:
    prefix = "${study}.${project}.${caller}"
    output_counts = "${prefix}.counts.tsv"
    output_aachange = "${prefix}.aachange.tsv"
    """
    variant-counts.R "${maf}" "${output_counts}" "${output_aachange}" "${study}" "${project}" "${caller}"
    """
}

// get all the tumor bam file uuid's from the table
variant_counts.map{ tsv, study, project, caller ->
    def tumor_bam_uuids = []
    tsv.withReader { reader ->
        while (line = reader.readLine()) {
            def tumor_bam_uuid = line.tokenize("\t")[1]
            tumor_bam_uuids.add(tumor_bam_uuid)
        }
    }
    return tumor_bam_uuids
}.collect().flatten().unique().set { all_tumor_bam_uuids }

def tcga_bed_urls = false
process get_tcga_bed_urls {
    tag "${prefix}"
    maxRetries 3
    errorStrategy 'retry'
    // maxForks 25
    storeDir "${tcgaDir}" // dont repeat this if its already done
    label "curl"

    input:
    val(bam_uuid) from all_tumor_bam_uuids

    output:
    set val(bam_uuid), file("${output_json}") into tumor_bam_bed_urls

    when:
    tcga_bed_urls == true

    script:
    prefix = "${bam_uuid}"
    output_json = "${prefix}.target_region_urls.json"
    """
    sleep "\$((RANDOM % 30 * ${task.attempt}))"
    get-gdc-bam-capture-kit.py "${bam_uuid}" "${output_json}"
    """
    // get-gdc-bam-capture-kit.py "492b84ee-c167-4a63-8dc7-0c5b43530094" "492b84ee-c167-4a63-8dc7-0c5b43530094.target_region_urls.json"
}


process calculate_tmb {
    tag "${prefix}"
    publishDir "${outputDir}/${prefix}", overwrite: true, mode: 'copy'

    input:
    set file(agg_tsv), val(study), val(project), val(caller) from variant_counts2

    output:
    set file("${output_tsv}"), val(study), val(project), val(caller) into tmb_tables

    script:
    prefix = "${study}.${project}.${caller}"
    output_tsv = "${prefix}.counts.tmb.tsv"
    // https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-017-0424-2
    // We used 38 Mb as the estimate of the exome size.
    """
    divide-col.py -i "${agg_tsv}" -f 3 --divisor 38 > "${output_tsv}"
    """
}
// tmb_tables.collectFile(name: "somefile.tsv", storeDir: "${outputDir}", keepHeader: true)

// get all the aggregate tables
tmb_tables.map{ tsv, study, project, caller ->
    return tsv
}.set { all_tmb_tables }

process collect_aggregate_tables {
    publishDir "${outputDir}", overwrite: true, mode: 'copy'

    input:
    file('t*') from all_tmb_tables.collect()

    output:
    file("${output_tsv}") into (all_agg_tsv, all_agg_tsv2)

    script:
    output_tsv = "all_variant_counts.tsv"
    """
    cat t* > "${output_tsv}"
    """
}

// get all the aa change aggregate tables
aa_counts.map {tsv, study, project, caller ->
    return tsv
}.set { all_aa_tables }
process collect_aachange_agg_tables {
    publishDir "${outputDir}", overwrite: true, mode: 'copy'

    input:
    file('t*') from all_aa_tables.collect()

    output:
    file("${output_tsv}") into (all_aachange_tsv, all_aachange_tsv2)

    script:
    output_tsv = "all_aa_counts.tsv"
    """
    cat t* > "${output_tsv}"
    """
}

process plot_results {
    publishDir "${outputDir}", overwrite: true, mode: 'copy'
    label "plot"

    input:
    set file(aa_tsv), file(agg_tsv) from all_aachange_tsv.combine(all_agg_tsv)

    output:
    file('*.pdf') into output_plots

    script:
    output_aa_name = "aa_changes"
    output_agg_counts_name = "variant_counts"
    """
    agg-plot.R "${aa_tsv}" "${agg_tsv}" "${output_aa_name}" "${output_agg_counts_name}"
    """
}

process gather_results {
    publishDir "${outputDir}", overwrite: true, mode: 'copy'
    stageInMode "copy"

    input:
    file('*') from output_plots.concat(all_aachange_tsv2, all_agg_tsv2).collect()

    output:
    file("${output_zip}")

    script:
    output_zip = "results.zip"
    """
    zip "${output_zip}" *
    """
}
