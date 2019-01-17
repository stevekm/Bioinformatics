params.runID = "170602_NB501073_0012_AHCKYCBGX2"
params.runDir = "data/${params.runID}"
params.inputDir = "${params.runDir}/variants/"
params.outputDir = "${params.runDir}/output/"

params.ref_dir = "/gpfs/scratch/kellys04/molecpathlab/ref"
params.ref_fa = "${params.ref_dir}/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
params.ref_fai = "${params.ref_dir}/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa.fai"
params.ref_dict = "${params.ref_dir}/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.dict"
params.ANNOVAR_DB_DIR = "${params.ref_dir}/annovar/db"

Channel.fromPath( file(params.ref_fa) ).into { ref_fasta; ref_fasta2 }
Channel.fromPath( file(params.ref_fai) ).into { ref_fai; ref_fai2 }
Channel.fromPath( file(params.ref_dict) ).into { ref_dict; ref_dict2 }
Channel.fromPath("${params.ANNOVAR_DB_DIR}").set { annovar_db_dir }

Channel.fromPath("${params.inputDir}/*.MuTect2.norm.vcf").map { item -> // Sample1-Tumor_Sample1-Normal.chr14.MuTect2.norm.vcf
        def itemID = item.baseName // Sample1-Tumor_Sample1-Normal.chr14.MuTect2.norm
        def sampleID = "${itemID}".replaceFirst(/.chr[0-9XYM]*.MuTect2.norm$/, "") // Sample1-Tumor_Sample1
        return([ sampleID, itemID, item ])
    }
    .set { input_vcfs }
// input_vcfs.subscribe { println "${it}" }

process filter_vcf {
    // filter the .vcf to only include 'PASS' entries
    tag "${itemID}"
    publishDir "${params.outputDir}", mode: 'copy', overwrite: true

    input:
    set val(sampleID), val(itemID), file(vcf) from input_vcfs

    output:
    set val(sampleID), val(itemID), file("${output_vcf}") into filtered_vcfs

    script:
    output_vcf = "${itemID}.filter.vcf"
    """
    # get the header
    grep '^#' "${vcf}" > "${output_vcf}"

    # get the 'PASS' entries
    grep -v '^#' "${vcf}" | grep 'PASS' >> "${output_vcf}" || :
    """
}

process vcf_to_tsv {
    // convert the .vcf to .tsv format
    tag "${itemID}"
    // publishDir "${params.outputDir}", overwrite: true, mode: 'copy'

    input:
    set val(sampleID), val(itemID), file(vcf), file(ref_fasta), file(ref_fai), file(ref_dict) from filtered_vcfs.combine(ref_fasta).combine(ref_fai).combine(ref_dict)

    output:
    set val(sampleID), val(itemID), file(vcf), file("${output_tsv}") into vcfs_tsvs

    script:
    output_tsv = "${itemID}.tsv"
    """
    # NOTE: automatically filters for only PASS entries
    gatk.sh -T VariantsToTable \
    -R "${ref_fasta}" \
    -V "${vcf}" \
    -F CHROM -F POS -F ID -F REF -F ALT -F FILTER -F QUAL -F AC -F AN -F NLOD -F TLOD \
    -GF AD -GF DP -GF AF \
    -o "${output_tsv}"
    """
}


process reformat_vcf_tsv {
    // reformat and adjust the TSV table for consistency downstream
    tag "${itemID}"
    publishDir "${params.outputDir}", mode: 'copy', overwrite: true

    input:
    set val(sampleID), val(itemID), file(vcf), file(tsv) from vcfs_tsvs

    output:
    set val(caller), val(sampleID), val(itemID), file(vcf), file("${output_tsv}") into vcfs_tsvs_reformat

    script:
    caller = 'MuTect2'
    output_tsv = "${itemID}.reformat.tsv"
    """
    reformat-vcf-table.py -c MuTect2 -s "${sampleID}" -i "${tsv}" | \
    paste-col.py --header "Sample" -v "${sampleID}"  | \
    paste-col.py --header "Run" -v "${params.runID}" | \
    paste-col.py --header "VariantCaller" -v "${caller}" > "${output_tsv}"
    """
}

// make sure there are variants in the TSV
import java.nio.file.Files;
vcfs_tsvs_reformat.filter { caller, sampleID, itemID, vcf, tsv ->
        long count = Files.lines(tsv).count()
        if (count <= 1) println ">>> WARNING: file ${tsv} does not have enough lines and will not be included"
        count > 1
                    }
                    .set { vcfs_tsvs_reformat_filtered }



process annotate {
    // annotate the VCF file
    tag "${itemID}"
    publishDir "${params.outputDir}", mode: 'copy', overwrite: true

    input:
    set val(caller), val(sampleID), val(itemID), file(vcf), file(tsv), file(annovar_db_dir) from vcfs_tsvs_reformat_filtered.combine(annovar_db_dir)

    output:
    set val(caller), val(sampleID), val(itemID), file(vcf), file(tsv), file("${annovar_output_txt}"), file("${avinput_tsv}") into vcfs_tsvs_annotations

    script:
    prefix = "${itemID}"
    avinput_file = "${prefix}.avinput"
    avinput_tsv = "${avinput_file}.tsv"
    annovar_output_txt = "${prefix}.${params.ANNOVAR_BUILD_VERSION}_multianno.txt"
    // annovar_output_vcf = "${prefix}.${params.ANNOVAR_BUILD_VERSION}_multianno.vcf"
    """
    table_annovar.pl "${vcf}" "${annovar_db_dir}" \
    --buildver "${params.ANNOVAR_BUILD_VERSION}" \
    --remove \
    --protocol "${params.ANNOVAR_PROTOCOL}" \
    --operation "${params.ANNOVAR_OPERATION}" \
    --nastring . \
    --vcfinput \
    --onetranscript \
    --outfile "${prefix}"

    printf "Chr\tStart\tEnd\tRef\tAlt\tCHROM\tPOS\tID\tREF\tALT\n" > "${prefix}.avinput.tsv"
    cut -f1-5,9-13 ${avinput_file} >>  "${avinput_tsv}"
    """
}


process merge_tables {
    // merge the annotation and vcf tables
    tag "${itemID}"
    publishDir "${params.outputDir}", mode: 'copy', overwrite: true

    input:
    set val(caller), val(sampleID), val(itemID), file(vcf), file(tsv), file(annovar_txt), file(avinput_tsv) from vcfs_tsvs_annotations

    output:
    file("${output_annotations}") into merged_tables2

    script:
    prefix = "${itemID}"
    output_annotations = "${prefix}.annotations.tsv"
    """
    merge-vcf-tables.R "${tsv}" "${annovar_txt}" "${avinput_tsv}" "${output_annotations}"
    """

}

// concatenate all the variant tables
merged_tables2.collectFile(name: "annotations.${params.runID}.tsv", storeDir: "${params.outputDir}", keepHeader: true, sort: false)
