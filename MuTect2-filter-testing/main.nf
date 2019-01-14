params.outputDir = "output"
def outputDir = "${params.outputDir}"
params.inputDir = "input"
params.ref_dir = "/gpfs/scratch/kellys04/molecpathlab/ref"
params.ref_fa = "${params.ref_dir}/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
params.ref_fa_bwa_dir = "${params.ref_dir}/BWA/hg19"
params.ref_fai = "${params.ref_dir}/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa.fai"
params.ref_dict = "${params.ref_dir}/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.dict"
params.ref_chrom_sizes = "${params.ref_dir}/Illumina/hg19/chrom.sizes"
params.microsatellites = "${params.ref_dir}/msisensor/hg19/microsatellites.list"
params.trimmomatic_contaminant_fa = "${params.ref_dir}/contaminants/trimmomatic.fa"
params.gatk_bundle_dir = "${params.ref_dir}/gatk-bundle"
params.gatk_1000G_phase1_indels_hg19_vcf = "${params.gatk_bundle_dir}/1000G_phase1.indels.hg19.vcf"
params.mills_and_1000G_gold_standard_indels_hg19_vcf = "${params.gatk_bundle_dir}/Mills_and_1000G_gold_standard.indels.hg19.vcf"
params.dbsnp_ref_vcf = "${params.gatk_bundle_dir}/dbsnp_138.hg19.vcf"
params.cosmic_ref_vcf = "${params.ref_dir}/hg19/CosmicCodingMuts_v73.hg19.vcf"

Channel.fromPath( file(params.ref_fa) ).into { ref_fasta; ref_fasta2 }
Channel.fromPath( file(params.ref_fai) ).into { ref_fai; ref_fai2 }
Channel.fromPath( file(params.ref_dict) ).into { ref_dict; ref_dict2 }

Channel.fromPath("variants-NGS580-nf/*.MuTect2.norm.vcf").set { norm_vcfs }
// .subscribe { println "${it}" }

process concat_vcf {
    publishDir "${outputDir}", overwrite: true, mode: 'copy'

    input:
    file(all_vcf: "*") from norm_vcfs.collect()
    set file(ref_fasta), file(ref_fai), file(ref_dict) from ref_fasta2.combine(ref_fai2).combine(ref_dict2)

    output:
    file("${concat_vcf}") into concatenated_vcf
    file("${concat_vcf_idx}") into concatenated_vcf_idx
    file("${concat_tsv}")

    script:
    concat_vcf = "concat_variants.vcf"
    concat_vcf_idx = "${concat_vcf}.idx"
    concat_tsv = "concat_variants.tsv"
    fixed_tsv = "concat_variants.fixed.tsv"
    """
    # print first header
    grep '^#' \$(echo "${all_vcf}" |  tr ' ' '\n' | head -1) > "${concat_vcf}"

    # print all file contents minus headers
    cat ${all_vcf} | grep -v '^#' | sort -V -k1,1 -k2,2n >> "${concat_vcf}"

    # convert VCF to TSV
    gatk.sh -T VariantsToTable \
    -R "${ref_fasta}" \
    -V "${concat_vcf}" \
    -F CHROM -F POS -F ID -F REF -F ALT -F FILTER -F QUAL -F AC -F AN -F NLOD -F TLOD \
    -GF AD -GF DP -GF AF \
    -o "${concat_tsv}"

    # fix the AD columns
    split.R "${concat_tsv}" "${fixed_tsv}"
    """
}

// messed up:
// tumor_AD_5x_normal_AD
// tumor_AD_gt_3pcnt
def args = [
["tumor_AD_gt_3pcnt", """ -select "(vc.getGenotype('TUMOR').getAD().1 / (vc.getGenotype('TUMOR').getAD().0 + vc.getGenotype('TUMOR').getAD().1) )  > 0.03" """],
["normal_AD_lt_5pcnt", """ -select "(vc.getGenotype('NORMAL').getAD().1 / (vc.getGenotype('NORMAL').getAD().0 + vc.getGenotype('NORMAL').getAD().1) )  < 0.05" """],
["tumor_AD_gt_5", """ -select "vc.getGenotype('TUMOR').getAD().1 > 5"  """],
["tumor_AD_5x_normal_AD", """-select "(vc.getGenotype('TUMOR').getAD().1 / (vc.getGenotype('TUMOR').getAD().0 + vc.getGenotype('TUMOR').getAD().1)) > ((vc.getGenotype('NORMAL').getAD().1 / (vc.getGenotype('NORMAL').getAD().0 + vc.getGenotype('NORMAL').getAD().1)) * 5)" """]
]
comb = []
1.upto(args.size()) {
    [args].multiply(it).eachCombination { list ->
      if(list.size() == 1 || (1..<list.size()).every { list[it - 1][0] < list[it][0] }) {
           comb << list
      }
    }
}

Channel.from(comb)
.set { comb_params }
// concatenated_vcf.combine(comb_params)
// .subscribe { println "$it" }
process filter_variants {
    tag "${label}"
    publishDir "${outputDir}", overwrite: true, mode: 'copy'

    input:
    set val(params), file(vcf), file(vcf_idx), file(ref_fasta), file(ref_fai), file(ref_dict) from comb_params.map{[it]}.combine(concatenated_vcf).combine(concatenated_vcf_idx).combine(ref_fasta).combine(ref_fai).combine(ref_dict)

    output:
    set file("${output_vcf}"), file("${output_tsv}"), val("${label}") into filtered_vcf

    script:
    label = params.collect { it[0] }.join('.')
    arg_list = params.collect { it[1] }.join(' ')
    output_vcf = "${label}.vcf"
    output_tsv = "${label}.tsv"
    """
    gatk.sh -T SelectVariants \
    -R "${ref_fasta}" \
    -V "${vcf}" \
    ${arg_list} > "${output_vcf}"

    # convert VCF to TSV
    gatk.sh -T VariantsToTable \
    -R "${ref_fasta}" \
    -V "${output_vcf}" \
    -F CHROM -F POS -F ID -F REF -F ALT -F FILTER -F QUAL -F AC -F AN -F NLOD -F TLOD \
    -GF AD -GF DP -GF AF \
    -o "${output_tsv}"
    """
}

process count_variants {
    tag "${label}"
    publishDir "${outputDir}", overwrite: true, mode: 'copy'

    input:
    set file(vcf), file(tsv), val(label) from filtered_vcf

    output:
    file("${count_vcf_file}")

    script:
    count_vcf_file = "${label}.counts.txt"
    """
    num_variants="\$(grep -v '^#' "${vcf}" | wc -l )"
    printf "\${num_variants}\t${label}\n" > "${count_vcf_file}"
    """
}
//
// # filter VCF
// # report if:
// # T frequency is more than 3%
// # N frequency is less than 5%
// # at least 5 variant call supporting reads
// # T frequency is sufficiently higher (5x) than N frequency
// # "we recommend applying post-processing filters, e.g. by hard-filtering calls with low minor allele frequencies"
// # only 'PASS' entries
// gatk.sh -T SelectVariants \
// -R "${ref_fasta}" \
// -V "${norm_vcf}" \
// -select "(vc.getGenotype('TUMOR').getAD().1 / (vc.getGenotype('TUMOR').getAD().0 + vc.getGenotype('TUMOR').getAD().1) )  > 0.03" \
// -select "(vc.getGenotype('NORMAL').getAD().1 / (vc.getGenotype('NORMAL').getAD().0 + vc.getGenotype('NORMAL').getAD().1) )  < 0.05" \
// -select "vc.getGenotype('TUMOR').getAD().1 > 5" \
// -select "(vc.getGenotype('TUMOR').getAD().1 / (vc.getGenotype('TUMOR').getAD().0 + vc.getGenotype('TUMOR').getAD().1) ) > (vc.getGenotype('NORMAL').getAD().1 / (vc.getGenotype('NORMAL').getAD().0 + vc.getGenotype('NORMAL').getAD().1) ) * 5" \
// -select 'vc.isNotFiltered()' \
// > "${filtered_vcf}"
//
// # vc.getGenotype('TUMOR').getAD().0 ; Tumor Allelic depth ref
// # vc.getGenotype('TUMOR').getAD().1 ; Tumor Allelic depth alt
// # vc.getGenotype('NORMAL').getAD().0 ; Normal Allelic depth ref
// # vc.getGenotype('NORMAL').getAD().1 ; Normal Allelic depth alt
//
// # ( Tumor Allelic depth alt / (Tumor Allelic depth ref + Tumor Allelic depth alt ) )  > 0.03
// # ( Normal Allelic depth alt / ( Normal Allelic depth ref + Normal Allelic depth alt ) )  < 0.05
// # Tumor Allelic depth alt > 5
// # ( Tumor Allelic depth alt / ( Tumor Allelic depth ref + Tumor Allelic depth alt ) ) > ( Normal Allelic depth alt / ( Normal Allelic depth ref + Normal Allelic depth alt ) ) * 5
//
// ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
// # CHROM	chr7
// # POS	2946342
// # ID	.
// # REF	A
// # ALT	G
// # QUAL	.
// # FILTER	clustered_events;homologous_mapping_event
// # INFO	ECNT=16;HCNT=3;MAX_ED=65;MIN_ED=5;NLOD=33.41;TLOD=23.04
// # FORMAT	GT:AD:AF:ALT_F1R2:ALT_F2R1:FOXOG:PGT:PID:QSS:REF_F1R2:REF_F2R1
// # TUMOR	0/1:1333,17:0.013:7:10:0.588:0|1:2946342_A_G:40125,535:641:689
// # NORMAL	0/0:137,0:0.00:0:0:.:0|1:2946342_A_G:3959,0:53:80


// 380115  normal_AD_lt_5pcnt
// 1       normal_AD_lt_5pcnt.tumor_AD_5x_normal_AD
// 1       normal_AD_lt_5pcnt.tumor_AD_5x_normal_AD.tumor_AD_gt_3pcnt
// 0       normal_AD_lt_5pcnt.tumor_AD_5x_normal_AD.tumor_AD_gt_3pcnt.tumor_AD_gt_5
// 0       normal_AD_lt_5pcnt.tumor_AD_5x_normal_AD.tumor_AD_gt_5
// 1       normal_AD_lt_5pcnt.tumor_AD_gt_3pcnt
// 0       normal_AD_lt_5pcnt.tumor_AD_gt_3pcnt.tumor_AD_gt_5
// 251089  normal_AD_lt_5pcnt.tumor_AD_gt_5
// 1       tumor_AD_5x_normal_AD
// 1       tumor_AD_5x_normal_AD.tumor_AD_gt_3pcnt
// 0       tumor_AD_5x_normal_AD.tumor_AD_gt_3pcnt.tumor_AD_gt_5
// 0       tumor_AD_5x_normal_AD.tumor_AD_gt_5
// 1       tumor_AD_gt_3pcnt
// 0       tumor_AD_gt_3pcnt.tumor_AD_gt_5
// 251089  tumor_AD_gt_5
// //
//
// ok:
// tumor_AD_gt_5
// normal_AD_lt_5pcnt
//
// messed up:
// tumor_AD_5x_normal_AD
// tumor_AD_gt_3pcnt
