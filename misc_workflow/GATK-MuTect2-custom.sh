#!/bin/bash

# run the custom variant calling using known positive and negative control samples with GATK MuTect2

timestamp () {
    printf "%s" "\$(date +"%Y-%m-%d-%H-%M-%S")"
}

input_files="
result1|/ifs/data/molecpathlab/NGS580_WES/170721_NB501073_0018_AH5C7GBGX3/results_2017-07-28_21-02-00/BAM-GATK-RA-RC/HapMap-250ng-1.dd.ra.rc.bam|/ifs/data/molecpathlab/NGS580_WES/170721_NB501073_0018_AH5C7GBGX3/results_2017-07-28_21-02-00/BAM-GATK-RA-RC/SeraCare-1to1-250ng-1.dd.ra.rc.bam
result2|/ifs/data/molecpathlab/NGS580_WES/170721_NB501073_0018_AH5C7GBGX3/results_2017-07-28_21-02-00/BAM-GATK-RA-RC/HapMap-250ng-2.dd.ra.rc.bam|/ifs/data/molecpathlab/NGS580_WES/170721_NB501073_0018_AH5C7GBGX3/results_2017-07-28_21-02-00/BAM-GATK-RA-RC/SeraCare-1to1-250ng-2.dd.ra.rc.bam
result3|/ifs/data/molecpathlab/NGS580_WES/170721_NB501073_0018_AH5C7GBGX3/results_2017-07-28_21-02-00/BAM-GATK-RA-RC/HapMap-250ng-3.dd.ra.rc.bam|/ifs/data/molecpathlab/NGS580_WES/170721_NB501073_0018_AH5C7GBGX3/results_2017-07-28_21-02-00/BAM-GATK-RA-RC/SeraCare-1to1-250ng-3.dd.ra.rc.bam
"
file_ext=".dd.ra.rc.bam"

for item in $input_files; do
    result_id="$(echo "$item" | cut -d '|' -f1)"

    normal_bam="$(echo "$item" | cut -d '|' -f2)"
    tumor_bam="$(echo "$item" | cut -d '|' -f3)"

    normal_name="$(basename "$normal_bam")"
    normal_name="${normal_name%%$file_ext}"

    tumor_name="$(basename "$tumor_bam")"
    tumor_name="${tumor_name%%$file_ext}"

    output_name="${tumor_name}-${normal_name}"
    output_dir="$(readlink -f $PWD)/${output_name}"
    output_vcf="${output_dir}/${output_name}.vcf"
    mkdir -p "$output_dir"

    mutect_cmd="$(
    printf "
java -Xms16G -Xmx16G -jar /ifs/home/id460/software/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar -T MuTect2 \
-dt NONE \
--logging_level WARN \
--standard_min_confidence_threshold_for_calling 30 \
--max_alt_alleles_in_normal_count 10 \
--max_alt_allele_in_normal_fraction 0.05 \
--max_alt_alleles_in_normal_qscore_sum 40 \
--reference_sequence /ifs/data/sequence/Illumina/igor/ref/hg19/genome.fa \
--dbsnp /ifs/home/id460/ref/hg19/gatk-bundle/dbsnp_138.hg19.vcf \
--cosmic /ifs/home/id460/ref/hg19/CosmicCodingMuts_v73.hg19.vcf \
--intervals /ifs/data/molecpathlab/NGS580_WES/170602_NB501073_0012_AHCKYCBGX2/results_2017-06-05_09-51-00/NGS580_targets.bed \
--interval_padding 10 \
--input_file:tumor %s \
--input_file:normal %s \
--out %s" "$tumor_bam" "$normal_bam" "$output_vcf"
    )"

    echo "$result_id"
    echo "$normal_name $normal_bam"
    echo "$tumor_name $tumor_bam"

    echo "$output_vcf"
    echo "$mutect_cmd"
    echo ""


    qsub -wd ${output_dir} -j y -N "${output_name}" -o :${output_dir}/ -e :${output_dir}/ <<E0F
set -x
printf "qsub job started at: %s" "$(timestamp)"
$mutect_cmd
E0F

done
