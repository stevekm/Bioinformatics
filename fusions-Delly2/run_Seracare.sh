#!/bin/bash
# set -x

# This script is a template for running Delly2 for fusion detection on Seracare samples, submitting each job to the HPC cluster


# https://github.com/dellytools/delly
# The SV type can be DEL, DUP, INV, TRA, or INS for deletions, tandem duplications, inversions, translocations and small insertions, respectively.

# Usage: delly call [OPTIONS] -g <ref.fa> <sample1.sort.bam> <sample2.sort.bam> ...
#
# Generic options:
#   -? [ --help ]                     show help message
#   -t [ --type ] arg (=DEL)          SV type (DEL, DUP, INV, BND, INS)


# ~~~~~ SETUP ~~~~~ #
proj_dir="/ifs/data/molecpathlab/NGS580_WES-development/sns-wes-fusion-analysis"

output_dir="${proj_dir}/output"
log_dir="${proj_dir}/logs"

mkdir -p "$output_dir"
mkdir -p "$log_dir"

cd "$proj_dir"

# find the sample IDs for the Seracare samples
analysis_results_dir="Seracare_results" # its a symlink
bam_dir="${analysis_results_dir}/BAM-GATK-RA-RC"
samplesheet="${proj_dir}/Seracare_samplesheet.tsv"


make_Seracare_samplesheet () {
    echo "Making the Seracare sample sheet"
    samples_list="$(cat "${analysis_results_dir}/samples.fastq-raw.csv" | cut -d ',' -f1 | sort -u | xargs)"
    for sample in $samples_list; do
        bam_file="$(find "$bam_dir" -name "${sample}*" -name "*.bam" -exec readlink -f {} \;)"
        printf "%s\t%s\n" "$sample" "$bam_file" >> "$samplesheet"
    done
}

[ ! -f "$samplesheet" ] && make_Seracare_samplesheet


# ~~~~~ FUNCTIONS ~~~~~ #
timestamp () {
    printf "%s" "$(date +"%Y-%m-%d-%H-%M-%S")"
}

bcf2vcf () {
    local input_bcf="$1"
    local output_vcf="$2"
    local delly_dir="/ifs/home/kellys04/software/delly"
    local delly_bcftools="${delly_dir}/src/bcftools/bcftools"
    # ./delly/src/bcftools/bcftools view delly.bcf > delly.vcf
    # [ -f "${sample_ID}.bcf" ] && $delly_bcftools view "${sample_ID}.bcf" > "${sample_ID}.vcf"
    [ ! -f "$output_vcf" ] && $delly_bcftools view "$input_bcf" > "$output_vcf"
}

run_delly2 () {
    local sample_bam="$1"
    local sample_ID="$2"
    local call_types="deletions|DEL duplications|DUP inversions|INV translocations|BND insertions|INS"
    local delly_dir="/ifs/home/kellys04/software/delly"
    local delly_bin="${delly_dir}/src/delly"
    local hg19_fa="/local/data/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"


    # [ ! -f "${sample_ID}.bcf" ] && $delly_bin call -t BND -g "$hg19_fa" -o "${output_dir}/${sample_ID}.bcf" "$sample_bam"
    for call_type in $call_types; do
        local type_command="$(echo "$call_type" | cut -d '|' -f2)"
        local type_ID="$(echo "$call_type" | cut -d '|' -f1)"
        local sample_output_SV_bcf="${output_dir}/${sample_ID}.${type_ID}.bcf"
        local sample_output_SV_vcf="${sample_output_SV_bcf%%.bcf}.vcf"
        local sample_log_file="${log_dir}/${sample_ID}.${type_ID}.$(timestamp).log"
        local sample_output_regenotyped_bcf="${output_dir}/${sample_ID}.${type_ID}.geno.bcf"
        (
        qsub -wd $PWD -o :${log_dir}/ -e :${log_dir}/ -j y -N "$sample_ID" <<E0F
        set -x
        timestamp () {
            printf "%s" "\$(date +"%Y-%m-%d-%H-%M-%S")"
        }

        bcf2vcf () {
            local input_bcf="\$1"
            local output_vcf="\$2"
            local delly_dir="/ifs/home/kellys04/software/delly"
            local delly_bcftools="\${delly_dir}/src/bcftools/bcftools"
            [ ! -f "\$output_vcf" ] && \$delly_bcftools view "\$input_bcf" > "\$output_vcf"
        }

        # SV calling
        # delly call -t DEL -x hg19.excl -o t1.bcf -g hg19.fa tumor1.bam control1.bam
        printf "%s\tStarting Delly2 SV calling\n" "\$(timestamp)"
        [ ! -f "${sample_output_SV_bcf}" ] && ${delly_bin} call -t ${type_command} -g "${hg19_fa}" -o "${sample_output_SV_bcf}" "${sample_bam}"

        printf "%s\tFinished Delly2 SV calling\n" "\$(timestamp)"
        # convert BCF to VCF
        printf "%s\tConverting BCF to VCF Delly2 SV calling\n" "\$(timestamp)"

        bcf2vcf "$sample_output_SV_bcf" "$sample_output_SV_vcf"
        printf "%s\tFinished Converting BCF to VCF Delly2 SV calling\n" "\$(timestamp)"

        # Re-genotype merged SV site list
        if [ "$type_command" == "BND" ]; then

            printf "%s\tStarting Delly2 Re-genotyping calling\n" "\$(timestamp)"
            [ ! -f "${sample_output_regenotyped_bcf}" ] && $delly_bin call -t BND -g "$hg19_fa" -v "${sample_output_SV_bcf}" -o "${sample_output_regenotyped_bcf}" "$sample_bam"
            printf "%s\tFinished Delly2 Re-genotyping calling\n" "\$(timestamp)"
        fi
E0F
sleep 1
        )
    done
}



# ~~~~~ RUN ~~~~~ #
cat "$samplesheet" | while read line; do
    if [ ! -z "$line" ]; then
        sample_ID="$(echo "$line" | cut -f1)"
        sample_bam="$(echo "$line" | cut -f2)"
        run_delly2 "$sample_bam" "$sample_ID"
    fi
done
