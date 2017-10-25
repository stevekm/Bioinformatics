#!/bin/bash

# run MuTect2 on intervals per chromosome in the reference genome

# get chrom coordinates
# reference_sequence="/ifs/data/sequence/Illumina/igor/ref/hg19/genome.fa"
# faidx "$reference_sequence" -i chromsizes > sizes.genome

# chrM	16571
# chr1	249250621
# chr2	243199373
# chr3	198022430
# chr4	191154276
# chr5	180915260
# chr6	171115067
# chr7	159138663
# chr8	146364022
# chr9	141213431
# chr10	135534747
# chr11	135006516
# chr12	133851895
# chr13	115169878
# chr14	107349540
# chr15	102531392
# chr16	90354753
# chr17	81195210
# chr18	78077248
# chr19	59128983
# chr20	63025520
# chr21	48129895
# chr22	51304566
# chrX	155270560
# chrY	59373566

cd "/ifs/home/kellys04/projects/Clinical_580_gene_panel/molecpathlab/NGS580_WES-development/MuTect2-split"

GATK_jar='/ifs/home/id460/software/GenomeAnalysisTK/GenomeAnalysisTK-3.8-0/GenomeAnalysisTK.jar'
qsub_logdir="/ifs/data/molecpathlab/NGS580_WES-development/MuTect2-split/logs"
reference_sequence="/ifs/data/sequence/Illumina/igor/ref/hg19/genome.fa"
intervals='targets.bed'
tumor_bam='Sample1.dd.ra.rc.bam'
normal_bam='HapMap-B17-1267.dd.ra.rc.bam'
dbsnp='/ifs/home/id460/ref/hg19/gatk-bundle/dbsnp_138.hg19.vcf'
cosmic='/ifs/home/id460/ref/hg19/CosmicCodingMuts_v73.hg19.vcf'



run_MuTect2 () {
    local chrom="$1"
    local end="$2"
    local interval="${chrom}:1-${end}"
    local job_name="MuTect2_${chrom}"
    printf 'interval is %s' "$interval"

    qsub -wd $PWD -o :${qsub_logdir}/ -e :${qsub_logdir}/ -j y -N "$job_name" <<E0F
set -x
date +"%Y-%m-%d-%H-%M-%S"
start=\$(date +%s)

java -Xms16G -Xmx16G -jar $GATK_jar -T MuTect2 \
-dt NONE \
--logging_level WARN \
--standard_min_confidence_threshold_for_calling 30 \
--max_alt_alleles_in_normal_count 10 \
--max_alt_allele_in_normal_fraction 0.05 \
--max_alt_alleles_in_normal_qscore_sum 40 \
--reference_sequence $reference_sequence \
--dbsnp $dbsnp \
--cosmic $cosmic \
--intervals $intervals \
--interval_padding 10 \
--input_file:tumor $tumor_bam \
--input_file:normal $normal_bam \
--out ${chrom}.vcf \
-L $interval

date +"%Y-%m-%d-%H-%M-%S"
echo "Duration: \$(((\$(date +%s)-\$start)/60)) minutes"
set +x
E0F

}


cat sizes.genome | while read line; do
    if [ ! -z "$line" ]; then
        # echo "$line"
        chrom="$(echo "$line" | cut -f1)"
        end="$(echo "$line" | cut -f2)"
        run_MuTect2 "$chrom" "$end"
    fi
    echo ''
done
