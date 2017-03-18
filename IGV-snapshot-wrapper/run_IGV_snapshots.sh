#/!bin/bash


proj_dir="/ifs/home/kellys04/projects/SmithLab_ChIpSeq_2017-02-11/project_notes/IGV-snapshots"
bam_outdir="${proj_dir}/bam_links"
snapshot_outdir="${proj_dir}/snapshots"
IGV_script="${proj_dir}/make_IGV_snapshots.py"
IGV_script="$(readlink -f "$IGV_script")"

regions_file="/ifs/home/kellys04/projects/SmithLab_ChIpSeq_2017-02-11/project_notes/IGV-snapshots/shella_regions.bed"

find "$bam_outdir" -mindepth 2 -type d | while read item; do
    sampledir="$item"
    # echo "$sampledir"
    sampleID="$(basename "$sampledir")"
    sample_group="$(basename $(dirname "$sampledir"))"
    printf "%s\t%s\n" "$sample_group" "$sampleID"
    sample_snapshot_outdir="${snapshot_outdir}/${sample_group}/${sampleID}"
    # echo "$sample_snapshot_outdir"
    mkdir -p "$sample_snapshot_outdir"
    bam_files="$(find "$sampledir" -name "*.bam")"
    # echo $bam_files
    # set -x
    python $IGV_script -r "$regions_file" -o "$sample_snapshot_outdir" $bam_files
    # set +x
    echo ""
done

