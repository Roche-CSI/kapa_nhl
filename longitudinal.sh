#!/usr/bin/env bash

# This is the Longitudinal module for the KAPA NHL Analysis Pipeline. 
# Uses the main module results of baseline, germline, and followup samples.
# Includes blocklist generation to longitudinal mutation analysis


# Generate background panel and blocklist 
# (performed only once for the analysis workflow setup)
# --------------------------------------
Rscript create_bg_panel.R \
    --bam_list umi_deduped_sorted_bams.csv \
    --panel_background panel_background.RDS \
    --target_bed DESIGN_capture_targets.bed \
    --blocklist_type variant \
    --reference BSgenome.Hsapiens.UCSC.hg38

Rscript create_blocklist.R \
    --panel_background panel_background.RDS \
    --blocklist blocklist.txt \
    --vaf_quantile 0.95 \
    --min_samples_one_read 12 \
    --min_samples_two_reads 10


# Select Baseline Reporter Variant Candidates
# --------------------------------------
Rscript select_reporters.R \
    --filter_reporters TRUE \
    --remove_snp TRUE \
    --reporters BASELINE_SAMPLE_vardict_annotated.vcf.txt \
    --germline GERMLINE_SAMPLE_umi_deduped.clipov.sorted.bam \
    --followup FOLLOWUP_SAMPLE_umi_deduped.clipov.sorted.bam \
    --retained selected_baseline_reporters.txt \
    --all all_reporters.txt \
    --blocklist blocklist.txt \
    --germline_cutoff 0.0005 \
    --min_af 0.005 \
    --max_af 0.35 \
    --min_dp 1000 \
    --min_vd 15 \
    --min_mq 55 \
    --min_qual 45 \
    --min_sbf 0.00001 \
    --max_nm 3 \
    --read_min_bq 30 \
    --read_min_mq 30 \
    --read_max_dp 20000


# Longitudinal Mutation Analysis
# --------------------------------------
Rscript longitudinal_analysis.R
    --reporters selected_baseline_reporters.txt \
    --sample_bam FOLLOWUP_SAMPLE_umi_deduped_clipov_sorted.bam \
    --target_bed target.bed \
    --blocklist blocklist.txt \
    --blocklist_type variant \
    --reads_threshold 1000 \
    --pvalue_threshold 0.001 \
    --vaf_threshold 0.1 \
    --read_min_bq 30 \
    --read_min_mq 40 \
    --read_max_dp 20000 \
    --n_sim 10000 \
    --output mrd.csv \
    --reference BSgenome.Hsapiens.UCSC.hg38 \

