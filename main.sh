#!/usr/bin/env bash

# This is the Main module for the KAPA NHL Analysis Pipeline. 
# Includes FASTQ to UMI-deduped BAM to variant calls.
# The Main module should be run for the baseline, germline, and followup samples prior to longitudinal analysis. 


# Convert FASTQ to BAM
# --------------------------------------
gatk FastqToSam \ 
    -F1 SAMPLE_R1.fastq.gz \
    -F2 SAMPLE_R2.fastq.gz \
    -O SAMPLE_unmapped.bam \
    -SM SAMPLE


# Extract UMIs from BAM
# --------------------------------------
fgbio ExtractUmisFromBam \
    -i SAMPLE_unmapped.bam \
    -o SAMPLE_unmapped_umi_extracted.bam \
    -r 3M3S+T 3M3S+T \
    -t RX \
    -a true


# Perform Adapter Trimming and Quality Filtering
# --------------------------------------
# Convert BAM to FASTQ
gatk SamToFastq \
    -I SAMPLE_unmapped_umi_extracted.bam \
    -F SAMPLE_umi_extracted_R1.fastq \
    -F2 SAMPLE_umi_extracted_R2.fastq \
    --CLIPPING_ATTRIBUTE XT \
    --CLIPPING_ACTION 2

# Perform Adapter and Quality Trimming
fastp \
    -i SAMPLE_umi_extracted_R1.fastq \
    -o SAMPLE_umi_extracted_trimmed_R1.fastq \
    -I SAMPLE_umi_extracted_R2.fastq \
    -O SAMPLE_umi_extracted_trimmed_R2.fastq \
    -g -W 5 -q 20 -u 40 -x -3 -l 75 -c \
    -j fastp.json \
    -h fastp.html \
    -w NumProcessors &> SAMPLE_fastp.log


# Select a Subsample of Reads from a FASTQ File
# --------------------------------------
# This will subset to 80 million reads, or 40 million read pairs total.
seqtk sample \
    -s 12345 \
    SAMPLE_umi_extracted_trimmed_R1.fastq 40000000 > SAMPLE_umi_extracted_trimmed_subset_R1.fastq
seqtk sample \
    -s 12345 \
    SAMPLE_umi_extracted_trimmed_R2.fastq 40000000 > SAMPLE_umi_extracted_trimmed_subset_R2.fastq


# Map Subsampled Reads to the Reference Genome
# --------------------------------------
bwa mem \
    -R "@RG\\tID:A\\tDS:KAPA_TE\\tPL:ILLUMINA\\tLB:SAMPLE\\tSM:SAMPLE" \
    -t NumProcessors \
    -Y \
    ref.fa \
    SAMPLE_umi_extracted_trimmed_subset_R1.fastq \
    SAMPLE_umi_extracted_trimmed_subset_R1.fastq \
    | \
    samtools view -f2 -Sb > SAMPLE_umi_aligned.bam


# Add UMI Information to the Reads in BAM
# --------------------------------------
# Sort the downsampled bam by queryname for MergeBamAlignment
gatk SortSam \
    I=SAMPLE_umi_aligned.bam \
    O=SAMPLE_umi_aligned_qnsorted.bam \
    SORT_ORDER="queryname"

# merge the two BAM files 
# containing the UMI information - SAMPLE_unmapped_umi_extracted.bam and 
# the alignment coordinate information - SAMPLE_umi_aligned_qnsorted.bam
gatk MergeBamAlignment \
    --ATTRIBUTES_TO_RETAIN X0 \
    --ATTRIBUTES_TO_REMOVE NM \
    --ATTRIBUTES_TO_REMOVE MD \
    --ALIGNED_BAM SAMPLE_umi_aligned_qnsorted.bam \
    --UNMAPPED_BAM SAMPLE_unmapped_umi_extracted.bam \
    --OUTPUT SAMPLE_umi_extracted_aligned_merged.bam \
    --REFERENCE_SEQUENCE ref.fa \
    --SORT_ORDER queryname \
    --ALIGNED_READS_ONLY true \
    --MAX_INSERTIONS_OR_DELETIONS -1 \
    --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
    --ALIGNER_PROPER_PAIR_FLAGS true \
    --CLIP_OVERLAPPING_READS false


# Identify and Group Reads Originating from the Same Source Molecule 
# --------------------------------------
# BAM file outputs by GroupReadsbyUMI will have template coordinate sort order.
# The output SAMPLE_umi_grouped.bam will have SO:unsorted GO:query SS:unsorted:template-coordinate
fgbio GroupReadsByUmi \
    --input=SAMPLE_umi_extracted_aligned_merged.bam \
    --output=SAMPLE_umi_grouped.bam \
    --strategy=adjacency \
    --edits=1 \
    --min-map-q=20 \
    -t RX \
    -f SAMPLE_umi_group_data.tsv


# Calculate Consensus Sequence 
# --------------------------------------
fgbio CallMolecularConsensusReads \
    --input=SAMPLE_umi_grouped.bam \
    --output=SAMPLE_umi_consensus_unmapped.bam \
    --error-rate-post-umi 40 \
    --error-rate-pre-umi 45 \
    --output-per-base-tags false \
    --min-reads 1 \
    --max-reads 100 \
    --min-input-base-quality 20 \
    --read-name-prefix='consensus'


# Filter Consensus Reads
# --------------------------------------
fgbio FilterConsensusReads \
    --input=SAMPLE_umi_consensus_unmapped.bam \
    --output=SAMPLE_umi_consensus_unmapped_filtered.bam \
    --ref=ref.fa \
    --max-read-error-rate=0.025 \
    --min-base-quality=20 \
    --max-base-error-rate=0.1 \
    --max-no-call-fraction=0.1 \
    --min-reads=1 \
    --reverse-per-base-tags=true 


# Convert BAM to FASTQ
# --------------------------------------
# sort by queryname for MergeBamAlignment
gatk SortSam \
    I=SAMPLE_umi_consensus_unmapped_filtered.bam \
    O=SAMPLE_umi_consensus_unmapped_filtered_qnsorted.bam \
    SORT_ORDER="queryname"

# bam to fastq
gatk SamToFastq \
    -I SAMPLE_umi_consensus_unmapped_filtered_qnsorted.bam \
    -F SAMPLE_umi_consensus_unmapped_R1.fastq \
    -F2 SAMPLE_umi_consensus_unmapped_R2.fastq \
    --CLIPPING_ATTRIBUTE XT \
    --CLIPPING_ACTION 2


# Map Consensus Reads to the Reference Genome
# --------------------------------------
bwa mem \
    -R "@RG\\tID:A\\tDS:KAPA_TE\\tPL:ILLUMINA\\tLB:SAMPLE\\tSM:SAMPLE" \
    -t NumProcessors \
    -C \
    -Y \
    ref.fa \
    SAMPLE_umi_consensus_unmapped_R1.fastq \
    SAMPLE_umi_consensus_unmapped_R2.fastq \
    | \
    samtools view -bh - > SAMPLE_umi_consensus_mapped.bam


# Add UMI Information to the Consensus Reads in BAM
# --------------------------------------
# Sort by queryname for MergeBamAlignment
# The output SAMPLE_umi_deduped_sorted.bam will be coordinate sorted.
gatk SortSam \
    I=SAMPLE_umi_consensus_mapped.bam \
    O=SAMPLE_umi_consensus_mapped_qnsorted.bam \
    SORT_ORDER="queryname"

gatk MergeBamAlignment \
    --ATTRIBUTES_TO_RETAIN X0 \
    --ATTRIBUTES_TO_RETAIN RX \
    --ALIGNED_BAM SAMPLE_umi_consensus_mapped_qnsorted.bam \
    --UNMAPPED_BAM SAMPLE_umi_consensus_unmapped_filtered_qnsorted.bam \
    --OUTPUT SAMPLE_umi_deduped_sorted.bam \
    --REFERENCE_SEQUENCE ref.fa \
    --SORT_ORDER coordinate \
    --ADD_MATE_CIGAR true \
    --MAX_INSERTIONS_OR_DELETIONS -1 \
    --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
    --ALIGNER_PROPER_PAIR_FLAGS true \
    --CLIP_OVERLAPPING_READS false \
    --ADD_PG_TAG_TO_READS false \
    --UNMAPPED_READ_STRATEGY COPY_TO_TAG \
    --ORIENTATIONS FR \
    --CREATE_INDEX true


# clip overhangs 
# --------------------------------------
gatk SortSam \
    I=SAMPLE_umi_deduped.bam \
    O=SAMPLE_umi_deduped_qnsorted.bam \
    SORT_ORDER="queryname" 

fgbio ClipBam \
    --clip-overlapping-reads \
    -i SAMPLE_umi_deduped_qnsorted.bam \
    -o SAMPLE_umi_deduped_qnsorted_clipov.bam \
    -m SAMPLE_clip_ov_metrics.txt \
    -r ref.fa


# Sort and Index the Deduped clip overhang BAM
# --------------------------------------
gatk SortSam \
    I=SAMPLE_umi_deduped_qnsorted_clipov.bam \
    O=SAMPLE_umi_deduped_clipov_sorted.bam \
    SORT_ORDER="coordinate" \
    CREATE_INDEX=true


# Sort and Index the Non-deduped BAM
# --------------------------------------
gatk SortSam \
    I=SAMPLE_umi_aligned.bam \
    O=SAMPLE_umi_aligned_sorted.bam \
    SORT_ORDER="coordinate" \
    CREATE_INDEX=true


# Vardict variant calling
# --------------------------------------
vardict-java \
    -G ref.fa \
    -f AF_CUTOFF \
    -N SAMPLE \
    -b SAMPLE_umi_deduped_clipov_sorted.bam \
    -c 1 -S 2 -E 3 -g 4 DESIGN_capture_targets.bed | \
    teststrandbias.R | \
    var2vcf_valid.pl -N SAMPLE -E -f 0.005 > SAMPLE_vardict.vcf


# Variant annotation with dbSNP
# --------------------------------------
SnpSift annotate \
    -name "DBSNP_" \
    -info CAF,COMMON \
    dbsnp.vcf.gz > SAMPLE_vardict_annotated.vcf


# VCF to table
# --------------------------------------
gatk VariantsToTable \
    -V SAMPLE_vardict_annotated.vcf \
    --show-filtered \
    -F CHROM -F POS -F REF -F ALT -F ID -F FILTER \
    -F ADJAF -F AF -F BIAS -F DP -F VD -F DUPRATE \
    -F END -F HIAF -F MQ -F HICNT -F HICOV -F LSEQ \
    -F RSEQ -F NM -F PSTD -F QSTD -F QUAL -F REFBIAS \
    -F VARBIAS -F SBF -F SN -F TYPE \
    -F DP4 -F HRUN -F SB -F DBSNP_CAF -F DBSNP_COMMON \
    -GF AD -GF ALD -GF GT -GF RD \
    -O SAMPLE_vardict_annotated_vcf.txt
