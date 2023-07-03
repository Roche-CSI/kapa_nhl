#!/usr/bin/env Rscript

# Info
# --------------------------------------------------
# This script will filter and remove candidate reporter variants based on filtering criteria based on Vardict output.
# It will then remove candidates that are present in germline.
#
# Input files are 
# 1) germline bam file.
# 2) followup bam file .
# 3) vcf table reporter list.
#
# Output files are 
# 1) a filtered list of final selected reporter candidates.
# 2) a filtered list of final selected reporter candidates in vcf format.
# 3) an annotated list of all reporter candidates. 

suppressPackageStartupMessages({
    library(purrr)
    library(tidyr)
    library(dplyr)
    library(Rsamtools)
    library(optparse)
})

# Options
# --------------------------------------------------
option_list <- list(
    make_option(c("--filter_reporters"), type="logical", default=TRUE, help="Filters out reporters using metrics from Vardict. Default=TRUE. Set to FALSE to turn off."),
    make_option(c("--remove_snp"), type="logical", default=TRUE, help="Filters out reporters if annotated with DBSNP_COMMON. Default=TRUE. Set to FALSE to turn off."),
    make_option(c("--reporters"), action="store", type="character", default=NULL, help="Input candidate reporter list"),
    make_option(c("--germline"), action="store", type="character", default=NULL, help="Input germline bam file"),
    make_option(c("--followup"), action="store", type="character", default=NULL, help="Input followup bam file"),
    make_option(c("--selected"), action="store", type="character", default=NULL, help="Output txt file for final selected reporters after filtering"),
    make_option(c("--selected_vcf"), action="store", type="character", default=NULL, help="Output vcf file for final selected reporters after filtering"),
    make_option(c("--all"), action="store", type="character", default=NULL, help="Output all reporters annotated with read count information"),
    make_option(c("--germline_cutoff"), action="store", type="numeric", default=NULL, help="Reporters with germline AF above this cutoff will be considered germline-positive"),
    make_option(c("--blocklist"), action="store", type="character", default=NULL, help="blocklist file, used for annotating reporters"),
    make_option(c("--min_af"), action="store", type="numeric", default=0.005, help="Minimum AF for variant filtering"),
    make_option(c("--max_af"), action="store", type="numeric", default=0.35, help="Maximum AF for variant filtering"),
    make_option(c("--min_dp"), action="store", type="numeric", default=1000, help="Minimum depth for variant filtering"),
    make_option(c("--min_vd"), action="store", type="numeric", default=15, help="Minimum alt depth for variant filtering"),
    make_option(c("--min_mq"), action="store", type="numeric", default=55, help="Minimum mapping quality for for variant filtering"),
    make_option(c("--min_qual"), action="store", type="numeric", default=45, help="Minimum average base quality for variant filtering"),
    make_option(c("--min_sbf"), action="store", type="numeric", default=0, help="Minimum Strand Bias Fisher p-value for variant filtering"),
    make_option(c("--max_nm"), action="store", type="numeric", default=6, help="Maximum mean mismatches in reads for variant filtering"),
    make_option(c("--read_min_bq"), action="store", type="numeric", default=30, help="Minimum base quality for read counts to be considered"),
    make_option(c("--read_min_mq"), action="store", type="numeric", default=30, help="Minimum mapping quality for read counts to be considered"),
    make_option(c("--read_max_dp"), action="store", type="numeric", default=20000, help="Maximum read depth above which sampling will happen")
)

parseobj = OptionParser(option_list=option_list)
opt = parse_args(parseobj)

if (is.null(opt$reporters))
    stop(print_help(parseobj))

# Main
# --------------------------------------------------
# Read input vcf table from GATK VariantsToTable
vcf <- read.table(opt$reporters, header=TRUE, sep='\t') 

# remove samplename from GT columns
names(vcf) <- gsub(".*\\.GT", "GT", names(vcf))
names(vcf) <- gsub(".*\\.AD", "AD", names(vcf))
names(vcf) <- gsub(".*\\.ALD", "ALD", names(vcf))
names(vcf) <- gsub(".*\\.RD", "RD", names(vcf))

vcf <- vcf %>% 
    dplyr::mutate_if(is.factor, as.character) %>%
    dplyr::filter(nchar(REF) == 1 & nchar(ALT) == 1) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
        AF = as.numeric(AF), 
        DP = as.numeric(DP),
        QUAL = as.numeric(QUAL),
        VD = as.numeric(VD),
        baseline_ref = DP - VD,
        baseline_alt = VD, 
        baseline_af = AF,
        variant_id = paste0(CHROM, "_", POS, "_", REF, "_", ALT)
        )

# Remove common SNP
if (opt$remove_snp) {
    vcf <- vcf %>% dplyr::filter(!DBSNP_COMMON %in% '1')
}

# Perform filtering
if (opt$filter_reporters) {
    required_cols <- c("AF","DP","VD","QUAL","MQ","SBF")
    if (length(setdiff(required_cols, colnames(vcf))) > 0) {
        stop("Missing required Vardict variant fields.")
    }
    vcf <- vcf %>% 
        dplyr::filter(
            FILTER %in% c("PASS"),
            AF > opt$min_af & AF < opt$max_af,
            DP > opt$min_dp,
            VD > opt$min_vd, 
            MQ > opt$min_mq,
            QUAL > opt$min_qual,
            SBF >= opt$min_sbf,
            NM <= opt$max_nm
            ) %>%
        dplyr::select(CHROM, POS, ID, REF, ALT, QUAL, FILTER, AF, DP, GT, 
            variant_id, baseline_ref, baseline_alt, baseline_af)
} else {
    vcf <- vcf %>% 
        dplyr::select(CHROM, POS, ID, REF, ALT, QUAL, FILTER, AF, DP, GT,
            variant_id, baseline_ref, baseline_alt, baseline_af)
}

# Use ctDNAtools get_read_counts function to obtain read counts for a given chrom-pos in a bam file.
# The same method from the ctDNAtools package will be used to obtain read counts in 
# germline bam file and followup bam file for each reporter variant.
# The read counts are then used to annotate the reporter summary file. 
get_read_counts <- function(chr, pos, bam, tag = "", min_base_quality = 20, max_depth = 1e+05, min_mapq = 30) {
    gr <- GenomicRanges::GRanges(chr, IRanges::IRanges(pos, pos))

    if (tag == "") {
        sbp <- Rsamtools::ScanBamParam(which = gr)
    } else {
        sbp <- Rsamtools::ScanBamParam(which = gr, tagFilter = list(RG = tag))
    }

    pileup_param <- Rsamtools::PileupParam(
    max_depth = max_depth, min_base_quality = min_base_quality,
    min_mapq = min_mapq, distinguish_strands = FALSE, include_deletions = FALSE, include_insertions = FALSE
    )

    p <- Rsamtools::pileup(bam, scanBamParam = sbp, pileupParam = pileup_param)

    cbase <- ifelse("C" %in% p$nucleotide, p[p$nucleotide == "C", "count"], 0)
    gbase <- ifelse("G" %in% p$nucleotide, p[p$nucleotide == "G", "count"], 0)
    abase <- ifelse("A" %in% p$nucleotide, p[p$nucleotide == "A", "count"], 0)
    tbase <- ifelse("T" %in% p$nucleotide, p[p$nucleotide == "T", "count"], 0)

    return(list(A = abase, C = cbase, G = gbase, T = tbase))
}

# For every vcf reporter variant, get germline sample read counts
germline_readcounts <- purrr::pmap_dfr(list(
    vcf$CHROM, vcf$POS, vcf$REF, vcf$ALT
    ), function(chr, pos, ref, alt) {
        counts <- get_read_counts(
            chr = chr, pos = pos, 
            bam = opt$germline, tag = "", min_base_quality = opt$read_min_bq,
            min_mapq = opt$read_min_mq, max_depth = opt$read_max_dp
            )
        dat <- data.frame(germline_ref = counts[[ref]], germline_alt = counts[[alt]]) %>% 
            rowwise() %>% 
            mutate(germline_dp = germline_alt + germline_ref,
                    germline_af = germline_alt / germline_dp) 
        return(dat)
    })

# For every vcf reporter variant, get followup sample read counts.
followup_readcounts <- purrr::pmap_dfr(list(
    vcf$CHROM, vcf$POS, vcf$REF, vcf$ALT
    ), function(chr, pos, ref, alt) {
        counts <- get_read_counts(
            chr = chr, pos = pos, 
            bam = opt$followup, tag = "", min_base_quality = opt$read_min_bq,
            min_mapq = opt$read_min_mq, max_depth = opt$read_max_dp
            )
        dat <- data.frame(followup_ref = counts[[ref]], followup_alt = counts[[alt]]) %>% 
            rowwise() %>% 
            mutate(followup_dp = followup_alt + followup_ref,
                    followup_af = followup_alt / followup_dp) 
        return(dat)
    })

# Annotation for reporters in blocklist
blist <- read.table(opt$blocklist, sep="\t", stringsAsFactors=F, header=F)

# Combine vcf entries with germline and followup readcounts
merged <- dplyr::bind_cols(vcf, germline_readcounts) %>% 
    dplyr::bind_cols(., followup_readcounts) %>% 
    dplyr::mutate(
        germline_status =  dplyr::case_when(
                            germline_af >= opt$germline_cutoff ~ "germline",
                            TRUE ~ NA_character_),
        blocklist_status =  dplyr::case_when(
                            variant_id %in% blist$V1 ~ "blocklist",
                            TRUE ~ NA_character_)
    )

# final clean reporters after germline subtraction
vcf_keep <- merged %>% 
    dplyr::filter(
        germline_af < opt$germline_cutoff,
        !blocklist_status %in% "blocklist"
    )

# Output final selected reporters
write.table(vcf_keep, opt$selected, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

# Save all the reporters (prior to germline subtraction)
write.table(merged, opt$all, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

# Prepare columns for output vcf file containing final selected reporters
# Include DP and AF INFO
vcf_dat <- vcf_keep %>% 
    dplyr::mutate(
        INFO = paste0("DP=", DP, ";", "AF=", AF),
        FORMAT = "GT:DP:AF",
        SAMPLE = paste0(GT, ":", DP, ":", AF)
        ) %>% 
    dplyr::select(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, SAMPLE) 

# vcf header for vardict
vcf_header <- paste0(
'##fileformat=VCFv4.2', "\n",
'##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">', "\n",
'##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">', "\n",
'##FILTER=<ID=q22.5,Description="Mean Base Quality Below 22.5">', "\n",
'##FILTER=<ID=Q10,Description="Mean Mapping Quality Below 10">', "\n",
'##FILTER=<ID=p8,Description="Mean Position in Reads Less than 8">', "\n",
'##FILTER=<ID=SN1.5,Description="Signal to Noise Less than 1.5">', "\n",
'##FILTER=<ID=Bias,Description="Strand Bias">', "\n",
'##FILTER=<ID=pSTD,Description="Position in Reads has STD of 0">', "\n",
'##FILTER=<ID=d3,Description="Total Depth < 3">', "\n",
'##FILTER=<ID=v2,Description="Var Depth < 2">', "\n",
'##FILTER=<ID=f0.001,Description="Allele frequency < 0.001">', "\n",
'##FILTER=<ID=MSI12,Description="Variant in MSI region with 12 non-monomer MSI or 13 monomer MSI">', "\n",
'##FILTER=<ID=NM5.25,Description="Mean mismatches in reads >= 5.25, thus likely false positive">', "\n",
'##FILTER=<ID=InGap,Description="The variant is in the deletion gap, thus likely false positive">', "\n",
'##FILTER=<ID=InIns,Description="The variant is adjacent to an insertion variant">', "\n",
'##FILTER=<ID=Cluster0bp,Description="Two variants are within 0 bp">', "\n",
'##FILTER=<ID=LongMSI,Description="The somatic variant is flanked by long A/T (>=14)">', "\n",
'##FILTER=<ID=AMPBIAS,Description="Indicate the variant has amplicon bias.">', "\n",
'##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">', "\n",
'##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total Depth">', "\n",
'##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">', "\n"
)

# write out vcf file with header
file.create(opt$selected_vcf)
cat(vcf_header, file = opt$selected_vcf, append = TRUE)
cat(paste0("#", paste(names(vcf_dat), collapse = "\t"), "\n"), file = opt$selected_vcf, append = TRUE)
write.table(vcf_dat, opt$selected_vcf, row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE, append = TRUE)
