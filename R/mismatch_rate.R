#!/usr/bin/env Rscript

# Info 
# --------------------------------------------------
# This script will calculate mismatch rate.
#
# Input files are
# 1) bam file
# 2) target.bed. 3 column with no headers.
# 3) blocklist.txt
# 
# Output files are
# 1) mismatch_rate.csv 

suppressPackageStartupMessages({
    library(purrr)
    library(tidyr)
    library(dplyr)
    library(optparse)
    library(ctDNAtools)
})

# Options
# --------------------------------------------------
option_list <- list(
    make_option(c("--sample_bam"), action="store", type="character", default=NULL, help="Input followup sample bam file"),
    make_option(c("--target_bed"), action="store", type="character", default=NULL, help="Input target bed file"),
    make_option(c("--blocklist"), action="store", type="character", default=NULL, help="Input blocklist file"),
    make_option(c("--blist_type"), action="store", type="character", default='variant', help="Loci or substitution-specific blocklist type. Avaialble options are 'loci' or 'variant'."),
    make_option(c("--vaf_threshold"), action="store", type="numeric", default=0.1, help="When calculating the background rate, the bases with higher than this VAF threshold will be ignored (real mutations/SNPs)."),
    make_option(c("--read_min_bq"), action="store", type="numeric", default=20, help="minimum base quality for a read to be counted"),
    make_option(c("--read_min_mq"), action="store", type="numeric", default=30, help="minimum mapping quality for a read to be counted"),
    make_option(c("--read_max_dp"), action="store", type="numeric", default=20000, help="Maximum read depth above which sampling will happen"),
    make_option(c("--output"), action="store", type="character", default=NULL, help="Output result file"),
    make_option(c("--reference"), action="store", type="character", default="BSgenome.Hsapiens.UCSC.hg38", help="BSgenome.Hsapiens reference")
)

parseobj = OptionParser(option_list=option_list)
opt = parse_args(parseobj)

if (is.null(opt$sample_bam))
    stop(print_help(parseobj))

# Main
# --------------------------------------------------

# human genome reference
library(opt$reference, character.only = TRUE)

# target bed
targets <- read.table(opt$target_bed, sep="\t", stringsAsFactors=F, header=F)
names(targets) <- c('chr', 'start', 'end')

if (is.null(opt$blocklist)) {
    print("No blacklist")
    mismatch_rate <- ctDNAtools::get_background_rate(
        bam = opt$sample_bam, 
        targets = targets, 
        reference = Hsapiens, 
        min_base_quality = opt$read_min_bq, 
        min_mapq = opt$read_min_mq,
        vaf_threshold = opt$read_max_dp)
} else {
    print("Use blacklist")
    # blocklist option
    if (opt$blist_type == "variant"){
        substitution <- TRUE
    } else if(opt$blist_type == "loci"){
        substitution <- FALSE
    }
    # blocklist
    blist <- read.table(opt$blocklist, sep="\t", stringsAsFactors=F, header=F) %>% pull(V1)
    mismatch_rate <- ctDNAtools::get_background_rate(
        bam = opt$sample_bam, 
        targets = targets, 
        reference = Hsapiens, 
        min_base_quality = opt$read_min_bq, 
        min_mapq = opt$read_min_mq,
        vaf_threshold = opt$read_max_dp,
        black_list = blist,
        substitution_specific = substitution
        )
}

write.table(mismatch_rate, opt$output, sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE)
