#!/usr/bin/env Rscript

# Info
# --------------------------------------------------
# This script will create blocklist.txt for use with ctDNAtools.

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
    make_option(c("--panel_background"), action="store", type="character", default=NULL, 
        help="RDS file panel background file containing the information for background mutations"),
    make_option(c("--blocklist"), action="store", type="character", default=NULL, 
        help="Output blocklist file used for filtering reporter candidates"),
    make_option(c("--vaf_quantile"), action="store", type="numeric", default=NULL, 
        help="The quantile of mean VAF above which the loci are considered noisy"),
    make_option(c("--min_samples_one_read"), action="store", type="numeric", default=NULL, 
        help="Loci that at least this number of samples exhibit at least one non-reference reads are considered noisy"),
    make_option(c("--min_samples_two_reads"), action="store", type="numeric", default=NULL, 
        help="Loci that at least this number of samples exhibit at least two non-reference reads are considered noisy"),
    make_option(c("--min_samples_n_reads"), action="store", type="numeric", default=NULL, 
        help="Loci that at least this number of samples exhibit at least n non-reference reads are considered noisy"),
    make_option(c("--n_reads"), action="store", type="numeric", default=NULL, 
        help="The number of reads to use in the option 'blocklist_min_samples_n_reads'")
)

parseobj = OptionParser(option_list=option_list)
opt = parse_args(parseobj)

if (is.null(opt$panel_background))
    stop(print_help(parseobj))

# Main
# --------------------------------------------------
# load panel background 
panel_bg <- readRDS(opt$panel_background)

# Creating blocklist
blocklist <- ctDNAtools::create_black_list(panel_bg,
          mean_vaf_quantile = opt$vaf_quantile,
          min_samples_one_read = opt$min_samples_one_read, 
          min_samples_two_reads = opt$min_samples_two_reads, 
          min_samples_n_reads = opt$min_samples_n_reads,
          n_reads = opt$n_reads)

write.table(blocklist, opt$blocklist,  sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
