#!/usr/bin/env Rscript

# Info 
# --------------------------------------------------
# This script will run longitudinal analysis.
#
# Input files are
# 1) final selected reporters
# 2) followup sample bam file and bai file
# 3) target.bed. 3 column with no headers.
# 4) blocklist.txt
# 
# Output files are
# 1) longitudinal mutation analysis csv file

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
    make_option(c("--reporters"), action="store", type="character", default=NULL, help="Input reporter list"),
    make_option(c("--sample_bam"), action="store", type="character", default=NULL, help="Input followup sample bam file"),
    make_option(c("--target_bed"), action="store", type="character", default=NULL, help="Input target bed file"),
    make_option(c("--blocklist"), action="store", type="character", default=NULL, help="Input blocklist"),
    make_option(c("--blist_type"), action="store", type="character", default='variant', help="Loci or substitution-specific blocklist type. Avaialble options are 'loci' or 'variant'."),
    make_option(c("--reads_threshold"), action="store", type="numeric", default=NULL, help="Informative reads threshold"),
    make_option(c("--pvalue_threshold"), action="store", type="numeric", default=NULL, help="Pvalue threshold for positivity call"),
    make_option(c("--vaf_threshold"), action="store", type="numeric", default=0.1, help="When calculating the background rate, the bases with higher than this VAF threshold will be ignored (real mutations/SNPs)."),
    make_option(c("--read_min_bq"), action="store", type="numeric", default=20, help="minimum base quality for a read to be counted"),
    make_option(c("--read_min_mq"), action="store", type="numeric", default=30, help="minimum mapping quality for a read to be counted"),
    make_option(c("--read_max_dp"), action="store", type="numeric", default=20000, help="Maximum read depth above which sampling will happen"),
    make_option(c("--n_sim"), action="store", type="numeric", default=10000, help="Number of Monte Carlo simulations"),
    make_option(c("--output"), action="store", type="character", default=NULL, help="Output result file"),
    make_option(c("--reference"), action="store", type="character", default="BSgenome.Hsapiens.UCSC.hg38", help="BSgenome.Hsapiens reference")
)

parseobj = OptionParser(option_list=option_list)
opt = parse_args(parseobj)

if (is.null(opt$reporters))
    stop(print_help(parseobj))

# Main
# --------------------------------------------------

# human genome reference
library(opt$reference, character.only = TRUE)

# blocklist option
if(opt$blist_type == "variant"){
    substitution <- TRUE
} else if(opt$blist_type == "loci"){
    substitution <- FALSE
}

# reporters
reporters <- read.table(opt$reporters, sep="\t", stringsAsFactors=F, header=T)
# blocklist
blist <- read.table(opt$blocklist, sep="\t", stringsAsFactors=F, header=F) %>% pull(V1)
# target bed
targets <- read.table(opt$target_bed, sep="\t", stringsAsFactors=F, header=F,
    col.names=c("chr", "start", "end"))

res <- ctDNAtools::test_ctDNA(
    reference = Hsapiens, 
    mutations = reporters,
    bam = opt$sample_bam,
    targets = targets,
    informative_reads_threshold = opt$reads_threshold,
    pvalue_threshold = opt$pvalue_threshold,
    vaf_threshold = opt$vaf_threshold,
    min_base_quality = opt$read_min_bq,
    min_mapq = opt$read_min_mq,
    max_depth = opt$read_max_dp,
    n_simulations = opt$n_sim,
    black_list = blist,
    substitution_specific = substitution)

write.table(res, opt$output, sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE)
