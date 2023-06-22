#!/usr/bin/env Rscript

# Info
# --------------------------------------------------
# This script will create panel_background.RDS for use with ctDNAtools.

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
    make_option(c("--bam_list"), action="store", type="character", default=NULL, 
        help="CSV file containing path to bam files, no column headers."),
    make_option(c("--panel_background"), action="store", type="character", default=NULL, 
        help="output RDS panel background file containing the information for background mutations."),
    make_option(c("--target_bed"), action="store", type="character", default=NULL, 
        help="Target bed file, tab separated, 3 columns no headers."),
    make_option(c("--blist_type"), action="store", type="character", default='variant', 
        help="Loci or substitution-specific blocklist type. Avaialble options are 'loci' or 'variant'."),
    make_option(c("--reference"), action="store", type="character", default="BSgenome.Hsapiens.UCSC.hg38", 
        help="BSgenome.Hsapiens reference.")
)

parseobj = OptionParser(option_list=option_list)
opt = parse_args(parseobj)

if (is.null(opt$bam_list))
    stop(print_help(parseobj))

# Main
# --------------------------------------------------
# blocklist option
if (opt$blist_type == "variant") {
    substitution <- TRUE
} else if (opt$blist_type == "loci") {
    substitution <- FALSE
}

# human genome reference
library(opt$reference, character.only = TRUE)

# Read input files
bams <- read.csv(opt$bam_list, header = FALSE) %>% pull(V1)
target_bed <- read.table(opt$target_bed, header = FALSE, sep="\t")
names(target_bed) <- c("chr", "start", "end")

# Create background panel
panel_bg <- ctDNAtools::create_background_panel(
    bam_list = bams, 
    targets = target_bed, 
    reference = Hsapiens, 
    substitution_specific = substitution)

saveRDS(panel_bg, opt$panel_background)
