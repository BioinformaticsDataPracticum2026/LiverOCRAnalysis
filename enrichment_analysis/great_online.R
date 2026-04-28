#!/usr/bin/env Rscript

# load packages quietly so the terminal output stays cleaner
suppressPackageStartupMessages({
  library(rGREAT)
  library(GenomicRanges)
})

# read command line inputs
args <- commandArgs(trailingOnly = TRUE)

# check that the user gave the expected arguments
if (length(args) != 3) {
  stop("Usage: Rscript enrichment_analysis/great_online.R <bed_file> <genome> <outdir>")
}

bed_file <- args[1]
genome <- args[2]
outdir <- args[3]

# make sure the BED file exists
if (!file.exists(bed_file)) {
  stop(paste("BED file not found:", bed_file))
}

# only allow the genomes supported in this script
if (!(genome %in% c("hg38", "mm10"))) {
  stop("Genome must be either hg38 or mm10")
}

# create output folder if needed
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# read BED file
d <- read.table(
  bed_file,
  sep = "\t",
  header = FALSE,
  quote = "",
  comment.char = "",
  stringsAsFactors = FALSE
)

# BED file needs at least chromosome, start, and end
if (ncol(d) < 3) {
  stop("BED file must have at least 3 columns: chr, start, end")
}

# keep only the main BED columns needed for GREAT
d <- d[, 1:3]
colnames(d) <- c("chr", "start", "end")

# make sure start and end are numeric positions
d$start <- as.integer(d$start)
d$end <- as.integer(d$end)

if (any(is.na(d$start)) || any(is.na(d$end))) {
  stop("BED start/end columns contain non-numeric values")
}

# chromosome sizes used by GREAT
if (genome == "hg38") {
  chr_sizes <- c(
    chr1=248956422, chr2=242193529, chr3=198295559, chr4=190214555,
    chr5=181538259, chr6=170805979, chr7=159345973, chr8=145138636,
    chr9=138394717, chr10=133797422, chr11=135086622, chr12=133275309,
    chr13=114364328, chr14=107043718, chr15=101991189, chr16=90338345,
    chr17=83257441, chr18=80373285, chr19=58617616, chr20=64444167,
    chr21=46709983, chr22=50818468, chrX=156040895, chrY=57227415
  )
} else {
  chr_sizes <- c(
    chr1=195471971, chr2=182113224, chr3=160039680, chr4=156508116,
    chr5=151834684, chr6=149736546, chr7=145441459, chr8=129401213,
    chr9=124595110, chr10=130694993, chr11=122082543, chr12=120129022,
    chr13=120421639, chr14=124902244, chr15=104043685, chr16=98207768,
    chr17=94987271, chr18=90702639, chr19=61431566, chrX=171031299,
    chrY=91744698
  )
}

# count regions before filtering
before_n <- nrow(d)

# remove regions with unsupported chromosomes or invalid coordinates
d <- d[d$chr %in% names(chr_sizes), ]
d <- d[d$start >= 0, ]
d <- d[d$end > d$start, ]
d <- d[d$end <= chr_sizes[d$chr], ]

# count regions after filtering
after_n <- nrow(d)

# print a small summary for checking
cat("Input BED:", bed_file, "\n")
cat("Genome:", genome, "\n")
cat("Original regions:", before_n, "\n")
cat("Filtered invalid regions:", before_n - after_n, "\n")
cat("Remaining regions:", after_n, "\n")

# stop if filtering removed everything
if (after_n == 0) {
  stop("No valid regions remain after filtering.")
}

# BED is 0-based half-open; GRanges is 1-based closed
gr <- GRanges(
  seqnames = d$chr,
  ranges = IRanges(start = d$start + 1, end = d$end)
)

# run online GREAT
job <- submitGreatJob(gr, genome = genome)

# get enrichment result tables from GREAT
tables <- getEnrichmentTables(job)

# make sure GO Biological Process results are available
if (!("GO Biological Process" %in% names(tables))) {
  stop("GO Biological Process table was not returned by GREAT.")
}

# save GO Biological Process enrichment results
write.csv(
  tables[["GO Biological Process"]],
  file.path(outdir, "gobp.csv"),
  row.names = FALSE
)

# save basic metadata about the run
writeLines(
  c(
    paste("bed_file:", bed_file),
    paste("genome:", genome),
    paste("original_regions:", before_n),
    paste("filtered_invalid_regions:", before_n - after_n),
    paste("remaining_regions:", after_n)
  ),
  file.path(outdir, "metadata.txt")
)

cat("DONE:", bed_file, "\n")