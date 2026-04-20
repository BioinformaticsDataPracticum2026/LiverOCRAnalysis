#!/usr/bin/env Rscript

required_pkgs <- c("rGREAT", "GenomicRanges", "readr")

missing_pkgs <- required_pkgs[
  !sapply(required_pkgs, requireNamespace, quietly = TRUE)
]

cat("DEBUG: .libPaths() in Rscript:\n")
print(.libPaths())
cat("DEBUG: required_pkgs = ", paste(required_pkgs, collapse = ", "), "\n")
cat("DEBUG: missing_pkgs = ", paste(missing_pkgs, collapse = ", "), "\n")

if (length(missing_pkgs) > 0) {
  stop(
    paste0(
      "Missing required R packages: ",
      paste(missing_pkgs, collapse = ", "),
      "\nInstall them first, e.g.:\n",
      "install.packages('BiocManager')\n",
      "BiocManager::install(c('rGREAT','GenomicRanges'))\n",
      "install.packages('readr')\n"
    )
  )
}

suppressPackageStartupMessages({
  library(rGREAT)
  library(GenomicRanges)
  library(readr)
})

# ---- simple argument parsing without optparse ----

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  if (flag %in% args) {
    i <- match(flag, args)
    if (i < length(args)) {
      return(args[i + 1])
    }
  }
  default
}

bed     <- get_arg("--bed")
species <- get_arg("--species", default = "hg38")
outdir  <- get_arg("--outdir")
prefix  <- get_arg("--prefix")

if (is.null(bed) || is.null(outdir) || is.null(prefix)) {
  stop("Missing required arguments: --bed, --outdir, --prefix")
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

bed_df <- read.delim(bed, header = FALSE, stringsAsFactors = FALSE)
if (ncol(bed_df) < 3) {
  stop("BED file must have at least 3 columns")
}
bed_df <- bed_df[, 1:3]
colnames(bed_df) <- c("chr", "start", "end")

gr <- GRanges(
  seqnames = bed_df$chr,
  ranges = IRanges(start = bed_df$start + 1, end = bed_df$end)
)

job <- submitGreatJob(gr, species = species)
tables <- getEnrichmentTables(job)

writeLines(
  c(
    paste("bed_file:", bed),
    paste("species:", species),
    paste("table_count:", length(tables))
  ),
  con = file.path(outdir, paste0(prefix, "_run_info.txt"))
)

if (length(tables) == 0) {
  warning("No enrichment tables returned by GREAT.")
} else {
  for (nm in names(tables)) {
    clean_name <- gsub("[^A-Za-z0-9_]+", "_", nm)
    out_file <- file.path(outdir, paste0(prefix, "_", clean_name, ".tsv"))
    write_tsv(tables[[nm]], out_file)
  }
}