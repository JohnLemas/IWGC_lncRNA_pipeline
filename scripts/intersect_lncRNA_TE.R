#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(rtracklayer)
  library(GenomicRanges)
  library(dplyr)
  library(stringr)
})

# 1. Define command-line options
option_list <- list(
  make_option(c("-g", "--gff3"), type = "character", default = NULL,
              help = "Path to input EDTA TE annotation file (.gff3) [REQUIRED]", metavar = "file"),
  make_option(c("-t", "--gtf"), type = "character", default = NULL,
              help = "Path to input lncRNA annotation file (.gtf) [REQUIRED]", metavar = "file"),
  make_option(c("-o", "--outdir"), type = "character", default = ".",
              help = "Output directory path [default: %default]", metavar = "dir"),
  make_option(c("--out_intersect"), type = "character", default = "TE_lncRNA_intersect.txt",
              help = "Filename for final intersection summary table [default: %default]", metavar = "filename"),
  make_option(c("--out_te_bed"), type = "character", default = "TE_full.bed",
              help = "Filename for converted full TE BED [default: %default]", metavar = "filename"),
  make_option(c("--out_lnc_bed"), type = "character", default = "lncRNA_full.bed",
              help = "Filename for converted full lncRNA BED [default: %default]", metavar = "filename"),
  make_option(c("-s", "--strand_specific"), action = "store_true", default = FALSE,
              help = "Enforce strand-specific overlap (default is strand-ignorant) [default: %default]")
)

opt_parser <- OptionParser(
  option_list = option_list,
  usage = "Usage: %prog -g <TE_annotation.gff3> -t <lncRNA_annotation.gtf> [options]"
)
opt <- parse_args(opt_parser)

# 2. Validate required arguments
if (is.null(opt$gff3) || is.null(opt$gtf)) {
  print_help(opt_parser)
  stop("Error: Both --gff3 (-g) and --gtf (-t) parameters are required.", call. = FALSE)
}

# Ensure output directory exists
if (!dir.exists(opt$outdir)) {
  dir.create(opt$outdir, recursive = TRUE)
}

# Construct full output paths
intersect_out_path <- file.path(opt$outdir, opt$out_intersect)
te_bed_path        <- file.path(opt$outdir, opt$out_te_bed)
lnc_bed_path       <- file.path(opt$outdir, opt$out_lnc_bed)

cat("=== Processing Input Files in R ===\n")
cat("TE GFF3:     ", opt$gff3, "\n")
cat("lncRNA GTF:  ", opt$gtf, "\n")
cat("Output Dir:  ", opt$outdir, "\n\n")

# 3. Import annotations as GRanges objects
cat("[1/4] Importing TE GFF3 and lncRNA GTF...\n")
te_gr     <- import(opt$gff3, format = "gff3")
lncrna_gr <- import(opt$gtf, format = "gtf")

# 4. Export Full BED files (rtracklayer converts to 0-based BED)
cat("[2/4] Exporting full BED files...\n")
export(te_gr, te_bed_path, format = "BED")
export(lncrna_gr, lnc_bed_path, format = "BED")

# 5. Extract metadata labels
te_id <- ifelse(!is.na(mcols(te_gr)$Classification), mcols(te_gr)$Classification,
         ifelse(!is.na(mcols(te_gr)$Name), mcols(te_gr)$Name, "."))
mcols(te_gr) <- DataFrame(name = te_id)

lnc_id <- ifelse(!is.na(mcols(lncrna_gr)$transcript_id), mcols(lncrna_gr)$transcript_id, ".")
mcols(lncrna_gr) <- DataFrame(name = lnc_id)

# 6. Compute overlaps
cat("[3/4] Finding intersections with GenomicRanges...\n")
ignore_strand_flag <- !opt$strand_specific
hits <- findOverlaps(lncrna_gr, te_gr, ignore.strand = ignore_strand_flag)

# Calculate exact overlap length in base pairs
overlaps_bp <- width(pintersect(lncrna_gr[queryHits(hits)], te_gr[subjectHits(hits)]))

# 7. Format and export summary table
cat("[4/4] Writing final summary table...\n")
intersect_df <- data.frame(
  Chr          = as.character(seqnames(lncrna_gr[queryHits(hits)])),
  LncRNA_start = start(lncrna_gr[queryHits(hits)]) - 1,  # Convert back to 0-based for BED alignment
  LncRNA_end   = end(lncrna_gr[queryHits(hits)]),
  LncRNA_ID    = mcols(lncrna_gr[queryHits(hits)])$name,
  Strand       = as.character(strand(lncrna_gr[queryHits(hits)])),
  TE_start     = start(te_gr[subjectHits(hits)]) - 1,      # Convert back to 0-based
  TE_end       = end(te_gr[subjectHits(hits)]),
  TE_ID        = mcols(te_gr[subjectHits(hits)])$name,
  Overlap      = overlaps_bp
) %>% 
  distinct()

write.table(intersect_df, intersect_out_path, sep = "\t", quote = FALSE, row.names = FALSE)

cat("Complete! Results written to:", intersect_out_path, "\n")