#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicRanges)
  library(dplyr)
  library(stringr)
})

# 1. Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript intersect_lncRNA_TE.R <TE_annotation.gff3> <lncRNA_annotation.gtf>")
}

te_file     <- args[1]
lncrna_file <- args[2]

cat("=== Processing Input Files in R ===\n")

# 2. Import annotations as GRanges objects
cat("[1/4] Importing TE GFF3 and lncRNA GTF...\n")
te_gr     <- import(te_file, format = "gff3")
lncrna_gr <- import(lncrna_file, format = "gtf")

# 3. Export Full BED files (rtracklayer automatically converts to 0-based BED format)
cat("[2/4] Exporting full BED files...\n")
export(te_gr, "TE_full.bed", format = "BED")
export(lncrna_gr, "lncRNA_full.bed", format = "BED")

# 4. Extract clean BED annotations & export BED4 / BED6
# Clean TE metadata (extract Classification or Name)
te_id <- ifelse(!is.na(mcols(te_gr)$Classification), mcols(te_gr)$Classification,
         ifelse(!is.na(mcols(te_gr)$Name), mcols(te_gr)$Name, "."))
mcols(te_gr) <- DataFrame(name = te_id)

# Clean lncRNA metadata (extract transcript_id)
lnc_id <- ifelse(!is.na(mcols(lncrna_gr)$transcript_id), mcols(lncrna_gr)$transcript_id, ".")
mcols(lncrna_gr) <- DataFrame(name = lnc_id)

export(te_gr, "TE_clean.bed", format = "BED")
export(lncrna_gr, "lncRNA_clean.bed", format = "BED")

# 5. Compute overlaps
cat("[3/4] Finding intersections with GenomicRanges...\n")
hits <- findOverlaps(lncrna_gr, te_gr, ignore.strand = TRUE)

# Calculate exact overlap in base pairs
overlaps_bp <- width(pintersect(lncrna_gr[queryHits(hits)], te_gr[subjectHits(hits)]))

# 6. Format and export the final text table (converting back to 0-based BED coords for display)
cat("[4/4] Generating final summary table...\n")
intersect_df <- data.frame(
  Chr          = as.character(seqnames(lncrna_gr[queryHits(hits)])),
  LncRNA_start = start(lncrna_gr[queryHits(hits)]) - 1,  # Adjust to 0-based for BED consistency
  LncRNA_end   = end(lncrna_gr[queryHits(hits)]),
  LncRNA_ID    = mcols(lncrna_gr[queryHits(hits)])$name,
  Strand       = as.character(strand(lncrna_gr[queryHits(hits)])),
  TE_start     = start(te_gr[subjectHits(hits)]) - 1,      # Adjust to 0-based
  TE_end       = end(te_gr[subjectHits(hits)]),
  TE_ID        = mcols(te_gr[subjectHits(hits)])$name,
  Overlap      = overlaps_bp
) %>% 
  distinct()  # Equivalent to sort -u

write.table(intersect_df, "TE_lncRNA_intersect.txt", sep = "\t", quote = FALSE, row.names = FALSE)

cat("Complete! Results saved to TE_lncRNA_intersect.txt\n")