#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(readr)
  library(stringr)
})

option_list <- list(
  make_option(c("-d", "--dir"), type = "character", default = ".",
              help = "Directory containing the 11_int_class output files [default: %default]", metavar = "dir"),
  make_option(c("-o", "--out"), type = "character", default = NULL,
              help = "Output path for final annotation TSV [REQUIRED]", metavar = "file"),
  make_option(c("-p", "--prefix"), type = "character", default = "EchCo",
              help = "Sample prefix used in filenames (e.g. EchCo) [default: %default]", metavar = "prefix")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$out)) {
  print_help(opt_parser)
  stop("Error: --out is required.", call. = FALSE)
}

d <- opt$dir
pfx <- opt$prefix

# ---- 1. Master list: final lncRNA transcript IDs ----
final_ids <- read_lines(file.path(d, paste0(pfx, "_final_ids.txt")))
base <- tibble(ID = final_ids)
cat("Loaded", nrow(base), "final lncRNA IDs\n")

# ---- 2. Per-tool pass/fail summary ----
summary_df <- read_csv(file.path(d, paste0(pfx, "_intersection_summary.csv")),
                        show_col_types = FALSE)
base <- left_join(base, summary_df, by = "ID")

# ---- 3. FEELnc classification (best hit per transcript) ----
classes_df <- read_tsv(file.path(d, paste0(pfx, "_lncRNA_classes.txt")),
                        show_col_types = FALSE) %>%
  filter(isBest == 1) %>%
  rename(ID = lncRNA_transcript) %>%
  select(ID, lncRNA_gene, partnerRNA_gene, partnerRNA_transcript,
         direction, type, subtype, location, distance)

base <- left_join(base, classes_df, by = "ID")

# ---- 4. TE overlap summary (collapsed to one row per lncRNA) ----
te_path <- file.path(d, "TE_lncRNA_intersect.txt")
if (file.exists(te_path)) {
  te_df <- read_tsv(te_path, show_col_types = FALSE)
  te_summary <- te_df %>%
    group_by(LncRNA_ID) %>%
    summarise(
      TE_overlap_count   = n(),
      TE_types           = paste(sort(unique(TE_ID)), collapse = ";"),
      TE_total_overlap_bp = sum(Overlap),
      .groups = "drop"
    ) %>%
    rename(ID = LncRNA_ID)

  base <- left_join(base, te_summary, by = "ID")
} else {
  cat("Warning:", te_path, "not found, skipping TE overlap columns\n")
}

# Fill in "no TE overlap" defaults rather than leaving NA
base <- base %>%
  mutate(
    TE_overlap_count    = ifelse(is.na(TE_overlap_count), 0, TE_overlap_count),
    TE_types            = ifelse(is.na(TE_types), "None", TE_types),
    TE_total_overlap_bp = ifelse(is.na(TE_total_overlap_bp), 0, TE_total_overlap_bp)
  )

# ---- 5. Write final annotation table ----
write_tsv(base, opt$out, na = "NA")
cat("Final annotation table written to:", opt$out, "\n")
cat("Rows:", nrow(base), " Columns:", ncol(base), "\n")
