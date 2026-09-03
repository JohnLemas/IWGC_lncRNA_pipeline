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

# ---- 2. Chr/Start/Stop from the lncRNA GTF (transcript-level rows) ----
gtf_path <- file.path(d, paste0(pfx, "_lncRNA.gtf"))
gtf_df <- read_tsv(gtf_path, col_names = FALSE, comment = "#",
                    col_types = cols(.default = "c"), show_col_types = FALSE)
colnames(gtf_df)[1:9] <- c("Chr", "source", "feature", "Start", "Stop",
                            "score", "strand", "frame", "attributes")

coords_df <- gtf_df %>%
  filter(feature == "transcript") %>%
  mutate(ID = str_match(attributes, 'transcript_id "([^"]+)"')[, 2]) %>%
  select(Chr, Start, Stop, ID) %>%
  mutate(Start = as.integer(Start), Stop = as.integer(Stop))

base <- left_join(base, coords_df, by = "ID")

# ---- 3. FEELnc classification (best hit per transcript) ----
classes_df <- read_tsv(file.path(d, paste0(pfx, "_lncRNA_classes.txt")),
                        show_col_types = FALSE) %>%
  filter(isBest == 1) %>%
  rename(ID = lncRNA_transcript) %>%
  select(ID, lncRNA_gene, partnerRNA_gene, direction, type, subtype, location)

base <- left_join(base, classes_df, by = "ID")

# ---- 4. TE overlap types (collapsed to one row per lncRNA) ----
te_path <- file.path(d, "TE_lncRNA_intersect.txt")
if (file.exists(te_path)) {
  te_df <- read_tsv(te_path, show_col_types = FALSE)
  te_summary <- te_df %>%
    group_by(LncRNA_ID) %>%
    summarise(TE_types = paste(sort(unique(TE_ID)), collapse = ";"), .groups = "drop") %>%
    rename(ID = LncRNA_ID)

  base <- left_join(base, te_summary, by = "ID")
} else {
  cat("Warning:", te_path, "not found, skipping TE_types column\n")
  base$TE_types <- NA
}

base$TE_types <- ifelse(is.na(base$TE_types), "None", base$TE_types)

# ---- 5. Final column order ----
final_df <- base %>%
  select(Chr, Start, Stop, lncRNA_gene, ID, partnerRNA_gene,
         direction, type, subtype, location, TE_types)

# ---- 6. Write final annotation table ----
write_tsv(final_df, opt$out, na = "NA")
cat("Final annotation table written to:", opt$out, "\n")
cat("Rows:", nrow(final_df), " Columns:", ncol(final_df), "\n")