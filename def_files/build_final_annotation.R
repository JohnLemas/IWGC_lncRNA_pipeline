#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(readr)
  library(stringr)
})

# This script now covers what used to be four separate steps:
#   1. build_final_annotation.R  - assemble Chr/Start/Stop, FEELnc classes, TE_types
#   2. classifications.R         - sort transcripts into structural categories
#                                   (previously written out as six standalone .txt
#                                   files that nothing downstream ever read; now a
#                                   single Category column instead)
#   3. add_tissue_annotation.R   - join in the tracking summary's Tissue_Type
#                                   (now optional, via --tracking)
#   4. parse_tpm_matrix.R output - join in mean/max TPM (and per-tissue mean TPM,
#                                   if present) from its *.summary.tsv/*.matrix.tsv
#                                   (now optional, via --tpm)

option_list <- list(
  make_option(c("-d", "--dir"), type = "character", default = ".",
              help = "Directory containing the 11_int_class output files [default: %default]", metavar = "dir"),
  make_option(c("-o", "--out"), type = "character", default = NULL,
              help = "Output path for final annotation TSV [REQUIRED]", metavar = "file"),
  make_option(c("-p", "--prefix"), type = "character", default = "EchCo",
              help = "Sample prefix used in filenames (e.g. EchCo) [default: %default]", metavar = "prefix"),
  make_option(c("-t", "--tracking"), type = "character", default = NULL,
              help = "Optional path to the tracking summary TSV from parse_gffcompare_tracking.R. If given, adds a Tissue_Type column.",
              metavar = "file"),
  make_option(c("--tpm"), type = "character", default = NULL,
              help = "Optional path to the *.summary.tsv (or *.matrix.tsv) from parse_tpm_matrix.R. If given, joins its TPM columns onto the final table.",
              metavar = "file")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$out)) {
  print_help(opt_parser)
  stop("Error: --out is required.", call. = FALSE)
}

if (!is.null(opt$tracking) && !file.exists(opt$tracking)) {
  stop(paste0("Error: --tracking file '", opt$tracking, "' does not exist."), call. = FALSE)
}

if (!is.null(opt$tpm) && !file.exists(opt$tpm)) {
  stop(paste0("Error: --tpm file '", opt$tpm, "' does not exist."), call. = FALSE)
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
classes_raw <- read_tsv(file.path(d, paste0(pfx, "_lncRNA_classes.txt")),
                         show_col_types = FALSE)

classes_df <- classes_raw %>%
  filter(isBest == 1) %>%
  mutate(distance = as.numeric(distance)) %>%
  rename(ID = lncRNA_transcript) %>%
  select(ID, lncRNA_gene, partnerRNA_gene, direction, type, subtype, location, distance)

# ---- 3b. Structural category (replaces classifications.R's six output files) ----
# Same rules classifications.R used to sort into
# LncRNA_antisense_exonic / intronic / upstream / downstream / intergenic /
# Bidirectional .txt files, kept here as one column since nothing downstream
# ever consumed those separate files directly.
classes_df <- classes_df %>%
  mutate(Category = case_when(
    direction == "antisense" & location == "exonic" ~ "Antisense_exonic",
    direction == "sense" & location == "intronic" ~ "Intronic",
    direction == "sense" & type == "intergenic" & distance <= 2000 & location == "upstream" ~ "Upstream",
    direction == "sense" & type == "intergenic" & distance <= 2000 & location == "downstream" ~ "Downstream",
    direction == "antisense" & type == "intergenic" & distance <= 2000 & location == "upstream" ~ "Bidirectional",
    type == "intergenic" & distance > 2000 ~ "Intergenic",
    TRUE ~ NA_character_
  )) %>%
  select(-distance)

n_uncat <- sum(is.na(classes_df$Category))
if (n_uncat > 0) {
  cat("Note:", n_uncat, "transcript(s) did not match any Category rule (left NA)\n")
}

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

# ---- 5. Optional Tissue_Type column (replaces add_tissue_annotation.R) ----
if (!is.null(opt$tracking)) {
  cat("Adding Tissue_Type column from", opt$tracking, "\n")

  trk <- read_tsv(opt$tracking, col_types = cols(.default = "c"),
                   na = c("", "NA"), show_col_types = FALSE)

  if (!all(c("Predicted_Gene_ID", "Tissue_Types") %in% colnames(trk))) {
    stop("Error: tracking summary must contain 'Predicted_Gene_ID' and 'Tissue_Types' columns.",
         call. = FALSE)
  }

  # Collapse duplicate Predicted_Gene_ID rows by unioning their tissue sets.
  # The gffcompare .tracking file can split one candidate transcript across
  # multiple lines, each carrying only a partial slice of its tissue evidence.
  union_tissues <- function(vals) {
    parts <- unlist(str_split(vals[!is.na(vals)], ","))
    parts <- str_trim(parts)
    parts <- parts[parts != ""]
    paste(sort(unique(parts)), collapse = ",")
  }

  agg <- trk %>%
    group_by(Predicted_Gene_ID) %>%
    summarise(Tissue_Type = union_tissues(Tissue_Types), .groups = "drop") %>%
    rename(ID = Predicted_Gene_ID)

  base <- left_join(base, agg, by = "ID")

  n_missing <- sum(is.na(base$Tissue_Type) | base$Tissue_Type == "")
  if (n_missing > 0) {
    cat("Warning:", n_missing,
        "annotation row(s) had no matching tracking data (Tissue_Type left blank)\n")
  }
} else {
  cat("No --tracking file given, skipping Tissue_Type column\n")
}

# ---- 5b. Optional TPM columns (from parse_tpm_matrix.R's output) ----
# parse_tpm_matrix.R's own header note says this join belongs here, "the same
# way it already joins in the tracking summary" - this does that. Works with
# either file it writes: pass the *.summary.tsv for mean_TPM/max_TPM/
# n_samples_detected (+ mean_TPM_<tissue> columns, if metadata had a tissue
# column) or the *.matrix.tsv for the full per-sample TPM values instead.
tpm_cols <- NULL
if (!is.null(opt$tpm)) {
  cat("Adding TPM columns from", opt$tpm, "\n")

  tpm_df <- read_tsv(opt$tpm, show_col_types = FALSE)

  if (!"transcript_id" %in% colnames(tpm_df)) {
    stop("Error: --tpm file must have a 'transcript_id' column (pass the *.summary.tsv or *.matrix.tsv from parse_tpm_matrix.R).",
         call. = FALSE)
  }

  tpm_df <- tpm_df %>% rename(ID = transcript_id)
  tpm_cols <- setdiff(colnames(tpm_df), "ID")

  base <- left_join(base, tpm_df, by = "ID")

  n_missing_tpm <- sum(!complete.cases(base[tpm_cols]))
  if (n_missing_tpm > 0) {
    cat("Warning:", n_missing_tpm,
        "annotation row(s) had no matching TPM data\n")
  }
} else {
  cat("No --tpm file given, skipping TPM columns\n")
}

# ---- 6. Final column order ----
select_cols <- c("Chr", "Start", "Stop", "lncRNA_gene", "ID", "partnerRNA_gene",
                  "direction", "type", "subtype", "location", "Category", "TE_types")
if ("Tissue_Type" %in% colnames(base)) select_cols <- c(select_cols, "Tissue_Type")
if (!is.null(tpm_cols)) select_cols <- c(select_cols, tpm_cols)

final_df <- base %>% select(all_of(select_cols))

if (nrow(final_df) != length(final_ids)) {
  cat("Warning: row count changed during processing (", length(final_ids), "->", nrow(final_df),
      "). Check for duplicate IDs.\n")
}

# ---- 7. Write final annotation table ----
write_tsv(final_df, opt$out, na = "NA")
cat("Final annotation table written to:", opt$out, "\n")
cat("Rows:", nrow(final_df), " Columns:", ncol(final_df), "\n")
