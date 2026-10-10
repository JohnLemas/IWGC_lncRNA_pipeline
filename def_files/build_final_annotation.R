#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(readr)
  library(stringr)
})

option_list <- list(
  make_option(c("-d", "--dir"), type = "character", default = ".",
              help = "Directory containing 11_int_class output files [default: %default]", metavar = "dir"),
  make_option(c("-o", "--out"), type = "character", default = NULL,
              help = "Output path for final annotation TSV [REQUIRED]", metavar = "file"),
  make_option(c("-p", "--prefix"), type = "character", default = "EchCo",
              help = "Sample prefix used in filenames [default: %default]", metavar = "prefix"),
  make_option(c("-t", "--tracking"), type = "character", default = NULL,
              help = "Optional path to tracking summary TSV.", metavar = "file"),
  make_option(c("--tpm"), type = "character", default = NULL,
              help = "Optional path to *.summary.tsv or *.matrix.tsv from parse_tpm_matrix.R.", metavar = "file")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$out)) stop("Error: --out is required.", call. = FALSE)

d <- opt$dir
pfx <- opt$prefix

# 1. Master IDs
final_ids <- read_lines(file.path(d, paste0(pfx, "_final_ids.txt")))
base <- tibble(ID = final_ids)
cat("Loaded", nrow(base), "final lncRNA IDs\n")

# 2. GTF Coordinates
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

# 3. FEELnc Classes
classes_raw <- read_tsv(file.path(d, paste0(pfx, "_lncRNA_classes.txt")), show_col_types = FALSE)

classes_df <- classes_raw %>%
  filter(isBest == 1) %>%
  mutate(distance = as.numeric(distance)) %>%
  rename(ID = lncRNA_transcript) %>%
  select(ID, lncRNA_gene, partnerRNA_gene, direction, type, subtype, location, distance) %>%
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

base <- left_join(base, classes_df, by = "ID")

# 4. TE Overlap
te_path <- file.path(d, "TE_lncRNA_intersect.txt")
if (file.exists(te_path)) {
  te_df <- read_tsv(te_path, show_col_types = FALSE)
  te_summary <- te_df %>%
    group_by(LncRNA_ID) %>%
    summarise(TE_types = paste(sort(unique(TE_ID)), collapse = ";"), .groups = "drop") %>%
    rename(ID = LncRNA_ID)

  base <- left_join(base, te_summary, by = "ID")
} else {
  base$TE_types <- NA
}
base$TE_types <- ifelse(is.na(base$TE_types), "None", base$TE_types)

# 5. Tracking Info
if (!is.null(opt$tracking)) {
  trk <- read_tsv(opt$tracking, col_types = cols(.default = "c"), na = c("", "NA"), show_col_types = FALSE)
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
}

# 6. TPM Summary/Matrix Info
tpm_cols <- NULL
if (!is.null(opt$tpm)) {
  tpm_df <- read_tsv(opt$tpm, show_col_types = FALSE)
  tpm_df <- tpm_df %>% rename(ID = transcript_id)
  tpm_cols <- setdiff(colnames(tpm_df), "ID")
  base <- left_join(base, tpm_df, by = "ID")
}

select_cols <- c("Chr", "Start", "Stop", "lncRNA_gene", "ID", "partnerRNA_gene",
                 "direction", "type", "subtype", "location", "Category", "TE_types")
if ("Tissue_Type" %in% colnames(base)) select_cols <- c(select_cols, "Tissue_Type")
if (!is.null(tpm_cols)) select_cols <- c(select_cols, tpm_cols)

final_df <- base %>% select(all_of(select_cols))
write_tsv(final_df, opt$out, na = "NA")
cat("Final annotation table written to:", opt$out, "\n")