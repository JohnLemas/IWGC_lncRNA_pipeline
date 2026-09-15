#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(stringr)
})

# 1. Define command-line options
option_list <- list(
  make_option(c("-a", "--annotation"), type = "character", default = NULL,
              help = "Path to the final annotation TSV from build_final_annotation.R [REQUIRED]", metavar = "file"),
  make_option(c("-t", "--tracking"), type = "character", default = NULL,
              help = "Path to the tracking summary TSV from parse_gffcompare_tracking.R [REQUIRED]", metavar = "file"),
  make_option(c("-o", "--out"), type = "character", default = NULL,
              help = "Path for the annotated output TSV [REQUIRED]", metavar = "file")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# 2. Validate input arguments
required <- c("annotation", "tracking", "out")
missing <- required[vapply(required, function(x) is.null(opt[[x]]), logical(1))]
if (length(missing) > 0) {
  print_help(opt_parser)
  stop(paste0("Error: missing required argument(s): ",
              paste(paste0("--", missing), collapse = ", ")), call. = FALSE)
}

for (f in c(opt$annotation, opt$tracking)) {
  if (!file.exists(f)) {
    stop(paste0("Error: input file '", f, "' does not exist."), call. = FALSE)
  }
}

# 3. Read inputs
# NOTE: read everything as character and treat ONLY "" / "NA" as missing.
# build_final_annotation.R deliberately writes the literal string "None" into
# TE_types to mean "no TE overlap" - it must NOT be coerced to a missing value.
ann <- read_tsv(opt$annotation, col_types = cols(.default = "c"),
                na = c("", "NA"), show_col_types = FALSE)
trk <- read_tsv(opt$tracking, col_types = cols(.default = "c"),
                na = c("", "NA"), show_col_types = FALSE)

if (!"ID" %in% colnames(ann)) {
  stop("Error: annotation file has no 'ID' column to join on.", call. = FALSE)
}
if (!all(c("Predicted_Gene_ID", "Tissue_Types") %in% colnames(trk))) {
  stop("Error: tracking summary must contain 'Predicted_Gene_ID' and 'Tissue_Types' columns.",
       call. = FALSE)
}

cat("Loaded", nrow(ann), "annotation rows and", nrow(trk), "tracking rows\n")

# 4. Collapse duplicate Predicted_Gene_ID rows by unioning their tissue sets.
# The gffcompare .tracking file can split one candidate transcript across
# multiple lines, each carrying only a partial slice of its tissue evidence.
# Taking the union here keeps every tissue and prevents the join below from
# duplicating annotation rows.
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

cat("Collapsed tracking summary to", nrow(agg), "unique transcript IDs\n")

# 5. Join tissue info onto the annotation table
final_df <- left_join(ann, agg, by = "ID")

n_missing <- sum(is.na(final_df$Tissue_Type) | final_df$Tissue_Type == "")
if (n_missing > 0) {
  cat("Warning:", n_missing,
      "annotation row(s) had no matching tracking data (Tissue_Type left blank)\n")
}

if (nrow(final_df) != nrow(ann)) {
  cat("Warning: row count changed during join (", nrow(ann), "->", nrow(final_df),
      "). Check for duplicate IDs in the annotation file.\n")
}

# 6. Write annotated table
write_tsv(final_df, opt$out, na = "")

cat("Annotated table written to:", opt$out, "\n")
cat("Rows:", nrow(final_df), " Columns:", ncol(final_df), "\n")
