#!/usr/bin/env Rscript
# parse_tpm_matrix.R
#
# Aggregates per-sample TPM values (from stringtie_quantify.sb's -e -B output) into a
# transcript x sample TPM matrix, restricted to the predicted lncRNA candidates, with
# summary columns (mean/max TPM, samples detected) and an optional per-tissue mean TPM
# table when the metadata file has a second column. Plays the same role for TPM that
# parse_gffcompare_tracking.R plays for tissue/sample tracking.
#
# Called as: singularity run rbio_tools.sif tpm -q <quant_dir> -c <candidates.gtf> \
#                -m <species_metadata.txt> -o <out_prefix>

suppressMessages({
  library(optparse)
  library(dplyr)
  library(purrr)
  library(readr)
  library(tidyr)
  library(tibble)
  library(rtracklayer)
})

option_list <- list(
  make_option(c("-q", "--quant_dir"), type = "character",
              help = "Directory containing <sample_id>/t_data.ctab subfolders (from stringtie -e -B)"),
  make_option(c("-c", "--candidates_gtf"), type = "character",
              help = "species_candidates.gtf from stringtie_merge.sb (defines which transcript_ids to keep)"),
  make_option(c("-m", "--metadata"), type = "character",
              help = "species_metadata.txt: sample_id in col 1, optional tissue/condition in col 2"),
  make_option(c("-o", "--out_prefix"), type = "character",
              help = "Prefix for output files: <prefix>.matrix.tsv and, if tissue given, <prefix>.summary.tsv")
)
opt <- parse_args(OptionParser(option_list = option_list))

for (req in c("quant_dir", "candidates_gtf", "metadata", "out_prefix")) {
  if (is.null(opt[[req]])) stop("Missing required argument: --", req, call. = FALSE)
}

# ---- 1. transcript_ids to keep, from the lncRNA candidates GTF ----
gtf <- import(opt$candidates_gtf)
candidate_ids <- unique(na.omit(mcols(gtf)$transcript_id))
if (length(candidate_ids) == 0) {
  stop("No transcript_id values found in ", opt$candidates_gtf, " - check the file.", call. = FALSE)
}
cat("Loaded", length(candidate_ids), "candidate lncRNA transcript IDs from",
    basename(opt$candidates_gtf), "\n")

# ---- 2. sample IDs (+ optional tissue) from metadata ----
meta_lines <- read_lines(opt$metadata)
meta_lines <- meta_lines[nzchar(trimws(meta_lines))]
meta_split <- strsplit(meta_lines, "\t")
sample_ids <- trimws(vapply(meta_split, `[`, character(1), 1))

has_tissue <- length(meta_split) > 0 && all(lengths(meta_split) >= 2)
tissue_map <- NULL
if (has_tissue) {
  tissue_map <- tibble(
    sample_id = sample_ids,
    tissue = trimws(vapply(meta_split, `[`, character(1), 2))
  )
}

# ---- 3. read each sample's t_data.ctab, keep only candidate transcripts ----
read_sample_tpm <- function(sample_id) {
  ctab_path <- file.path(opt$quant_dir, sample_id, "t_data.ctab")
  if (!file.exists(ctab_path)) {
    warning("No t_data.ctab found for ", sample_id, ", skipping.", call. = FALSE)
    return(NULL)
  }
  read_tsv(ctab_path, show_col_types = FALSE) %>%
    transmute(transcript_id = t_name, !!sample_id := TPM) %>%
    filter(transcript_id %in% candidate_ids) %>%
    distinct(transcript_id, .keep_all = TRUE)
}

sample_frames <- compact(map(sample_ids, read_sample_tpm))
if (length(sample_frames) == 0) {
  stop("No t_data.ctab files were found - check --quant_dir and sample IDs in --metadata.",
       call. = FALSE)
}

# ---- 4. build the wide matrix: every candidate transcript_id x every sample found ----
matrix <- reduce(sample_frames, full_join, by = "transcript_id") %>%
  right_join(tibble(transcript_id = sort(candidate_ids)), by = "transcript_id") %>%
  mutate(across(-transcript_id, ~ replace_na(.x, 0)))

sample_cols <- setdiff(names(matrix), "transcript_id")

matrix <- matrix %>%
  mutate(
    mean_TPM           = rowMeans(across(all_of(sample_cols))),
    max_TPM             = do.call(pmax, as.list(select(matrix, all_of(sample_cols)))),
    n_samples_detected = rowSums(across(all_of(sample_cols), ~ .x > 0))
  )

matrix_out <- paste0(opt$out_prefix, ".matrix.tsv")
write_tsv(matrix, matrix_out)
cat("Wrote per-sample TPM matrix:", matrix_out, "\n")

# ---- 5. optional per-tissue mean TPM, for joining into the final annotation table ----
if (!is.null(tissue_map)) {
  present_samples <- intersect(sample_cols, tissue_map$sample_id)
  tissues <- unique(tissue_map$tissue[tissue_map$sample_id %in% present_samples])

  tissue_means <- map_dfc(tissues, function(t) {
    cols <- tissue_map$sample_id[tissue_map$tissue == t & tissue_map$sample_id %in% present_samples]
    tibble(!!paste0("mean_TPM_", t) := rowMeans(select(matrix, all_of(cols))))
  })

  summary <- bind_cols(
    select(matrix, transcript_id, mean_TPM, max_TPM, n_samples_detected),
    tissue_means
  )

  summary_out <- paste0(opt$out_prefix, ".summary.tsv")
  write_tsv(summary, summary_out)
  cat("Wrote per-tissue TPM summary:", summary_out, "\n")
} else {
  cat("No tissue/condition column detected in metadata (needs >= 2 columns) - ",
      "skipping per-tissue summary; only the sample-level matrix was written.\n", sep = "")
}
