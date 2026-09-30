#!/bin/usr/env Rscript
# parse_tpm_matrix.R

suppressMessages({
  library(rtracklayer)
  library(optparse)
  library(dplyr)
  library(purrr)
  library(readr)
  library(tidyr)
  library(tibble)
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

# ---- 3. read each sample's quantification GTF, keep only candidate transcripts ----
read_sample_tpm <- function(sample_id) {
  gtf_path <- file.path(opt$quant_dir, sample_id, paste0(sample_id, "_quant.gtf"))
  if (!file.exists(gtf_path)) {
    warning("No quantification GTF found for ", sample_id, " at ", gtf_path,
            ", skipping.", call. = FALSE)
    return(NULL)
  }
  gr <- import(gtf_path)
  if (!"TPM" %in% names(mcols(gr))) {
    warning("Quantification GTF for ", sample_id, " has no TPM attribute, skipping.",
            call. = FALSE)
    return(NULL)
  }
  
  as_tibble(mcols(gr)) %>%
    dplyr::filter(type == "transcript") %>%
    dplyr::mutate(TPM = as.numeric(TPM)) %>%
    dplyr::transmute(transcript_id = transcript_id, !!sample_id := TPM) %>%
    dplyr::filter(transcript_id %in% candidate_ids) %>%
    dplyr::distinct(transcript_id, .keep_all = TRUE)
}

sample_frames <- purrr::compact(map(sample_ids, read_sample_tpm))
if (length(sample_frames) == 0) {
  stop("No quantification GTF files were found - check --quant_dir and sample IDs in --metadata.",
       call. = FALSE)
}

# ---- 4. build the wide matrix: every candidate transcript_id x every sample found ----
matrix <- purrr::reduce(sample_frames, dplyr::full_join, by = "transcript_id") %>%
  dplyr::right_join(tibble(transcript_id = sort(candidate_ids)), by = "transcript_id") %>%
  dplyr::mutate(dplyr::across(-transcript_id, ~ tidyr::replace_na(.x, 0)))

sample_cols <- setdiff(names(matrix), "transcript_id")

# Explicit base indexing for pmax avoids dplyr::select S4 method collisions
matrix <- matrix %>%
  dplyr::mutate(
    mean_TPM           = rowMeans(matrix[, sample_cols, drop = FALSE]),
    max_TPM             = do.call(pmax, matrix[, sample_cols, drop = FALSE]),
    n_samples_detected = rowSums(matrix[, sample_cols, drop = FALSE] > 0)
  )

matrix_out <- paste0(opt$out_prefix, ".matrix.tsv")
write_tsv(matrix, matrix_out)
cat("Wrote per-sample TPM matrix:", matrix_out, "\n")

# ---- 5. optional per-tissue mean TPM, for joining into the final annotation table ----
if (!is.null(tissue_map)) {
  present_samples <- intersect(sample_cols, tissue_map$sample_id)
  tissues <- unique(tissue_map$tissue[tissue_map$sample_id %in% present_samples])

  tissue_means <- purrr::map_dfc(tissues, function(t) {
    cols <- tissue_map$sample_id[tissue_map$tissue == t & tissue_map$sample_id %in% present_samples]
    tibble(!!paste0("mean_TPM_", t) := rowMeans(matrix[, cols, drop = FALSE]))
  })

  summary <- bind_cols(
    dplyr::select(matrix, transcript_id, mean_TPM, max_TPM, n_samples_detected),
    tissue_means
  )

  summary_out <- paste0(opt$out_prefix, ".summary.tsv")
  write_tsv(summary, summary_out)
  cat("Wrote per-tissue TPM summary:", summary_out, "\n")
} else {
  cat("No tissue/condition column detected in metadata (needs >= 2 columns) - ",
      "skipping per-tissue summary; only the sample-level matrix was written.\n", sep = "")
}