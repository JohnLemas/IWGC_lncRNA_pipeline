#!/usr/bin/env Rscript
# parse_tpm_matrix.R
#
# Aggregates per-sample TPM values (read from each sample's StringTie -e quantification GTF)
# into a transcript x sample TPM matrix, restricted to candidate lncRNAs.

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
              help = "Directory containing <sample_id>/<sample_id>_quant.gtf subfolders (from stringtie -e)"),
  make_option(c("-c", "--candidates_gtf"), type = "character",
              help = "species_candidates.gtf from stringtie_merge.sb (defines candidate transcript_ids)"),
  make_option(c("-m", "--metadata"), type = "character",
              help = "species_metadata.txt: col 1 = orig ID, col 2 = sample_id, col 3 = tissue_type"),
  make_option(c("-o", "--out_prefix"), type = "character",
              help = "Prefix for output files: <prefix>.matrix.tsv and <prefix>.summary.tsv"),
  make_option("--allow_missing", action = "store_true", default = TRUE,
              help = "Warn instead of stopping when candidate transcript_ids are absent from quantification GTFs; zero-fills absent candidates.")
)
opt <- parse_args(OptionParser(option_list = option_list))

for (req in c("quant_dir", "candidates_gtf", "metadata", "out_prefix")) {
  if (is.null(opt[[req]])) stop("Missing required argument: --", req, call. = FALSE)
}

# 1. Transcript IDs to keep
gtf <- import(opt$candidates_gtf)
candidate_ids <- unique(na.omit(mcols(gtf)$transcript_id))
if (length(candidate_ids) == 0) {
  stop("No transcript_id values found in ", opt$candidates_gtf, call. = FALSE)
}
cat("Loaded", length(candidate_ids), "candidate lncRNA transcript IDs from", basename(opt$candidates_gtf), "\n")

# 2. Sample IDs + Tissue mapping
meta_lines <- read_lines(opt$metadata)
meta_lines <- meta_lines[nzchar(trimws(meta_lines))]
meta_split <- strsplit(meta_lines, "\t")

if (!all(lengths(meta_split) >= 3)) {
  stop("Expected at least 3 tab-separated columns in --metadata.", call. = FALSE)
}

sample_ids <- trimws(vapply(meta_split, `[`, character(1), 2))
tissue_map <- tibble(
  sample_id = sample_ids,
  tissue = trimws(vapply(meta_split, `[`, character(1), 3))
)

# 3. Read sample GTFs
read_sample_tpm <- function(sample_id) {
  gtf_path <- file.path(opt$quant_dir, sample_id, paste0(sample_id, "_quant.gtf"))
  if (!file.exists(gtf_path)) {
    warning("No quantification GTF found for ", sample_id, ", skipping.", call. = FALSE)
    return(NULL)
  }
  gr <- import(gtf_path)
  if (!"TPM" %in% names(mcols(gr))) {
    warning("Quantification GTF for ", sample_id, " has no TPM attribute, skipping.", call. = FALSE)
    return(NULL)
  }
  as_tibble(as.data.frame(mcols(gr))) %>%
    filter(type == "transcript") %>%
    mutate(TPM = as.numeric(TPM)) %>%
    transmute(transcript_id = transcript_id, !!sample_id := TPM) %>%
    filter(transcript_id %in% candidate_ids) %>%
    distinct(transcript_id, .keep_all = TRUE)
}

sample_frames <- purrr::compact(purrr::map(sample_ids, read_sample_tpm))
if (length(sample_frames) == 0) {
  stop("No quantification GTFs were found. Check --quant_dir and --metadata.", call. = FALSE)
}

# 4. Build wide TPM matrix
wide <- purrr::reduce(sample_frames, dplyr::full_join, by = "transcript_id") %>%
  mutate(across(-transcript_id, ~ replace_na(.x, 0)))

missing_ids <- setdiff(candidate_ids, wide$transcript_id)
if (length(missing_ids) > 0) {
  msg <- sprintf("%d of %d candidate transcript_ids (%.1f%%) are absent from quantification GTFs.",
                 length(missing_ids), length(candidate_ids),
                 100 * length(missing_ids) / length(candidate_ids))
  if (isTRUE(opt$allow_missing)) {
    cat("Warning:", msg, "Zero-filling unexpressed candidates.\n")
  } else {
    stop(msg, " Re-run with --allow_missing to override.", call. = FALSE)
  }
}

matrix <- wide %>%
  right_join(tibble(transcript_id = sort(candidate_ids)), by = "transcript_id") %>%
  mutate(across(-transcript_id, ~ replace_na(.x, 0)))

sample_cols <- setdiff(names(matrix), "transcript_id")

matrix <- matrix %>%
  mutate(
    mean_TPM           = rowMeans(across(all_of(sample_cols))),
    max_TPM            = do.call(pmax, as.list(select(matrix, all_of(sample_cols)))),
    n_samples_detected = rowSums(across(all_of(sample_cols), ~ .x > 0))
  )

matrix_out <- paste0(opt$out_prefix, ".matrix.tsv")
write_tsv(matrix, matrix_out)
cat("Wrote per-sample TPM matrix:", matrix_out, "\n")

# 5. Per-tissue mean TPM
present_samples <- intersect(sample_cols, tissue_map$sample_id)
tissues <- unique(tissue_map$tissue[tissue_map$sample_id %in% present_samples])

tissue_means <- purrr::map_dfc(tissues, function(t) {
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