#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(stringr)
})

option_list <- list(
  make_option(c("-t", "--tracking"), type = "character", default = NULL,
              help = "Path to the gffcompare candidate.tracking file [REQUIRED]", metavar = "file"),
  make_option(c("-l", "--list"), type = "character", default = NULL,
              help = "Path to the assembly file list used as gffcompare -i input [REQUIRED]", metavar = "file"),
  make_option(c("-f", "--fasta"), type = "character", default = NULL,
              help = "Path to FASTA of transcript IDs to track [REQUIRED]", metavar = "file"),
  make_option(c("-m", "--metadata"), type = "character", default = NULL,
              help = "Path to sample metadata (col1 = read ID, col2 = sample name, col3 = tissue) [REQUIRED]", metavar = "file"),
  make_option(c("-o", "--out"), type = "character", default = NULL,
              help = "Path for output tracking summary TSV [REQUIRED]", metavar = "file")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

required <- c("tracking", "list", "fasta", "metadata", "out")
missing <- required[vapply(required, function(x) is.null(opt[[x]]), logical(1))]
if (length(missing) > 0) {
  print_help(opt_parser)
  stop(paste0("Error: missing required argument(s): ", paste(paste0("--", missing), collapse = ", ")), call. = FALSE)
}

for (f in c(opt$tracking, opt$list, opt$fasta, opt$metadata)) {
  if (!file.exists(f)) stop(paste0("Error: input file '", f, "' does not exist."), call. = FALSE)
}

meta <- read_tsv(opt$metadata, col_names = FALSE, col_types = cols(.default = "c"), show_col_types = FALSE)
if (ncol(meta) < 3) {
  stop("Error: metadata file must have at least three tab-separated columns.", call. = FALSE)
}

meta <- meta %>%
  select(Read_ID = 1, Sample_Name = 2, Tissue = 3) %>%
  mutate(
    Read_ID     = str_trim(Read_ID),
    Sample_Name = str_trim(Sample_Name),
    Tissue      = na_if(str_trim(Tissue), "")
  ) %>%
  filter(!is.na(Sample_Name), Sample_Name != "")

sample_to_read_id <- setNames(meta$Read_ID, meta$Sample_Name)
sample_to_tissue <- setNames(meta$Tissue, meta$Sample_Name)

cat("Loaded tissue/sample metadata for", nrow(meta), "samples\n")

fasta_lines <- read_lines(opt$fasta)
headers <- fasta_lines[str_detect(fasta_lines, "^>")]
first_tokens <- str_remove(str_trim(headers), "^>")
first_tokens <- vapply(str_split(first_tokens, "\\s+"), function(x) x[1], character(1))

validated_ids <- unique(vapply(first_tokens, function(tok) {
  if (str_detect(tok, fixed("|"))) {
    parts <- str_split(tok, fixed("|"))[[1]]
    parts[length(parts)]
  } else {
    tok
  }
}, character(1), USE.NAMES = FALSE))

cat("Loaded", length(validated_ids), "transcript IDs to target\n")

list_entries <- str_trim(read_lines(opt$list))
list_entries <- list_entries[list_entries != ""]
sample_names <- str_remove(basename(list_entries), "_assembled\\.(gtf|gff3?|gff)$")

cat("Loaded", length(sample_names), "sample identifiers from the assembly list\n")

track_lines <- read_lines(opt$tracking)
parts_list <- str_split(track_lines, "\t")
parts_list <- parts_list[vapply(parts_list, length, integer(1)) >= 5]

validated_env <- new.env(hash = TRUE, parent = emptyenv(), size = length(validated_ids) * 2L)
for (vid in validated_ids) assign(vid, TRUE, envir = validated_env)
n_samples <- length(sample_names)

rows <- lapply(parts_list, function(parts) {
  tcons_raw <- str_split(parts[1], fixed("|"), simplify = FALSE)[[1]][1]
  ref_raw <- parts[3]
  sample_columns <- parts[5:length(parts)]

  if (!is.na(ref_raw) && ref_raw != "-" && str_trim(ref_raw) != "") {
    ref_parts <- str_split(ref_raw, fixed("|"), simplify = FALSE)[[1]]
    candidate_id <- if (length(ref_parts) > 1) ref_parts[2] else ref_parts[1]
  } else {
    candidate_id <- tcons_raw
  }

  if (is.na(candidate_id) || !nzchar(candidate_id)) return(NULL)
  if (!exists(candidate_id, envir = validated_env, inherits = FALSE)) return(NULL)

  n <- min(length(sample_columns), n_samples)
  if (n == 0) return(NULL)

  present <- which(str_trim(sample_columns[seq_len(n)]) != "-")
  if (length(present) == 0) return(NULL)

  samples_present <- sample_names[present]
  read_ids <- sort(unique(unname(sample_to_read_id[samples_present])))
  tissues <- sort(unique(unname(sample_to_tissue[samples_present])))
  read_ids <- read_ids[!is.na(read_ids)]
  tissues <- tissues[!is.na(tissues)]

  tibble(
    Predicted_Gene_ID = candidate_id,
    Sample_Count = length(samples_present),
    Sample_IDs = paste(sort(samples_present), collapse = ","),
    Read_IDs = paste(read_ids, collapse = ","),
    Tissue_Count = length(tissues),
    Tissue_Types = paste(tissues, collapse = ",")
  )
})

rows <- rows[!vapply(rows, is.null, logical(1))]

if (length(rows) == 0) {
  stop("Error: no rows extracted. Verify FASTA IDs match tracking file IDs.", call. = FALSE)
}

final_df <- bind_rows(rows) %>%
  arrange(desc(Tissue_Count), desc(Sample_Count))

write_tsv(final_df, opt$out, na = "NA")
cat("Tracking summary written to:", opt$out, "\n")