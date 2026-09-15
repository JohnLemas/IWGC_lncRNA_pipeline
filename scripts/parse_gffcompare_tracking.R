#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(stringr)
})

# 1. Define command-line options
option_list <- list(
  make_option(c("-t", "--tracking"), type = "character", default = NULL,
              help = "Path to the gffcompare candidate.tracking file [REQUIRED]", metavar = "file"),
  make_option(c("-l", "--list"), type = "character", default = NULL,
              help = "Path to the assembly file list used as gffcompare -i input [REQUIRED]", metavar = "file"),
  make_option(c("-f", "--fasta"), type = "character", default = NULL,
              help = "Path to FASTA of transcript IDs to track (full candidate set or a filtered subset) [REQUIRED]", metavar = "file"),
  make_option(c("-m", "--metadata"), type = "character", default = NULL,
              help = "Path to sample metadata (col1 = read ID, col2 = sample name) [REQUIRED]", metavar = "file"),
  make_option(c("-o", "--out"), type = "character", default = NULL,
              help = "Path for the output tracking summary TSV [REQUIRED]", metavar = "file")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# 2. Validate input arguments
required <- c("tracking", "list", "fasta", "metadata", "out")
missing <- required[vapply(required, function(x) is.null(opt[[x]]), logical(1))]
if (length(missing) > 0) {
  print_help(opt_parser)
  stop(paste0("Error: missing required argument(s): ",
              paste(paste0("--", missing), collapse = ", ")), call. = FALSE)
}

for (f in c(opt$tracking, opt$list, opt$fasta, opt$metadata)) {
  if (!file.exists(f)) {
    stop(paste0("Error: input file '", f, "' does not exist."), call. = FALSE)
  }
}

# 3. Tissue/stage code translation table
# (e.g. "sdl_rt_2" -> "Seedling Root", "4-6tls_lvs_1" -> "4-6 True Leaf Stage Leaves")
# Extend this if new codes appear in future metadata files.
TISSUE_CODE_MAP <- c(
  "sdl"    = "Seedling",
  "sdls"   = "Seedlings",
  "lvs"    = "Leaves",
  "rt"     = "Root",
  "rts"    = "Roots",
  "4-6tls" = "4-6 True Leaf Stage",
  "etio"   = "Etiolated",
  "etl"    = "Etiolated",
  "flbd"   = "Flower Bud",
  "fldb"   = "Flower Bud",
  "fls"    = "Flowers",
  "exp"    = "Expanded"
)

sample_name_to_tissue <- function(sample_name) {
  # Strip the trailing replicate number, then translate each remaining
  # underscore-delimited code into a human-readable label.
  stripped <- str_replace(sample_name, "_\\d+$", "")
  tokens <- str_split(stripped, "_", simplify = FALSE)[[1]]
  translated <- ifelse(tokens %in% names(TISSUE_CODE_MAP),
                       unname(TISSUE_CODE_MAP[tokens]),
                       tokens)
  paste(translated, collapse = " ")
}

# 4. Load sample metadata (col1 = read/sample ID, col2 = sample name)
meta <- read_tsv(opt$metadata, col_names = FALSE,
                 col_types = cols(.default = "c"), show_col_types = FALSE)
if (ncol(meta) < 2) {
  stop("Error: metadata file must have at least two tab-separated columns.", call. = FALSE)
}

meta <- meta %>%
  rename(Read_ID = X1, Sample_Name = X2) %>%
  filter(!is.na(Sample_Name), str_trim(Sample_Name) != "") %>%
  mutate(
    Read_ID     = str_trim(Read_ID),
    Sample_Name = str_trim(Sample_Name),
    Tissue      = vapply(Sample_Name, sample_name_to_tissue, character(1), USE.NAMES = FALSE)
  )

sample_to_read_id <- setNames(meta$Read_ID, meta$Sample_Name)
sample_to_tissue <- setNames(meta$Tissue, meta$Sample_Name)

cat("Loaded tissue/sample metadata for", nrow(meta), "samples\n")

# 5. Target transcript IDs from the FASTA
fasta_lines <- read_lines(opt$fasta)
headers <- fasta_lines[str_detect(fasta_lines, "^>")]

# Take the first whitespace-delimited token of each header, then if it is
# pipe-delimited keep the last field (e.g. "gene|MSTRG.1.1" -> "MSTRG.1.1").
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

# 6. Map tracking file column positions to sample identifiers
# Strip the "_assembled.gtf" suffix but keep the full sample name intact,
# e.g. "sdl_rt_2_assembled.gtf" -> "sdl_rt_2". Splitting on the first
# underscore would collapse distinct tissue/replicate samples together.
list_entries <- str_trim(read_lines(opt$list))
list_entries <- list_entries[list_entries != ""]
sample_names <- str_remove(basename(list_entries), "_assembled\\.(gtf|gff3?|gff)$")

cat("Loaded", length(sample_names), "sample identifiers from the assembly list\n")

# 7. Scan the tracking matrix
cat("Extracting per-transcript tracking info for target sequences...\n")

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
  stop("Error: no rows extracted. Verify that your FASTA IDs exactly match the tracking file IDs.",
       call. = FALSE)
}

# 8. Compile, sort (tissue breadth first, then sample count), and write
final_df <- bind_rows(rows) %>%
  arrange(desc(Tissue_Count), desc(Sample_Count))

write_tsv(final_df, opt$out, na = "NA")

cat("Tracking summary written to:", opt$out, "\n")
cat("Rows:", nrow(final_df), " Columns:", ncol(final_df), "\n")
cat("\nTop 10 candidates by tissue breadth:\n")
print(as.data.frame(head(final_df[, c("Predicted_Gene_ID", "Sample_Count", "Tissue_Count", "Tissue_Types")], 10)),
      row.names = FALSE)
