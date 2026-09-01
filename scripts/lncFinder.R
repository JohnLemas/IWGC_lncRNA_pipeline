#!/usr/bin/env Rscript

library(optparse)
library(seqinr)
library(LncFinder)

# ==========================================
# DEFINE EXPECTED COMMAND-LINE FLAGS
# ==========================================
option_list <- list(
  make_option(c("-i", "--fasta"), type = "character", default = NULL,
              help = "Path to the raw input candidate FASTA file", metavar = "character"),
  make_option(c("-o", "--out"), type = "character", default = NULL,
              help = "Path where the output results CSV should be saved", metavar = "character"),
  make_option(c("-p", "--cores"), type = "integer", default = 1,
              help = "Number of parallel CPU cores to use [default %default]", metavar = "integer"),
  make_option(c("-m", "--model"), type = "character", default = "wheat",
              help = "Model type: 'wheat' (monocot), 'dicot', 'human', 'mouse', or 'other' [default %default]", metavar = "character")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$fasta) || is.null(opt$out)) {
  stop("CRITICAL ERROR: Both flags -i/--fasta and -o/--out must be provided.\n", call. = FALSE)
}

# Map user-friendly model flags to LncFinder internal parameters
model_type <- tolower(opt$model)

if (model_type == "wheat") {
  freq_file <- "wheat"
  svm_mod   <- "wheat"
} else if (model_type == "dicot") {
  freq_file <- "A.thaliana"
  svm_mod   <- "A.thaliana"
} else if (model_type == "human") {
  freq_file <- "human"
  svm_mod   <- "human"
} else if (model_type == "mouse") {
  freq_file <- "mouse"
  svm_mod   <- "mouse"
} else if (model_type == "other") {
  freq_file <- "A.thaliana"
  svm_mod   <- "human"
} else {
  stop(sprintf("CRITICAL ERROR: Invalid model '%s'. Choose from: wheat, dicot, human, mouse, other.\n", opt$model), call. = FALSE)
}

# Ensure output directory exists
output_dir <- dirname(opt$out)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ==========================================
# STEP 1: LOAD & FILTER FASTA IN R
# ==========================================
cat(">>> [Step 1] Loading and filtering FASTA dataset in memory...\n")
raw_fasta <- read.fasta(opt$fasta, seqtype = "DNA", as.string = TRUE)

seq_lengths <- sapply(raw_fasta, nchar)
filtered_fasta <- raw_fasta[seq_lengths >= 200 & seq_lengths <= 5000]

cat(sprintf(">>> Filtered out %d outlier sequences. Processing %d valid candidates.\n", 
            length(raw_fasta) - length(filtered_fasta), length(filtered_fasta)))

filtered_fasta_path <- file.path(output_dir, "filtered_candidates.fasta")
write.fasta(
  sequences = filtered_fasta, 
  names     = names(filtered_fasta), 
  file.out  = filtered_fasta_path
)

# ==========================================
# STEP 2: RUN NATIVE RNAFOLD VIA LNCFINDER
# ==========================================
cat(">>> [Step 2] Calling run_RNAfold() natively on filtered sequences...\n")

SS_Raw <- run_RNAfold(
  Sequences      = filtered_fasta,
  RNAfold.path   = "RNAfold",
  parallel.cores = opt$cores
)

SS_List <- as.data.frame(SS_Raw, stringsAsFactors = FALSE)
cat(sprintf(">>> Successfully generated structured matrix with %d sequences.\n", nrow(SS_List)))

# ==========================================
# STEP 3: RUN LNCFINDER SVM PREDICTION
# ==========================================
cat(sprintf(">>> [Step 3] Running LncFinder SVM model using selected model profile: '%s'...\n", model_type))

results <- lnc_finder(
  Sequences        = SS_List,     
  SS.features      = TRUE,        
  format           = "SS",        
  frequencies.file = freq_file,
  svm.model        = svm_mod,
  parallel.cores   = opt$cores
)

# ==========================================
# SAVE RESULTS
# ==========================================
# Convert row names (sequence IDs) into an explicit column
results$Sequence_ID <- rownames(results)

# Reorder columns so Sequence_ID is first
results <- results[, c("Sequence_ID", setdiff(colnames(results), "Sequence_ID"))]

# Save without R's automatic row index column
write.csv(results, opt$out, row.names = FALSE)
cat(">>> Success! Results cleanly saved to:\n\t", opt$out, "\n")