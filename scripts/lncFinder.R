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
              help = "Model type: 'wheat' (only plant-trained model LncFinder ships), 'human', or 'mouse' [default %default]",
              metavar = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$fasta) || is.null(opt$out)) {
  stop("CRITICAL ERROR: Both flags -i/--fasta and -o/--out must be provided.\n", call. = FALSE)
}

# ==========================================
# VALIDATE MODEL CHOICE
# ==========================================
# LncFinder only ships three pre-built models: human, mouse, wheat.
# There is no "dicot", "A.thaliana", or "other" model in this package.
# For dicots (e.g. soybean / g.max), "wheat" is the closest available option
# (it's the only plant-trained model), but it was trained on a monocot, so
# expect reduced accuracy relative to a model trained on dicot sequences.
# For real accuracy on a dicot species, train a custom model with
# make_frequencies() + build_model() using a dicot reference transcriptome.
model_type <- tolower(opt$model)
valid_models <- c("wheat", "human", "mouse")
if (!(model_type %in% valid_models)) {
  stop(sprintf(paste0(
    "CRITICAL ERROR: Invalid model '%s'. LncFinder only ships pre-built models for: %s.\n",
    "  There is no dicot/A.thaliana/other model in this package.\n",
    "  For dicot species, use 'wheat' as the closest available option, or train\n",
    "  a custom model with make_frequencies()/build_model() on your own reference data.\n"
  ), opt$model, paste(valid_models, collapse = ", ")), call. = FALSE)
}
freq_file <- model_type
svm_mod   <- model_type

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

if (length(filtered_fasta) == 0) {
  stop("CRITICAL ERROR: No sequences survived the length filter (200-5000 nt).\n", call. = FALSE)
}

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

# IMPORTANT: do NOT wrap this in as.data.frame(). run_RNAfold() returns a
# list structured for direct use by lnc_finder(format = "SS") - each element
# is a per-sequence vector (raw seq, dot-bracket structure, MFE, etc).
# Coercing that list to a data.frame transposes it: sequences become columns
# and the handful of per-sequence fields become rows, so nrow() ends up
# reporting the field count (e.g. 3) instead of the sequence count.
cat(sprintf(">>> Successfully folded %d sequences with RNAfold.\n", length(SS_Raw)))

# ==========================================
# STEP 3: RUN LNCFINDER SVM PREDICTION
# ==========================================
cat(sprintf(">>> [Step 3] Running LncFinder SVM model using selected model profile: '%s'...\n", model_type))

results <- lnc_finder(
  Sequences        = SS_Raw,
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