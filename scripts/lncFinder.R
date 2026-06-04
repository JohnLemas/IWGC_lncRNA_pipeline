#!/usr/bin/env Rscript
library(optparse)
library(seqinr)
library(LncFinder)

# ==========================================
# DEFINE EXPECTED COMMAND-LINE FLAGS
# ==========================================
option_list <- list(
  make_option(c("-i", "--fasta"), type = "character", default = NULL,
              help = "Path to the full candidate FASTA file", metavar = "character"),
  make_option(c("-o", "--out"), type = "character", default = NULL,
              help = "Path where the output results CSV should be saved", metavar = "character"),
  make_option(c("-p", "--cores"), type = "integer", default = 1,
              help = "Number of parallel CPU cores to use [default %default]", metavar = "integer")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$fasta) || is.null(opt$out)) {
  stop("CRITICAL ERROR: Both flags -i/--fasta and -o/--out must be provided.\n", call. = FALSE)
}

# ==========================================
# STEP 1: LOAD & FILTER FASTA IN R
# ==========================================
cat(">>> [Step 1] Loading and filtering FASTA dataset natively in memory...\n")
raw_fasta <- read.fasta(opt$fasta, seqtype = "DNA", as.string = TRUE)

# Filter out sequences longer than 50,000 nt to prevent O(n^3) computational explosion
seq_lengths <- sapply(raw_fasta, nchar)
filtered_fasta <- raw_fasta[seq_lengths >= 200 & seq_lengths <= 50000]

cat(sprintf(">>> Filtered out %d outlier sequences. Processing %d valid candidates.\n",
            length(raw_fasta) - length(filtered_fasta), length(filtered_fasta)))

# ==========================================
# STEP 2 & 3: INTERNAL FOLDING & PREDICTION
# ==========================================
cat(">>> [Step 2] Computing secondary structures and running SVM model natively...\n")

# By passing format = "DNA" and SS.features = TRUE, LncFinder calls its internal
# RNAfold algorithm directly. This ensures perfect matrix alignment.
results <- lnc_finder(
  Sequences        = filtered_fasta,
  SS.features      = TRUE,
  format           = "DNA",
  frequencies.file = "wheat",
  svm.model        = "wheat",
  parallel.cores   = opt$cores
)

# ==========================================
# SAVE RESULTS
# ==========================================
output_dir <- dirname(opt$out)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

write.csv(results, opt$out, row.names = TRUE)
cat(">>> Success! Results cleanly saved to:\n\t", opt$out, "\n")