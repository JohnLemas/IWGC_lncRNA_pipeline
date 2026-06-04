#!/usr/bin/env Rscript

# Load required libraries
if (!requireNamespace("optparse", quietly = TRUE)) {
  stop("CRITICAL ERROR: The 'optparse' package is required but not installed in this container.", call. = FALSE)
}
library(optparse)
library(LncFinder)

# ==========================================
# DEFINE EXPECTED COMMAND-LINE FLAGS
# ==========================================
option_list <- list(
  make_option(c("-s", "--struct"), type = "character", default = NULL,
              help = "Path to the raw input folded_structures.txt file from RNAfold", 
              metavar = "character"),
  make_option(c("-o", "--out"), type = "character", default = NULL,
              help = "Path where the output results CSV should be saved", 
              metavar = "character"),
  make_option(c("-p", "--cores"), type = "integer", default = 1,
              help = "Number of parallel CPU cores to use [default %default]", 
              metavar = "integer")
)

# Parse the incoming arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Validate mandatory flags
if (is.null(opt$struct) || is.null(opt$out)) {
  print_help(opt_parser)
  stop("CRITICAL ERROR: Both flags -s/--struct and -o/--out must be provided.\n", call. = FALSE)
}

cat("==========================================================\n")
cat("Running LncFinder with Named Flags:\n")
cat(sprintf("  Input Structure File (-s): %s\n", opt$struct))
cat(sprintf("  Output Results CSV   (-o): %s\n", opt$out))
cat(sprintf("  Parallel Cores       (-p): %d\n", opt$cores))
cat("==========================================================\n")

# ==========================================
# LOAD STRUCTURE DATA VIA LNCFINDER NATIVE PARSER
# ==========================================
cat(">>> [Step 2] Loading secondary structures using native read_SS parser...\n")

if (!file.exists(opt$struct)) {
  stop(sprintf("CRITICAL ERROR: The structure file '%s' does not exist!", opt$struct), call. = FALSE)
}

# Parse the single combined file 
SS_Raw <- read_SS(oneFile.loc = opt$struct, separateFile = FALSE, withMFE = TRUE)

# Convert to a clean Data Frame to prevent dimensions errors
SS_List <- as.data.frame(SS_Raw, stringsAsFactors = FALSE)

cat(sprintf(">>> Successfully loaded structured dataframe with %d rows and %d columns.\n", 
            nrow(SS_List), ncol(SS_List)))

# ==========================================
# [Step 3] RUN LNCFINDER PREDICTIONS
# ==========================================
cat(">>> [Step 3] Computing LncFinder structural features and running SVM model...\n")

results <- lnc_finder(
  Sequences        = SS_List,     
  SS.features      = TRUE,        
  format           = "SS",        
  frequencies.file = "wheat",
  svm.model        = "wheat",
  parallel.cores   = opt$cores
)

# ==========================================
# SAVE RESULTS
# ==========================================
output_dir <- dirname(opt$out)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

write.csv(results, opt$out, row.names = TRUE)
cat(">>> Success! Results cleanly saved to:\n\t", opt$out, "\n")