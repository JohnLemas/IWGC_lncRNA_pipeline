#!/usr/bin/env Rscript
library(LncFinder)
library(seqinr)   # read input fasta sequences
library(parallel)  # run in parallel
library(e1071)     # for SVM model

# ==========================================
# SET PATHS & CONFIGURATION
# ==========================================
input_fasta   <- "/mnt/scratch/lemasjoh/BASF/IWGC_lncRNA_pipeline/species/03_outputs/04_gffread/species_candidates.fasta"
fold_output   <- "/mnt/scratch/lemasjoh/BASF/IWGC_lncRNA_pipeline/species/03_outputs/08_lncFinder/folded_structures.txt"
output_dir    <- "/mnt/scratch/lemasjoh/BASF/IWGC_lncRNA_pipeline/species/03_outputs/08_lncFinder"
output_csv    <- file.path(output_dir, "species_lncFinder_wheat_model_results.csv")
num_cores     <- 32

# ==========================================
# LOAD STRUCTURE DATA VIA LNCFINDER NATIVE PARSER
# ==========================================
cat(">>> [Step 2] Loading secondary structures using native read_SS parser...\n")

# Read the structures
SS_Raw <- read_SS(fold_output)

# FORCE conversion to a clean Data Frame to prevent the "incorrect number of dimensions" error
SS_List <- as.data.frame(SS_Raw, stringsAsFactors = FALSE)

# Log the dimensions for your slurm .out log to verify it has rows and columns
cat(sprintf(">>> Successfully loaded structured dataframe with %d rows and %d columns.\n",
            nrow(SS_List), ncol(SS_List)))

# ==========================================
# [Step 3] RUN LNCFINDER PREDICTIONS
# ==========================================
cat(">>> [Step 3] Computing LncFinder features and running SVM model...\n")

results <- lnc_finder(
  Sequences        = SS_List,     # Pass the explicitly formatted data frame
  SS.features      = TRUE,
  format           = "SS",
  frequencies.file = "wheat",
  svm.model        = "wheat",
  parallel.cores   = num_cores
)
# ==========================================
# SAVE RESULTS
# ==========================================
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

write.csv(results, output_csv, row.names = TRUE)
cat(">>> Success! Results saved to:\n\t", output_csv, "\n")