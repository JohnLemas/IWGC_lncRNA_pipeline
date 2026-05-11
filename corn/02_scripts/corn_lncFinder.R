#!/usr/bin/env Rscript
library(LncFinder)
library(seqinr)
library(parallel)
library(e1071)

# ==========================================
# SET PATHS
# ==========================================
input_fasta   <- "/mnt/scratch/lemasjoh/BASF/RNA_central/corn/03_output/04_gffread/corn_candidates.fasta"
output_dir    <- "/mnt/scratch/lemasjoh/BASF/RNA_central/corn/03_output/08_lncFinder"
output_csv    <- file.path(output_dir, "corn_lncFinder_wheat_model_results.csv")

cat(">>>LncFinder:\ninput data:\n\t", input_fasta,
    "\n output file:\n\t", output_csv,
    "\nRunning predictions with built-in wheat model..\n"
    )

# ==========================================
# READ INPUT SEQUENCES
# ==========================================

Seqs <- read.fasta(input_fasta, seqtype = "DNA")

# ==========================================
# RUN COMMAND
# ==========================================

results <- lnc_finder(
  Seqs, 
  SS.features = FALSE, 
  format = "DNA", 
  frequencies.file = "wheat", 
  svm.model = "wheat", 
  parallel.cores = 32
)

# ==========================================
# SAVE OUTPUT
# ==========================================

write.csv(results, output_csv, row.names = TRUE)
cat(">>>LncFinder:\n\tDone!\n\tResults saved to", output_csv, "\n")