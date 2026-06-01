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

cat(">>> [Step 2] Parsing RNAfold structures into memory...\n")

# ==========================================
# PARSE STRUCTURE TEXT FILE
# ==========================================
lines <- readLines(fold_output)

# Identify lines mathematically by position (3 lines per entry in RNAfold output)
header_idx <- seq(1, length(lines), by = 3)
struct_idx <- seq(3, length(lines), by = 3)

# Extract and clean headers
headers <- gsub("^>", "", lines[header_idx])

# Extract dot-bracket notation, dropping the MFE tracking value at the end
raw_structs <- lines[struct_idx]
structures  <- sapply(strsplit(raw_structs, " "), `[`, 1)
names(structures) <- headers

# ==========================================
# LOAD FASTA & ALIGN
# ==========================================
Seqs <- read.fasta(input_fasta, seqtype = "DNA")
fasta_names <- names(Seqs)

# Identify which sequences were successfully folded
folded_names <- names(structures)

# Count how many were skipped for your logs
skipped_count <- sum(!fasta_names %in% folded_names)
if (skipped_count > 0) {
  cat(sprintf(">>> WARNING: %d ultra-long sequences were skipped by RNAfold. Filtering them out...\n", skipped_count))
  
  # Keep only the fasta entries that match what successfully folded
  Seqs <- Seqs[fasta_names %in% folded_names]
  fasta_names <- names(Seqs)
}

# Reorder structure vector to perfectly match the filtered FASTA arrangement
Structures <- structures[fasta_names]

# ==========================================
# FORMAT FOR LNCFINDER SS MODULE
# ==========================================
# LncFinder's internal 'read_SS' logic expects a list structure where each element 
# contains the secondary structure string vector (Dot-Bracket Notation).
# We convert our named vector into the list format LncFinder expects.

SS_List <- lapply(Structures, function(x) {
  # Split the dot-bracket string into individual character elements
  strsplit(x, "")[[1]]
})
names(SS_List) <- fasta_names

# ==========================================
# [Step 3] RUN LNCFINDER PREDICTIONS
# ==========================================
cat(">>> [Step 3] Computing LncFinder structural features and running SVM model...\n")

results <- lnc_finder(
  Sequences        = SS_List,     # Pass the structure list directly here!
  SS.features      = TRUE,        # Keep secondary structure features active
  format           = "SS",        # CHANGED: Must be "SS" when feeding structures
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