#!/usr/bin/env Rscript

library(optparse)
library(readr)
library(dplyr)

# 1. Define command-line options
option_list <- list(
    make_option(c("-i", "--input"), type = "character", default = NULL, 
                help = "Path to the FEELnc classification result file (e.g., lncRNA_classes.txt)", metavar = "character"),
    make_option(c("-o", "--outdir"), type = "character", default = ".", 
                help = "Directory where output files should be saved [default: current directory]", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# 2. Validate input argument
if (is.null(opt$input)) {
    print_help(opt_parser)
    stop("Error: Input file must be specified using -i or --input.", call. = FALSE)
}

if (!file.exists(opt$input)) {
    stop(paste("Error: Input file '", opt$input, "' does not exist."), call. = FALSE)
}

# 3. Handle output directory creation if it doesn't exist
if (!dir.exists(opt$outdir)) {
    dir.create(opt$outdir, recursive = TRUE)
}

message("Reading ", opt$input, "...")

# Read the file without headers (X1, X2... map to $1, $2...)
df <- read_tsv(opt$input, col_names = FALSE, show_col_types = FALSE)

# Helper function to join output directory path and write file
write_output <- function(data_frame, base_filename) {
    target_path <- file.path(opt$outdir, base_filename)
    write_delim(data_frame, target_path, delim = "\t", col_names = FALSE)
}

# 4. Filter and write out files using dplyr logic
message("Processing data and writing to: ", opt$outdir)

# Antisense_exonic-lncRNA
df %>% 
    filter(X1 == 1, X6 == "antisense", X10 == "exonic") %>% 
    write_output("LncRNA_antisense_exonic.txt")

# Intronic-lncRNA
df %>% 
    filter(X1 == 1, X6 == "sense", X10 == "intronic") %>% 
    write_output("LncRNA_intronic.txt")

# Upstream-lncRNA
df %>% 
    filter(X1 == 1, X6 == "sense", X7 == "intergenic", X8 <= 2000, X10 == "upstream") %>% 
    write_output("LncRNA_upstream.txt")

# Downstream-lncRNA
df %>% 
    filter(X1 == 1, X6 == "sense", X7 == "intergenic", X8 <= 2000, X10 == "downstream") %>% 
    write_output("LncRNA_downstream.txt")

# Intergenic-lncRNA
df %>% 
    filter(X1 == 1, X7 == "intergenic", X8 > 2000) %>% 
    write_output("LncRNA_intergenic.txt")

# Bidirectional-lncRNA
df %>% 
    filter(X1 == 1, X6 == "antisense", X7 == "intergenic", X8 <= 2000, X10 == "upstream") %>% 
    write_output("LncRNA_Bidirectional.txt")

message("All classification files generated successfully!")