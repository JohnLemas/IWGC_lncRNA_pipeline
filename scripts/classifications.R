#!/usr/bin/env Rscript

library(readr)
library(dplyr)

input_file <- "lncRNA_classes.txt"

if (!file.exists(input_file)) {
    stop(paste("Error: Input file", input_file, "not found in current working directory."))
}

message("Reading ", input_file, "...")

# Read the file without headers. 
# read_tsv automatically names columns X1, X2, X3... which maps perfectly to $1, $2, $3...
df <- read_tsv(input_file, col_names = FALSE, show_col_types = FALSE)

# Function to write outputs uniformly without headers or quotes, just like awk
write_output <- function(data_frame, file_name) {
    write_delim(data_frame, file_name, delim = "\t", col_names = FALSE)
}

# 1. Antisense_exonic-lncRNA
message("Processing: Antisense exonic...")
df %>% 
    filter(X1 == 1, X6 == "antisense", X10 == "exonic") %>% 
    write_output("LncRNA_antisense_exonic.txt")

# 2. Intronic-lncRNA
message("Processing: Intronic...")
df %>% 
    filter(X1 == 1, X6 == "sense", X10 == "intronic") %>% 
    write_output("LncRNA_intronic.txt")

# 3. Upstream-lncRNA
message("Processing: Upstream...")
df %>% 
    filter(X1 == 1, X6 == "sense", X7 == "intergenic", X8 <= 2000, X10 == "upstream") %>% 
    write_output("LncRNA_upstream.txt")

# 4. Downstream-lncRNA
message("Processing: Downstream...")
df %>% 
    filter(X1 == 1, X6 == "sense", X7 == "intergenic", X8 <= 2000, X10 == "downstream") %>% 
    write_output("LncRNA_downstream.txt")

# 5. Intergenic-lncRNA
message("Processing: Intergenic...")
df %>% 
    filter(X1 == 1, X7 == "intergenic", X8 > 2000) %>% 
    write_output("LncRNA_intergenic.txt")

# 6. Bidirectional-lncRNA
message("Processing: Bidirectional...")
df %>% 
    filter(X1 == 1, X6 == "antisense", X7 == "intergenic", X8 <= 2000, X10 == "upstream") %>% 
    write_output("LncRNA_Bidirectional.txt")

message("All classification files generated successfully!")