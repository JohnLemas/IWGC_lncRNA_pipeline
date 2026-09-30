#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tidyverse)
  library(optparse)
})

# This custom script accepts outputs from multiple lncRNA prediction tools, extracts candidate IDs, and identifies the consensus set of candidates that pass all filters. It also generates a summary table indicating which candidates passed each tool's criteria.
# This was designed to run with the intersect_ids.sif image.

# 1. DEFINE ARGUMENTS: arguments for each input and output used in the script.
option_list <- list(
  make_option(c("-f", "--feelnc"), type="character", help="FEELnc GTF or TXT list"),
  make_option(c("-c", "--cpat"), type="character", help="CPAT output (TSV)"),
  make_option(c("-b", "--boost"), type="character", help="lncBoost output (TSV)"),
  make_option(c("-l", "--lncfinder"), type="character", help="LncFinder output (CSV)"),
  make_option(c("-d", "--diamond"), type="character", help="DIAMOND output (TSV)"),
  make_option(c("-o", "--out"), type="character", default="consensus", help="Output prefix")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Case-insensitive matching: different tools can report the same transcript
# with different ID casing (e.g. "rna-XR_..." vs "RNA-XR_..."), which used to
# make reduce(provided, intersect) and setdiff() silently treat them as
# different transcripts. norm() is used ONLY to decide which IDs are "the
# same" - we still write out the original-cased ID everywhere, so exact-match
# downstream steps (grep -Fwf, seqkit grep, build_final_annotation.R's joins)
# keep working unchanged. If those still miss matches, the GTF itself is using
# yet a third casing convention and those steps need case-insensitive matching too.
norm <- function(x) tolower(trimws(as.character(x)))

# 2. DATA LOADING LOGIC (With Memory Optimization)
get_ids <- function(file, tool_name) {
  if (is.null(file) || !file.exists(file)) {
    cat(tool_name, ": [Skipped]\n")
    return(NULL)
  }
  cat("Loading", tool_name, "... ")
  
  if (tool_name == "FEELnc") {
    if (grepl("\\.gtf$", file, ignore.case = TRUE)) {
      # GTF Parser: Grabs transcript_id "MSTRG.1.1" -> MSTRG.1.1
      lines <- read_lines(file)
      ids <- str_extract(lines, 'transcript_id "[^"]+"') %>%
             str_replace_all('transcript_id |"', '') %>%
             unique()
      ids <- ids[!is.na(ids)]
    } else {
      ids <- read_lines(file)
    }
    
  } else if (tool_name == "CPAT") {
    # Fix header shift: skip header, name columns manually
    df <- read_delim(file, delim = "\t", skip = 1, col_names = FALSE, show_col_types = FALSE)
    ids <- df %>% filter(X6 < 0.46) %>% pull(X1)
    rm(df)
    
  } else if (tool_name == "LncFinder") {
    # Handles "NonCoding" filter and unnamed ID columns
    df <- read_csv(file, show_col_types = FALSE)
    ids <- df %>% filter(Pred == "NonCoding") %>% pull(1)
    rm(df)
    
  } else if (tool_name == "lncBoost") {
    df <- read_delim(file, delim = "\t", show_col_types = FALSE)
    ids <- df %>% filter(Predicted_label == 1) %>% pull(1)
    rm(df)
    
  } else if (tool_name == "DIAMOND") {
    # Negative filter: pulls unique IDs from first column
    df <- read_delim(file, delim = "\t", col_names = FALSE, show_col_types = FALSE)
    ids <- df %>% filter(X11 < 1e-5) %>% pull(1) %>% unique()
    rm(df)
  }
  
  ids <- trimws(as.character(ids))
  cat(length(ids), "candidates identified.\n")
  gc() # Force garbage collection to free RAM
  return(ids)
}

# 3. RUN PIPELINE
feelnc_ids <- get_ids(opt$feelnc, "FEELnc")
cpat_ids   <- get_ids(opt$cpat, "CPAT")
lncf_ids   <- get_ids(opt$lncfinder, "LncFinder")
boost_ids  <- get_ids(opt$boost, "lncBoost")
prot_hits  <- get_ids(opt$diamond, "DIAMOND")

# INTERSECTION
pos_list <- list(FEELnc=feelnc_ids, CPAT=cpat_ids, LncFinder=lncf_ids, Boost=boost_ids)
provided <- pos_list[!sapply(pos_list, is.null)]

if(length(provided) < 1) stop("No input files found!", call.=FALSE)

# Lookup from normalized ID -> first-seen original-cased ID, built across
# every tool's output plus DIAMOND, so normalized intersections/setdiffs can
# be mapped back to a real, on-disk ID string for the files we write out.
canonical_lookup <- c(provided, list(DIAMOND = prot_hits)) %>%
  compact() %>%
  unlist(use.names = FALSE) %>%
  { tibble(key = norm(.), id = .) } %>%
  distinct(key, .keep_all = TRUE) %>%
  deframe()

provided_norm <- map(provided, norm)
final_ids_norm <- reduce(provided_norm, intersect)

# Apply DIAMOND filter
if(!is.null(prot_hits)) {
  final_ids_norm <- setdiff(final_ids_norm, norm(prot_hits))
}

final_ids <- unname(canonical_lookup[final_ids_norm])

# 4. SAVE OUTPUTS
# Save the clean list first (safety in case of OOM on summary)
write_lines(final_ids, paste0(opt$out, "_final_ids.txt"))
cat("\n--- FINAL CONSENSUS COUNT:", length(final_ids), "---\n")

# Build and save Summary Table
cat("Building summary table...\n")
all_ids_norm <- unique(unlist(provided_norm, use.names = FALSE))
all_ids <- unname(canonical_lookup[all_ids_norm])
summary_tab <- data.frame(ID = all_ids, stringsAsFactors = FALSE)

if(!is.null(feelnc_ids)) summary_tab$Pass_FEELnc  <- all_ids_norm %in% norm(feelnc_ids)
if(!is.null(cpat_ids))   summary_tab$Pass_CPAT    <- all_ids_norm %in% norm(cpat_ids)
if(!is.null(lncf_ids))   summary_tab$Pass_LncF    <- all_ids_norm %in% norm(lncf_ids)
if(!is.null(boost_ids))  summary_tab$Pass_Boost   <- all_ids_norm %in% norm(boost_ids)
if(!is.null(prot_hits))  summary_tab$Is_Protein   <- all_ids_norm %in% norm(prot_hits)
summary_tab$Final_Set <- all_ids_norm %in% final_ids_norm

write_csv(summary_tab, paste0(opt$out, "_intersection_summary.csv"))
cat("Process complete. Results saved with prefix:", opt$out, "\n")