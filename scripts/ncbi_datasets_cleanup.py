#!/usr/bin/env python3
import os
import glob
import sys
import re

# 1. Capture the species name from the command-line argument
if len(sys.argv) < 2:
    print("Error: Missing species name argument!")
    print("Usage: python3 ncbi_datasets_cleanup.py [species_name]")
    print("Example: python3 ncbi_datasets_cleanup.py g.max")
    sys.exit(1)

species_prefix = sys.argv[1]

print("Searching for NCBI dataset directory structure...")

# 2. Find the local GCF dataset directory
gcf_dirs = glob.glob("data/GCF*")

if not gcf_dirs:
    print("Error: Could not find any directory matching 'data/GCF*'!")
    print(f"Current working directory is: {os.getcwd()}")
    print("Are you standing inside the 'ncbi_dataset' folder?")
    sys.exit(1)

target_dir = gcf_dirs[0]
print(f"Found active genome directory: {target_dir}")

# 3. Dynamically resolve file paths inside that specific GCF folder
fasta_matches = glob.glob(os.path.join(target_dir, "*_genomic.fna"))
input_fasta = fasta_matches[0] if fasta_matches else None
input_gff = os.path.join(target_dir, "genomic.gff")

# Safety check
missing_files = []
if not input_fasta: missing_files.append("*_genomic.fna")
if not os.path.exists(input_gff): missing_files.append("genomic.gff")

if missing_files:
    print(f"Error: Missing required files in {target_dir}: {', '.join(missing_files)}")
    sys.exit(1)

# Define custom output paths using your naming specifications
output_fasta = os.path.join(target_dir, f"{species_prefix}.chromosomes.fa")
output_gff = os.path.join(target_dir, f"{species_prefix}.gene.gff")

print(f"-> Found FASTA:          {os.path.basename(input_fasta)}")
print(f"-> Found GFF:            {os.path.basename(input_gff)}")
print(f"-> Target FASTA Output:  {output_fasta}")
print(f"-> Target GFF Output:    {output_gff}")

# --- Step 1: Extract Chromosome Map directly from FASTA headers ---
print("\nStep 1: Parsing FASTA headers to map primary chromosomes...")
chr_mapping = {}

# Regex to find things like "chromosome 1" or "chromosome X" but ignore "unlocalized genomic scaffold"
chrom_regex = re.compile(r"chromosome\s+([A-Za-z0-9]+)")

with open(input_fasta, 'r') as infile:
    for line in infile:
        if line.startswith('>'):
            header = line.strip()
            # Skip unlocalized, unplaced, or organelle scaffolds explicitly
            if "unlocalized" in header.lower() or "unplaced" in header.lower() or "mitochondrion" in header.lower() or "chloroplast" in header.lower():
                continue
                
            # Search for the chromosome keyword and number
            match = chrom_regex.search(header)
            if match:
                accession = header.split()[0].replace('>', '')
                chrom_num = match.group(1)
                chr_mapping[accession] = f"chr{chrom_num}"

if not chr_mapping:
    print("Error: No primary chromosomes could be parsed from the FASTA headers!")
    sys.exit(1)
else:
    print(f"Success! Detected {len(chr_mapping)} main chromosomes.")
    # Show a brief preview of what was found
    preview = list(chr_mapping.items())[:3]
    for acc, chrom in preview:
        print(f"  Mapping: {acc} -> {chrom}")

# --- Step 2: Clean FASTA ---
print("\nStep 2: Filtering and renaming FASTA file...")
with open(input_fasta, 'r') as infile, open(output_fasta, 'w') as outfile:
    write_sequence = False
    for line in infile:
        if line.startswith('>'):
            accession = line.split()[0].replace('>', '')
            if accession in chr_mapping:
                outfile.write(f">{chr_mapping[accession]}\n")
                write_sequence = True
            else:
                write_sequence = False
        else:
            if write_sequence:
                outfile.write(line)

print(f"Saved: {os.path.basename(output_fasta)}")

# --- Step 3: Clean GFF (Extracting transcript structures only) ---
print("\nStep 3: Filtering GFF to primary chromosomes and extracting transcript structures...")
with open(input_gff, 'r') as infile, open(output_gff, 'w') as outfile:
    for line in infile:
        if line.startswith('#'):
            outfile.write(line)
            continue
            
        columns = line.split('\t')
        seqid = columns[0]
        feature_type = columns[2]
        
        if seqid in chr_mapping:
            # Captures standard structural elements of coding and non-coding transcripts for alignment accuracy
            if feature_type in ['gene', 'mRNA', 'exon', 'CDS', 'five_prime_UTR', 'three_prime_UTR', 'lnc_RNA']:
                columns[0] = chr_mapping[seqid]
                outfile.write('\t'.join(columns))

print(f"Saved: {os.path.basename(output_gff)}")
print(f"\nProcessing complete! Outputs generated successfully for string ID: '{species_prefix}'")