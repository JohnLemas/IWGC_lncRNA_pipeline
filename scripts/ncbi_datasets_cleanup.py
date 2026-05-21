#!/usr/bin/env python3
import json
import os
import glob
import sys

# 1. Capture the species name from the command-line argument
if len(sys.argv) < 2:
    print("Error: Missing species name argument!")
    print("Usage: python3 universal_clean.py [species_name]")
    print("Example: python3 universal_clean.py g.max")
    sys.exit(1)

species_prefix = sys.argv[1]

print("Searching for NCBI dataset directory inside /data/GCF*...")

# 2. Automatically find the GCF dataset directory
gcf_dirs = glob.glob("data/GCF*")

if not gcf_dirs:
    print("Error: Could not find any directory matching '/data/GCF*'!")
    sys.exit(1)

target_dir = gcf_dirs[0]
print(f"Found active genome directory: {target_dir}")

# 3. Dynamically resolve file paths inside that directory
fasta_matches = glob.glob(os.path.join(target_dir, "*_genomic.fna"))
input_fasta = fasta_matches[0] if fasta_matches else None
input_gff = os.path.join(target_dir, "genomic.gff")
seq_report_path = os.path.join(target_dir, "sequence_report.jsonl")

# Safety check
missing_files = []
if not input_fasta: missing_files.append("*_genomic.fna")
if not os.path.exists(input_gff): missing_files.append("genomic.gff")
if not os.path.exists(seq_report_path): missing_files.append("sequence_report.jsonl")

if missing_files:
    print(f"Error: Missing required files in {target_dir}: {', '.join(missing_files)}")
    sys.exit(1)

# Define custom output paths using your exact naming specs
output_fasta = os.path.join(target_dir, f"{species_prefix}.chromosomes.fa")
output_gff = os.path.join(target_dir, f"{species_prefix}_mRNA.gff")

print(f"-> Using FASTA:          {os.path.basename(input_fasta)}")
print(f"-> Using GFF:            {os.path.basename(input_gff)}")
print(f"-> Target FASTA Output:  {os.path.basename(output_fasta)}")
print(f"-> Target GFF Output:    {os.path.basename(output_gff)}")

# --- Step 1: Build Chromosome Map ---
print("\nStep 1: Parsing sequence report to map chromosomes...")
chr_mapping = {}

with open(seq_report_path, 'r') as json_file:
    for line in json_file:
        data = json.loads(line)
        if data.get("assignedMoleculeLocationType") == "assigned-molecule" and data.get("role") == "assembled-molecule":
            refseq_acc = data.get("refseqAccession")
            chr_name = data.get("chromosomeName")
            
            if refseq_acc and chr_name:
                chr_mapping[refseq_acc] = f"chr{chr_name}"

if not chr_mapping:
    print("Warning: No primary chromosomes detected in the report data.")
else:
    print(f"Success! Detected {len(chr_mapping)} main chromosomes.")

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

print(f"Saved: {output_fasta}")

# --- Step 3: Clean GFF (Extracting mRNA lines only) ---
print("\nStep 3: Filtering GFF to primary chromosomes and extracting mRNA features...")
with open(input_gff, 'r') as infile, open(output_gff, 'w') as outfile:
    for line in infile:
        if line.startswith('#'):
            outfile.write(line)
            continue
            
        columns = line.split('\t')
        seqid = columns[0]
        feature_type = columns[2]
        
        # Keep the line ONLY if it belongs to a core chromosome AND represents mRNA structures
        # Note: We keep features of type 'mRNA', plus 'exon', 'CDS', 'five_prime_UTR', etc. 
        # to ensure the transcript structures remain fully intact.
        if seqid in chr_mapping:
            if feature_type in ['gene', 'mRNA', 'exon', 'CDS', 'five_prime_UTR', 'three_prime_UTR']:
                columns[0] = chr_mapping[seqid]
                outfile.write('\t'.join(columns))

print(f"Saved: {output_gff}")
print(f"\nProcessing complete! Outputs generated for string ID: '{species_prefix}'")