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

# --- Step 3: Clean and Normalize to Valid GFF3 ---
print("\nStep 3: Filtering GFF to primary chromosomes and formatting for GFF3 specification...")

# Helper function to parse GFF attribute string into a python dictionary
def parse_attributes(attr_string):
    attrs = {}
    for item in attr_string.strip().split(';'):
        if '=' in item:
            key, val = item.split('=', 1)
            attrs[key.strip()] = val.strip()
    return attrs

# Track mappings so child features (exons) inherit their parent gene's true name
tx_to_gene_name = {}
gene_id_to_name = {}

# Pass 1: Scan file to pre-populate parent gene names
with open(input_gff, 'r') as infile:
    for line in infile:
        if line.startswith('#') or not line.strip():
            continue
        cols = line.split('\t')
        if len(cols) < 9 or cols[0] not in chr_mapping:
            continue
            
        feature_type = cols[2]
        attrs = parse_attributes(cols[8])
        
        if feature_type == 'gene':
            g_id = attrs.get('ID')
            g_name = attrs.get('gene', attrs.get('Name', g_id))
            if g_id:
                gene_id_to_name[g_id] = g_name
        elif feature_type in ['mRNA', 'lnc_RNA', 'transcript']:
            tx_id = attrs.get('ID')
            parent_gene = attrs.get('Parent')
            if tx_id and parent_gene:
                # Find the gene name associated with this parent ID
                g_name = gene_id_to_name.get(parent_gene, parent_gene)
                tx_to_gene_name[tx_id] = g_name

# Pass 2: Filter lines and build clean GFF3 records
with open(input_gff, 'r') as infile, open(output_gff, 'w') as outfile:
    # Standard GFF3 Header
    outfile.write("##gff-version 3\n")
    
    infile.seek(0) # Reset file read pointer to the beginning
    for line in infile:
        if line.startswith('#'):
            continue
            
        columns = line.strip().split('\t')
        if len(columns) < 9:
            continue
            
        seqid = columns[0]
        feature_type = columns[2]
        
        if seqid in chr_mapping:
            # Captures standard structural elements of coding and non-coding transcripts
            if feature_type in ['gene', 'mRNA', 'exon', 'CDS', 'five_prime_UTR', 'three_prime_UTR', 'lnc_RNA']:
                attrs = parse_attributes(columns[8])
                new_attrs = []
                
                # Re-inject/Keep crucial structural tags
                if 'ID' in attrs:
                    new_attrs.append(f"ID={attrs['ID']}")
                if 'Parent' in attrs:
                    new_attrs.append(f"Parent={attrs['Parent']}")
                
                # Explicitly resolve parent gene mapping for gffread's tracking mechanism
                resolved_gene = "unknown_gene"
                if feature_type == 'gene':
                    resolved_gene = attrs.get('gene', attrs.get('Name', attrs.get('ID')))
                elif feature_type in ['mRNA', 'lnc_RNA']:
                    resolved_gene = tx_to_gene_name.get(attrs.get('ID'), attrs.get('Parent', 'unknown_gene'))
                else: # Child components (exon, CDS, UTR)
                    resolved_gene = tx_to_gene_name.get(attrs.get('Parent'), 'unknown_gene')
                
                new_attrs.append(f"gene={resolved_gene}")
                new_attrs.append(f"Name={resolved_gene}")
                
                # Reassemble line
                columns[0] = chr_mapping[seqid]
                columns[8] = ";".join(new_attrs)
                outfile.write('\t'.join(columns) + '\n')

print(f"Saved clean GFF3: {os.path.basename(output_gff)}")
print(f"\nProcessing complete! Outputs generated successfully for string ID: '{species_prefix}'")