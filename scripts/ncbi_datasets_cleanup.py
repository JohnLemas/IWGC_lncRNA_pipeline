#!/usr/bin/env python3
import os
import glob
import sys
import json

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
# NOTE: NCBI datasets also produces a "cds_from_genomic.fna" file alongside
# the real genome FASTA, and a naive "*_genomic.fna" glob matches BOTH --
# glob.glob() does not guarantee ordering, so taking [0] could silently grab
# the wrong file. Filter explicitly to the genome assembly file, which is
# always named "<GCF accession>_<assembly name>_genomic.fna".
all_fna_matches = glob.glob(os.path.join(target_dir, "*_genomic.fna"))
fasta_matches = [f for f in all_fna_matches if not os.path.basename(f).startswith("cds_from_")]
if len(fasta_matches) > 1:
    print(f"Warning: multiple candidate genome FASTA files found in {target_dir}:")
    for f in fasta_matches:
        print(f"  {os.path.basename(f)}")
    print(f"Using the first one: {os.path.basename(fasta_matches[0])}")
input_fasta = fasta_matches[0] if fasta_matches else None
input_gff = os.path.join(target_dir, "genomic.gff")

report_matches = glob.glob(os.path.join(target_dir, "sequence_report.jsonl"))
input_report = report_matches[0] if report_matches else None

# Safety check
missing_files = []
if not input_fasta: missing_files.append("*_genomic.fna")
if not os.path.exists(input_gff): missing_files.append("genomic.gff")
if not input_report: missing_files.append("sequence_report.jsonl")

if missing_files:
    print(f"Error: Missing required files in {target_dir}: {', '.join(missing_files)}")
    sys.exit(1)

# Define custom output paths using your naming specifications
output_fasta = os.path.join(target_dir, f"{species_prefix}.chromosomes.fa")
# Output direct to GTF format now
output_gtf = os.path.join(target_dir, f"{species_prefix}_mRNA.gtf")

print(f"-> Found FASTA:          {os.path.basename(input_fasta)}")
print(f"-> Found GFF:            {os.path.basename(input_gff)}")
print(f"-> Found Sequence Report:{os.path.basename(input_report)}")
print(f"-> Target FASTA Output:  {output_fasta}")
print(f"-> Target GTF Output:    {output_gtf}")

# --- Step 1: Build Chromosome Map from the Sequence Report (authoritative) ---
print("\nStep 1: Parsing sequence report to map primary chromosomes...")
chr_mapping = {}

# sequence_report.jsonl has one JSON object per line, one line per sequence
# (chromosome, scaffold, organelle, etc.) in the assembly. We only want
# records where role == "assembled-molecule" AND
# assignedMoleculeLocationType == "Chromosome" -- this naturally excludes
# unplaced/unlocalized scaffolds, mitochondria, chloroplasts, plasmids, etc.
# without any fragile keyword matching on free-text headers.
with open(input_report, 'r') as infile:
    for line in infile:
        line = line.strip()
        if not line:
            continue
        try:
            row = json.loads(line)
        except json.JSONDecodeError:
            continue

        role = row.get('role', '')
        mol_type = row.get('assignedMoleculeLocationType', '')
        accession = row.get('refseqAccession', '')
        chrom_name = row.get('chrName', '')

        if role == 'assembled-molecule' and mol_type == 'Chromosome':
            if accession:
                chr_mapping[accession] = f"chr{chrom_name}"

if not chr_mapping:
    print("Error: No primary chromosomes could be parsed from the sequence report!")
    print("Check that records have 'role' and 'assignedMoleculeLocationType' fields.")
    sys.exit(1)
else:
    print(f"Success! Detected {len(chr_mapping)} main chromosomes.")
    preview = list(chr_mapping.items())[:3]
    for acc, chrom in preview:
        print(f"  Mapping: {acc} -> {chrom}")

# Defensive check: make sure we don't have two accessions mapping to the same
# output chromosome name -- this is exactly the bug class that caused the
# original "identical sequence IDs" failure. Fail loudly instead of silently
# collapsing two sequences into one.
seen_targets = {}
collisions = []
for acc, chrom in chr_mapping.items():
    if chrom in seen_targets:
        collisions.append((seen_targets[chrom], acc, chrom))
    else:
        seen_targets[chrom] = acc

if collisions:
    print("\nError: Duplicate target chromosome names detected! Refusing to proceed.")
    for acc1, acc2, chrom in collisions:
        print(f"  Both {acc1} and {acc2} would map to {chrom}")
    sys.exit(1)

# Also confirm every accession we plan to pull actually exists in the FASTA,
# so a mismatched assembly report vs FASTA pairing fails loudly rather than
# silently producing an incomplete/empty output.
fasta_accessions = set()
with open(input_fasta, 'r') as infile:
    for line in infile:
        if line.startswith('>'):
            fasta_accessions.add(line[1:].split()[0])

missing_in_fasta = [acc for acc in chr_mapping if acc not in fasta_accessions]
if missing_in_fasta:
    print("\nError: The following accessions from the sequence report were not found in the FASTA:")
    for acc in missing_in_fasta:
        print(f"  {acc} (expected -> {chr_mapping[acc]})")
    print("This usually means the sequence report and FASTA are from different assembly versions.")
    sys.exit(1)

# --- Step 2: Clean FASTA ---
print("\nStep 2: Filtering and renaming FASTA file...")
with open(input_fasta, 'r') as infile, open(output_fasta, 'w') as outfile:
    write_sequence = False
    for line in infile:
        if line.startswith('>'):
            accession = line[1:].split()[0]
            if accession in chr_mapping:
                outfile.write(f">{chr_mapping[accession]}\n")
                write_sequence = True
            else:
                write_sequence = False
        else:
            if write_sequence:
                outfile.write(line)

print(f"Saved: {os.path.basename(output_fasta)}")

# Final verification: confirm the written FASTA has exactly as many sequences
# as chromosomes we mapped, and that every header is unique. This catches any
# remaining edge case before EDTA (or anything else) ever sees the file.
written_headers = []
with open(output_fasta, 'r') as outfile:
    for line in outfile:
        if line.startswith('>'):
            written_headers.append(line.strip()[1:])

if len(written_headers) != len(chr_mapping):
    print(f"Error: Expected {len(chr_mapping)} sequences in output FASTA, found {len(written_headers)}.")
    sys.exit(1)
if len(written_headers) != len(set(written_headers)):
    print("Error: Output FASTA contains duplicate sequence IDs after writing. Aborting.")
    sys.exit(1)
print(f"Verified: {len(written_headers)} unique chromosome sequences written.")

# --- Step 3: Direct Clean-to-GTF Generation ---
print("\nStep 3: Generating FEELnc-compliant GTF directly from database...")

def parse_attributes(attr_string):
    attrs = {}
    for item in attr_string.strip().split(';'):
        if '=' in item:
            key, val = item.split('=', 1)
            attrs[key.strip()] = val.strip()
    return attrs

tx_to_gene_name = {}
gene_id_to_name = {}

# Pass 1: Build the global relationship dictionary
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
        elif feature_type in ['mRNA', 'lnc_RNA', 'transcript', 'tRNA', 'rRNA']:
            tx_id = attrs.get('ID')
            parent_gene = attrs.get('Parent')
            if tx_id and parent_gene:
                g_name = gene_id_to_name.get(parent_gene, parent_gene)
                tx_to_gene_name[tx_id] = g_name

# Pass 2: Print explicitly structured GTF features directly
with open(input_gff, 'r') as infile, open(output_gtf, 'w') as outfile:
    outfile.write("#gtf-version 2.2\n")
    infile.seek(0)

    for line in infile:
        if line.startswith('#') or not line.strip():
            continue

        columns = line.strip().split('\t')
        if len(columns) < 9:
            continue

        seqid = columns[0]
        feature_type = columns[2]

        if seqid in chr_mapping:
            # Standard structural features required by FEELnc
            if feature_type in ['mRNA', 'exon', 'CDS', 'five_prime_UTR', 'three_prime_UTR', 'lnc_RNA']:
                attrs = parse_attributes(columns[8])

                # 1. Resolve the clean gene_id (The exact locus name like LOC103641878)
                if feature_type in ['mRNA', 'lnc_RNA']:
                    tx_id = attrs.get('ID', 'unknown_tx')
                    resolved_gene = tx_to_gene_name.get(tx_id, attrs.get('Parent', 'unknown_gene'))
                else:
                    tx_id = attrs.get('Parent', 'unknown_tx')
                    resolved_gene = tx_to_gene_name.get(tx_id, 'unknown_gene')

                # 2. Re-label features to comply with standard GTF nomenclature
                if feature_type == 'lnc_RNA':
                    feature_type = 'mRNA'

                # 3. Formulate the explicit GTF block layout for column 9
                # Every line gets both gene_id and transcript_id, without exception.
                gtf_attributes = f'gene_id "{resolved_gene}"; transcript_id "{tx_id}"; gene_name "{resolved_gene}";'

                # 4. Modify columns and save
                columns[0] = chr_mapping[seqid]
                columns[2] = feature_type
                columns[8] = gtf_attributes

                outfile.write('\t'.join(columns) + '\n')

print(f"Saved compliant GTF: {os.path.basename(output_gtf)}")
print(f"\nProcessing complete! Direct output generated successfully for: '{species_prefix}'")