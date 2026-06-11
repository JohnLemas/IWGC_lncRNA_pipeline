#!/usr/bin/env python3
import argparse
import os
import sys
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description="Map input SRA datasets to high-confidence validated lncRNAs.")
    parser.add_argument("-t", "--tracking", required=True, help="Path to candidate.tracking file")
    parser.add_argument("-l", "--list", required=True, help="Path to sra_file_list.txt")
    parser.add_argument("-f", "--fasta", required=True, help="Path to final high-confidence filtered fasta file")
    parser.add_argument("-o", "--output", required=True, help="Path to output dataset contributions summary CSV")
    
    args = parser.parse_args()

    TRACKING_FILE = args.tracking
    SRA_LIST_FILE = args.list
    FASTA_FILE = args.fasta
    OUTPUT_CSV = args.output

    # 1. Parse your final high-confidence transcript IDs from your filtered FASTA file
    validated_ids = set()
    if not os.path.exists(FASTA_FILE):
        print(f"Error: High-confidence FASTA file not found at {FASTA_FILE}")
        sys.exit(1)

    with open(FASTA_FILE, "r") as f:
        for line in f:
            if line.startswith(">"):
                # Grab the primary sequence identifier token up to the first space
                tx_id = line.strip().split()[0].replace(">", "")
                if "|" in tx_id:
                    tx_id = tx_id.split("|")[-1]
                validated_ids.add(tx_id)

    print(f"Loaded {len(validated_ids)} high-confidence validated transcript IDs to target.")

    # 2. Read SRA list to map tracking file column positions to SRR identifiers
    sra_names = []
    if not os.path.exists(SRA_LIST_FILE):
        print(f"Error: SRA list file not found at {SRA_LIST_FILE}")
        sys.exit(1)

    with open(SRA_LIST_FILE, "r") as f:
        for line in f:
            if line.strip():
                sra_id = os.path.basename(line.strip()).split("_")[0]
                sra_names.append(sra_id)

    print(f"Loaded {len(sra_names)} SRA identifiers from your dataset list.")

    # 3. Initialize contribution counter matrix
    sra_counts = {sra: 0 for sra in sra_names}
    matched_transcripts = 0

    if not os.path.exists(TRACKING_FILE):
        print(f"Error: Tracking file not found at {TRACKING_FILE}")
        sys.exit(1)

    print("Scanning tracking matrix specifically against your high-confidence list...")
    with open(TRACKING_FILE, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 5:
                continue
                
            # Extract names based on GffCompare joint tracking structure
            tcons_raw = parts[0].split("|")[0]
            ref_raw = parts[2]
            sra_columns = parts[4:]  # Sample evaluation strings run from index 4 onward
            
            # Unpack correct candidate transcript identifier matching your FASTA file
            if ref_raw != "-" and ref_raw.strip() != "":
                ref_parts = ref_raw.split("|")
                candidate_id = ref_parts[1] if len(ref_parts) > 1 else ref_parts[0]
            else:
                candidate_id = tcons_raw
                
            # If this transcript row matches one of your 2,570 high-confidence IDs, tally its sources
            if candidate_id in validated_ids:
                matched_transcripts += 1
                for idx, cell in enumerate(sra_columns):
                    if idx >= len(sra_names):
                        break
                    # If the matrix cell contains content other than a dash, this dataset contained it
                    if cell.strip() != "-":
                        sra_counts[sra_names[idx]] += 1

    # 4. Compile, sort, and save into your final summary table
    df_summary = pd.DataFrame([
        {"SRA_Accession": sra, "Validated_lncRNA_Count": count}
        for sra, count in sra_counts.items()
    ])

    df_summary = df_summary.sort_values(by="Validated_lncRNA_Count", ascending=False)
    df_summary.to_csv(OUTPUT_CSV, index=False)

    print(f"\nSuccess! Successfully processed your sequence mapping matrix.")
    print(f"Identified {matched_transcripts} matching structure rows inside the tracking table.")
    print(f"Saved dataset contribution matrix directly to: {OUTPUT_CSV}\n")
    print("Top 10 SRA datasets that contributed the most sequences to your final pool:")
    print(df_summary.head(10).to_string(index=False))

if __name__ == "__main__":
    main()