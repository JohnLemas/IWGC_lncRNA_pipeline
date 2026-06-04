#!/usr/bin/env python3
import argparse
import os
import sys
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description="Parse GffCompare tracking files and map to SRA lists.")
    parser.add_argument("-t", "--tracking", required=True, help="Path to candidate_tracking.tracking file")
    parser.add_argument("-l", "--list", required=True, help="Path to sra_file_list.txt")
    parser.add_argument("-o", "--output", required=True, help="Path to output summary CSV")
    
    args = parser.parse_args()

    TRACKING_FILE = args.tracking
    SRA_LIST_FILE = args.list
    OUTPUT_CSV = args.output

    # 1. Read SRA list to map column indices to SRR identifiers
    sra_names = []
    if not os.path.exists(SRA_LIST_FILE):
        print(f"Error: SRA list file not found at {SRA_LIST_FILE}")
        sys.exit(1)

    with open(SRA_LIST_FILE, "r") as f:
        for line in f:
            if line.strip():
                # Extract the SRR string from the filename path
                sra_id = os.path.basename(line.strip()).split("_")[0]
                sra_names.append(sra_id)

    print(f"Loaded {len(sra_names)} SRA identifiers from your dataset list.")

    # 2. Parse the tracking file line by line
    rows = []
    print("Processing tracking matrix...")

    if not os.path.exists(TRACKING_FILE):
        print(f"Error: Tracking file not found at {TRACKING_FILE}")
        sys.exit(1)

    with open(TRACKING_FILE, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 5:
                continue
                
            # GffCompare tracking columns:
            # Col 2 = Reference Transcript ID (your candidate MSTRG/STRG id)
            # Col 4 = Class code (=, j, etc.)
            candidate_id = parts[2]
            class_code = parts[4]
            
            # Extract sample columns (everything after column 4)
            sra_columns = parts[5:]
            
            # We look for strict matches (=) or compatible splicing structures (j)
            if class_code in ["=", "j"]:
                sras_present = []
                
                # Zip ensures we don't index past the number of loaded SRA names
                for idx, cell in enumerate(sra_columns):
                    if idx >= len(sra_names):
                        break
                    # If the cell contains data (not a dash), this SRA has the sequence
                    if cell.strip() != "-":
                        sras_present.append(sra_names[idx])
                
                if sras_present:
                    rows.append({
                        "Transcript_ID": candidate_id,
                        "Match_Type": class_code,
                        "SRA_Count": len(sras_present),
                        "SRA_IDs": ",".join(sorted(sras_present))
                    })

    # 3. Safe guard empty datasets, then export
    if not rows:
        print("\n[WARNING]: No transcripts matched with class code '=' or 'j'.")
        print("This usually means your candidates file is using slightly different structural IDs than the tracking file's reference anchor column.")
        sys.exit(1)

    df = pd.DataFrame(rows)
    df = df.sort_values(by="SRA_Count", ascending=False)
    df.to_csv(OUTPUT_CSV, index=False)

    print(f"\nSuccess! Saved tracking summary to {OUTPUT_CSV}")
    print(f"Total novel candidates matched: {len(df)}")
    print("\nTop 10 most widely conserved candidates:")
    print(df[["Transcript_ID", "Match_Type", "SRA_Count"]].head(10).to_string(index=False))

if __name__ == "__main__":
    main()