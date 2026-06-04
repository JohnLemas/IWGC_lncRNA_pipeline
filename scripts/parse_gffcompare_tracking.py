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
                sra_id = os.path.basename(line.strip()).split("_")[0]
                sra_names.append(sra_id)

    print(f"Loaded {len(sra_names)} SRA identifiers from your dataset list.")

    # 2. Parse the tracking file line by line
    rows = []
    print("Processing tracking matrix...")

    if not os.path.exists(TRACKING_FILE):
        print(f"Error: Tracking file not found at {TRACKING_FILE}")
        sys.exit(1)

    # Biological codes to capture structural overlap/containment
    TARGET_CODES = {"=", "j", "c", "m"}

    with open(TRACKING_FILE, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 4:
                continue
                
            # Adjusted Indices based on your file layout:
            # Col 0: TCONS ID
            # Col 1: XLOC ID
            # Col 2: MSTRG ID (Reference column)
            # Col 3: Class Code
            # Col 4 onwards: SRA sample columns
            tcons_raw = parts[0].split("|")[0]
            ref_raw = parts[2]
            class_code = parts[3].strip()  # Fix: shifted to index 3
            sra_columns = parts[4:]        # Fix: shifted to index 4
            
            # Extract clean candidate identifier (e.g., from MSTRG.64|MSTRG.64.4 -> MSTRG.64.4)
            if ref_raw != "-" and ref_raw.strip() != "":
                ref_parts = ref_raw.split("|")
                candidate_id = ref_parts[1] if len(ref_parts) > 1 else ref_parts[0]
            else:
                candidate_id = tcons_raw
            
            if class_code in TARGET_CODES:
                sras_present = []
                
                for idx, cell in enumerate(sra_columns):
                    if idx >= len(sra_names):
                        break
                    if cell.strip() != "-":
                        sras_present.append(sra_names[idx])
                
                if sras_present:
                    rows.append({
                        "Transcript_ID": candidate_id,
                        "Match_Type": class_code,
                        "SRA_Count": len(sras_present),
                        "SRA_IDs": ",".join(sorted(sras_present))
                    })

    if not rows:
        print("\n[ERROR]: No rows extracted. Please double check column positions or codes.")
        sys.exit(1)

    df = pd.DataFrame(rows)
    df = df.sort_values(by="SRA_Count", ascending=False)
    df.to_csv(OUTPUT_CSV, index=False)

    print(f"\nSuccess! Saved tracking summary to {OUTPUT_CSV}")
    print(f"Total novel candidates matched: {len(df)}")
    print("\nTop 10 most widely conserved candidates across accessions:")
    print(df[["Transcript_ID", "Match_Type", "SRA_Count"]].head(10).to_string(index=False))

if __name__ == "__main__":
    main()