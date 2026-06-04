import pandas as pd

TRACKING_FILE = "candidate_tracking.tracking"
SRA_LIST_FILE = "sra_file_list.txt"
OUTPUT_CSV = "z.mays_lncRNA_publication_counts.csv"

# 1. Read your SRA list to know which column belongs to which SRR identifier
sra_names = []
with open(SRA_LIST_FILE, "r") as f:
    for line in f:
        if line.strip():
            # Get just the SRRXXXXXX part out of the full path
            sra_id = line.strip().split("/")[-1].split("_")[0]
            sra_names.append(sra_id)

print(f"Loaded {len(sra_names)} SRA identifiers from your dataset list.")

# 2. Parse the tracking file
rows = []
print("Processing tracking matrix...")

with open(TRACKING_FILE, "r") as f:
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) < 5:
            continue
            
        # GffCompare tracking format:
        # Col 0: Unique tracking ID
        # Col 1: Associated Candidate Gene/Transcript ID
        # Col 2: Reference transcript ID (Your MSTRG ID)
        # Col 4: Class code (=, j, x, u, etc.)
        # Col 5 onwards: Individual SRA columns
        
        candidate_id = parts[2]
        class_code = parts[4]
        sra_columns = parts[5:]
        
        # We only want to count robust structural matches (Exact '=' or Splicing compatible 'j')
        if class_code in ["=", "j"]:
            sras_present = []
            
            for idx, sra_cell in enumerate(sra_columns):
                # If the cell is not "-", it means this SRA contained the transcript
                if sra_cell != "-":
                    sras_present.append(sra_names[idx])
            
            if sras_present:
                rows.append({
                    "Transcript_ID": candidate_id,
                    "Match_Type": class_code,
                    "SRA_Count": len(sras_present),
                    "SRA_IDs": ",".join(sras_present)
                })

# 3. Collapse multiple isoforms if necessary and save
df = pd.DataFrame(rows)
# Group by transcript ID to aggregate counts if a transcript hits multiple tracking categories
df = df.sort_values(by="SRA_Count", ascending=False)
df.to_csv(OUTPUT_CSV, index=False)

print(f"\nSaved successful mapping to {OUTPUT_CSV}!")
print(df[["Transcript_ID", "Match_Type", "SRA_Count"]].head(15).to_string(index=False))
