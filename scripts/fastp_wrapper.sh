#!/bin/bash

# Function to submit jobs
submit_batch() {
    LIST_FILE=$1
    FOLDER_NAME=$2
    
    echo "Submitting $FOLDER_NAME batch..."
    while read -r accession; do
        # 1. Skip empty lines
        [ -z "$accession" ] && continue
        
        # 2. Only submit if it looks like a real SRA ID (SRR123, ERR123, etc.)
        if [[ $accession =~ ^[SED]RR ]]; then
            sbatch fastp_template.sb "$accession" "$FOLDER_NAME"
        else
            echo "Skipping non-accession line: $accession"
        fi
    done < "$LIST_FILE"
}

# Run for Soy
submit_batch "/mnt/scratch/lemasjoh/BASF/RNA_central/soy/sra_accessions/g.max_all_accessions.txt" "soy"

# Run for Corn
submit_batch "/mnt/scratch/lemasjoh/BASF/RNA_central/corn/sra_accessions/z.mays_all_accessions.txt" "corn"