#!/bin/bash

TEMP="/mnt/scratch/lemasjoh/BASF/IWGC_lncRNA_pipeline/species/02_scripts/fastp_template.sb"

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
            sbatch "$TEMP" "$accession" "$FOLDER_NAME"
        else
            echo "Skipping non-accession line: $accession"
        fi
    done < "$LIST_FILE"
}

submit_batch "/mnt/scratch/lemasjoh/BASF/IWGC_lncRNA_pipeline/species/01_inputs/species_sra_IDs.txt" "species"