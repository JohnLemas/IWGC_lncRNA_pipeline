#!/bin/bash

# Function to submit jobs
submit_batch() {
    LIST_FILE=$1
    FOLDER_NAME=$2
    
    echo "Submitting $FOLDER_NAME batch..."
    while read -r accession; do
        [ -z "$accession" ] && continue
        # Pass both the ID and the folder name to the slurm script
        sbatch fq_dump_template.sb "$accession" "$FOLDER_NAME"
    done < "$LIST_FILE"
}

# Run for Soy
submit_batch "/mnt/scratch/lemasjoh/BASF/RNA_central/soy/sra_accessions/g.max_all_accessions.txt" "soy"

# Run for Corn
submit_batch "/mnt/scratch/lemasjoh/BASF/RNA_central/corn/sra_accessions/z.mays_all_accessions.txt" "corn"