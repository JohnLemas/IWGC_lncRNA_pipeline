#!/bin/bash --login
#SBATCH --job-name=soy_cd-hit
#SBATCH --cpus-per-task=8
#SBATCH --mem=1G
#SBATCH --time=02:00:00
#SBATCH --output=/mnt/scratch/lemasjoh/BASF/logs/%x-%j.out

# Paths
INPUT="/mnt/scratch/lemasjoh/BASF/RNA_central/soy/truth_data/PLncDB/G.max_PLncDB_all.fasta"
OUTPUT="/mnt/scratch/lemasjoh/BASF/RNA_central/soy/truth_data/PLncDB/G.max_PLncDB_all_cd-hit.0.9.fasta"

# module
ml purge
ml load CD-HIT

# command
# -c sequence identity threshold, default 0.9
# -M memory limit (in MB) for the program, default 800; 0 for unlimitted
# -T number of threads, default 1; with 0, all CPUs will be used
# -l length of throw_away_sequences, default 10
cmd="cd-hit -i ${INPUT} -o ${OUTPUT} -c 0.9 -M 800 -T 8 -l 200"

echo -e "run command:\n\t${cmd}\n"

$cmd

echo "done"