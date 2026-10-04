#!/bin/bash
# This is the config file that references the locationsos of all inputs and outputs throughout the pipeline
# Sourcing this file will set all the necessary paths and variables for each script in the pipeline

# Working directory location
BASE="/path/to/project/dir"

# Species prefix for both the working directory and output files
# IWGC standard is first three letters of the family name and the two first letters of the species name
# i.e. Amaranthus tuberculatus: AmaTu or Echinocloa colona: EchCo
SPEC="species"

# All necessary input file paths- they need to live in BASE/SPEC/01_inputs/
GENOME="${BASE}/${SPEC}/01_inputs/${SPEC}_genome.fasta" # Just the chromosome contigs- no unplaced sequences or mito/plastid assemblies
REF_ANN="${BASE}/${SPEC}/01_inputs/${SPEC}_mRNA.gtf" # this has to be a gtf file. convert using gffread -T
TE_GFF3="${BASE}/${SPEC}/01_inputs/${SPEC}_TEanno.gff3" # Direct output from EDTA- rename for convenience
METADATA="${BASE}/${SPEC}/01_inputs/${SPEC}_metadata.txt" # Sequencing ID in column 1. a label designating tissue type and replicate number should be in column 2 (ie sdl_rt_2 or exp_lvs_1)
RAW_FASTQ="${BASE}/${SPEC}/01_inputs/raw_fastq" # Directory containing raw FASTQ files

# First batch of outputs
