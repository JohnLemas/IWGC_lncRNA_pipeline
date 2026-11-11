#!/usr/bin/bash


# --- public images for the lncRNA annotation pipeline --- #
# These images were pulled from public sources. This script
# serves as a catalogue of tool versions for the pipeline,
# and can be used to pull all necessary public tools.
# Usage: bash README.sh
# ------------------------------------------------------- #

# 1. STAR
singularity pull star_2.7.11b.sif docker://quay.io/biocontainers/star:2.7.11b--h5ca1c30_4

# 2. GFFcompare
singularity pull gffcompare_0.12.6.sif docker://quay.io/biocontainers/gffcompare:0.12.6--h9f5acd7_0

# 3. FastP
singularity pull fastp_0.23.4.sif docker://quay.io/biocontainers/fastp:0.23.4--h5f740d0_0

# 4. FastQC
singularity pull fastqc_0.12.1.sif docker://quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0

# 5. MultiQC
singularity pull docker://quay.io/biocontainers/multiqc:1.14--pyhdfd78af_0

# 6. SAMtools
singularity pull samtools_1.18.sif docker://quay.io/biocontainers/samtools:1.18--h50ea8bc_1

# 7. StringTie
singularity pull stringtie_2.2.1.sif docker://quay.io/biocontainers/stringtie:2.2.1--hecb563c_2

# 8. SeqKit
singularity pull seqkit_2.3.0.sif docker://quay.io/biocontainers/seqkit:2.3.0--h9ee0642_0

# 9. AGAT
singularity pull agat_1.5.0.sif docker://quay.io/biocontainers/agat:1.5.0--pl5321hdfd78af_0
