#!/bin/bash
# =============================================================================
# config.sh - single source of truth for every path in the lncRNA pipeline.
#
# Every step script (00-16) and every submit_*.sh sources this file, so a path is
# defined exactly once and the step that WRITES a file and the step that READS it
# can never drift apart. Change BASE/SPEC (and the input files) here, nowhere else.
#
# Step scripts find this file through $CONFIG (exported by the submit scripts) and
# fall back to ./config.sh in the directory you ran sbatch from.
# =============================================================================

# ---- 1. Project identity ----------------------------------------------------
# Working directory location
BASE="/mnt/scratch/lemasjoh/BASF/lncRNA_met_dev"

# Species prefix for both the working directory and output files
# IWGC standard is first three letters of the family name and the two first letters of the species name
# i.e. Amaranthus tuberculatus: AmaTu or Echinocloa colona: EchCo
SPEC="A.tha"

# Where the pipeline scripts live (auto-detected; the submit scripts set it first)
PIPELINE_DIR="${PIPELINE_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"

# ---- 2. Directory layout ----------------------------------------------------
PROJ="${BASE}/${SPEC}"
INPUT_ROOT="${PROJ}/01_inputs"
OUTPUT_ROOT="${PROJ}/03_outputs"
LOG_DIR="${PROJ}/logs"                  # submit scripts write every SLURM log under here

# Shared resources that live at BASE level (not per species)
IMG_DIR="${BASE}/images"                # singularity images
SCRIPTS_DIR="${BASE}/scripts"           # FEELnc, Plant-LncRNA-pipeline-v2, CPAT logit file
DB_DIR="${BASE}/databases"              # SortMeRNA default db, UniProt/DIAMOND db, barrnap custom db

# Singularity bind mount used by the container steps (10, 12, 13).
# Was hard-coded to /mnt/scratch/lemasjoh; widen this if a tool needs paths outside BASE.
BIND_DIR="${BASE}"

# ---- 3. Input files (live in BASE/SPEC/01_inputs/) -------------------------
GENOME="${INPUT_ROOT}/${SPEC}_genome.fasta"   # Just the chromosome contigs- no unplaced sequences or mito/plastid assemblies
REF_ANN="${INPUT_ROOT}/${SPEC}_mRNA.gtf"      # this has to be a gtf file. convert using gffread -T
TE_GFF3="${INPUT_ROOT}/${SPEC}_TEanno.gff3"   # Direct output from EDTA- rename for convenience
METADATA="${INPUT_ROOT}/${SPEC}_metadata.txt" # TAB-delimited. Sequencing ID in column 1. a label designating tissue type and replicate number should be in column 2 (ie sdl_rt_2 or exp_lvs_1)
RAW_FASTQ="${INPUT_ROOT}/raw_fastq"           # Directory containing raw FASTQ files
# Genome that barrnap scans for rRNA genes (step 02). Point this at the full NCBI datasets assembly
# (nuclear + organelles + unplaced) so mito/plastid rRNA ends up in the custom SortMeRNA db.
# A single fasta; .gz is fine. Defaults to GENOME until you set it, e.g.:
#   BARRNAP_GENOME="${INPUT_ROOT}/ncbi_dataset/data/GCF_XXXXXXXXX.X/GCF_XXXXXXXXX.X_genomic.fna"
BARRNAP_GENOME="${GENOME}"

# ---- 4. Containers / databases ----------------------------------------------
SMR_SIF="${IMG_DIR}/sortmerna.sif"            # also runs barrnap (02) and sortmerna (02, 03)
RBIO_SIF="${IMG_DIR}/rbio_tools.sif"          # R helpers: tpm, track, intersect, intersect_te, finalize, lnc_predict
FEELNC_SIF="${IMG_DIR}/FEELnc2.sif"
CPAT_SIF="${IMG_DIR}/cpat_plant.sif"
LNCBOOST_SIF="${IMG_DIR}/plant_lnc_boost.sif"

# ---- 5. Step outputs, in pipeline order -------------------------------------
# 00 STAR genomeGenerate  -> index (read by 05)
STAR_INDEX="${INPUT_ROOT}/${SPEC}_indexed"
INTRON_LIMITS="${INPUT_ROOT}/${SPEC}_intron_limits.txt"          # 00 also writes this (KEY=value lines, sourceable); 05 reads INTRON_LIMIT from it

# 01 fastp                -> trimmed reads (read by 03, 04)
FASTP_DIR="${OUTPUT_ROOT}/01_fastp"

# 02 barrnap + cd-hit-est -> custom rRNA db + shared SortMeRNA index (read by 03)
SMR_DEFAULT_DB="${DB_DIR}/smr_v4.3_default_db.fasta"   # must be the SAME default db 03 uses
BARRNAP_DIR="${DB_DIR}/${SPEC}_barrnap"
SMR_CUSTOM_DB="${BARRNAP_DIR}/custom_db.fasta"
SMR_IDX_DIR="${BARRNAP_DIR}/idx"

# 03 SortMeRNA            -> rRNA-depleted reads (read by 04, 05)
SMR_DIR="${OUTPUT_ROOT}/02_sortmerna"
SMR_CLEAN_DIR="${SMR_DIR}/clean_reads"
SMR_RRNA_DIR="${SMR_DIR}/rRNA_reads"
SMR_FAILED_DIR="${SMR_DIR}/FAILED_workdirs"

# 04 FastQC + MultiQC     -> one report per stage
QC_RAW_DIR="${RAW_FASTQ}/QC"
QC_TRIM_DIR="${FASTP_DIR}/QC"
QC_CLEAN_DIR="${SMR_DIR}/QC"

# 05 STAR alignment       -> sorted BAMs (read by 06, 08)
STAR_DIR="${OUTPUT_ROOT}/03_STAR"

# 06-09 StringTie
STRINGTIE_DIR="${OUTPUT_ROOT}/04_stringtie"
ASSEMBLY_DIR="${STRINGTIE_DIR}/sample_assemblies"                # 06 writes ${ID}_assembled.gtf; 07 and 15 read
MERGED_GTF="${STRINGTIE_DIR}/${SPEC}_denovo_trans.gtf"           # 07 -> read by 08
CANDIDATES_GTF="${STRINGTIE_DIR}/${SPEC}_candidates.gtf"         # 07 (MSTRG-only) -> read by 09,10,15,16
CANDIDATES_FASTA="${STRINGTIE_DIR}/${SPEC}_candidates.fasta"     # 07 (gffread -w) -> read by 11,12,13,14,15,16
QUANT_DIR="${STRINGTIE_DIR}"                                     # 08 writes ${QUANT_DIR}/${ID}/t_data.ctab; 09 reads
TPM_PREFIX="${STRINGTIE_DIR}/${SPEC}_lncRNA_TPM"                 # 09 -> .matrix.tsv / .summary.tsv
TPM_MATRIX="${TPM_PREFIX}.matrix.tsv"
TPM_SUMMARY="${TPM_PREFIX}.summary.tsv"                          # read by 16 (needs a tissue column in METADATA)

# 10 FEELnc filter
FEELNC_HOME="${SCRIPTS_DIR}/feelnc"                              # FEELnc install (lib/ and scripts/)
FEELNC_OUT_DIR="${OUTPUT_ROOT}/05_FEELnc"
FEELNC_GTF="${FEELNC_OUT_DIR}/${SPEC}_FEELnc.gtf"

# 11 CPAT-plant
CPAT_HEXAMER="${SCRIPTS_DIR}/Plant-LncRNA-pipeline-v2/Model/Plant_Hexamer.tsv"  # hexamer table from the Plant-LncPipe GitHub repo
CPAT_LOGIT="${SCRIPTS_DIR}/Plant.logit.v2.RData"                 # also from Plant-LncPipe, but re-saved to work with cpat.sif- use this one
CPAT_DIR="${OUTPUT_ROOT}/06_cpat"
CPAT_OUT="${CPAT_DIR}/${SPEC}_cpat.output"

# 12 PlantLncBoost
LNCBOOST_SCRIPT_DIR="${SCRIPTS_DIR}/Plant-LncRNA-pipeline-v2/PlantLncBoost/Script"
LNCBOOST_MODEL="${SCRIPTS_DIR}/Plant-LncRNA-pipeline-v2/Model/PlantLncBoost_model.cb"
LNCBOOST_DIR="${OUTPUT_ROOT}/07_lncBoost"
LNCBOOST_FEATURES="${LNCBOOST_DIR}/${SPEC}_features.csv"
LNCBOOST_PREDS="${LNCBOOST_DIR}/${SPEC}_predictions.txt"

# 13 LncFinder
LNCFINDER_MODEL="wheat"                                          # wheat, human, or mouse
LNCFINDER_DIR="${OUTPUT_ROOT}/08_lncFinder"
LNCFINDER_CSV="${LNCFINDER_DIR}/${SPEC}_lncFinder_results.csv"

# 14 DIAMOND blastx vs UniProt
DIAMOND_DB="${DB_DIR}/uniprot_db/uniprot_out.dmnd"               # diamond makedb --in uniprot_sprot.fasta -d uniprot_out.dmnd
DIAMOND_DIR="${OUTPUT_ROOT}/09_diamond"
DIAMOND_OUT="${DIAMOND_DIR}/${SPEC}_diamond_out.txt"

# 15 gffcompare tracking (tissue support)
TRACE_DIR="${OUTPUT_ROOT}/10_trace"
TRACE_LIST="${TRACE_DIR}/file_list.txt"
TRACKING_SUMMARY="${TRACE_DIR}/${SPEC}_tracking_summary.tsv"     # read by 16

# 16 intersect + classify + finalize
INT_CLASS_DIR="${OUTPUT_ROOT}/11_int_class"
INT_PFX="${INT_CLASS_DIR}/${SPEC}"                               # file prefix used by every 16 sub-step
FINAL_ANNOTATION="${INT_PFX}_final_annotation.tsv"               # the end product

# ---- 6. Side utility: truth set clustering (cd-hit.sb, not part of the DAG) --
TRUTH_DIR="${BASE}/truth_data/${SPEC}"
TRUTH_LNCRNA_FA="${TRUTH_DIR}/${SPEC}_PLncDB_lncRNA.fa"          # all PLncDB lncRNA sequences for this species
TRUTH_CDHIT_FA="${TRUTH_DIR}/${SPEC}_PLncDB_cd-hit.0.9.fa"       # representatives after clustering at 90% identity
