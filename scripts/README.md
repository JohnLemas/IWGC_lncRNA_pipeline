# Template script directory
### John Lemas
### 5-8-2026
---
This direcotry contains all scripts used to predict and verify lncRNAs in a given species. Scripts can be copied into a species' sub-directory for editing.
```bash
# use a wild card to correct file paths and file names quickly for each sbatch script. Example:
cp /your/file/path/scripts/* /your/file/path/A.thaliana/02_scripts
sed -i 's/[species]/A.thaliana/g' /your/file/path/A.thaliana/02_scripts*sb
```
---
### Contents
- [Alignment array](#alignment-array)
  - Array `SLURM` job that uses an sra list to align trimmed `fastq` files
- [Bam indexing](#bam-indexing)
  - Array `SLURM` job to index `bam` files
- [CD-HIT](#cd-hit)
  - removes duplicate sequences from truth dataset
- [CPAT plant](#cpat-plant)
  - Logistic regression Ml tool for lncRNA prediction
- [DIAMOND](#diamond)
  - protein database homology search to eliminate candidates
- [fastp](#fastp)
  - template `sbatch` script called by `fastp_wrapper.sh`
  - wrapper script to submit all raw `fastq` files for trimming
- [FEELnc](#feelnc)
  - Filter: decision tree ensemble ML tool for lncRNA candidate identification
  - Classifier: edits gtf entries for predicted lncRNAs to include classifications
- [Intersect IDs](#intersect-ids)
  - Custom R script and sbatch script to combine outputs from up to four Ml tools and from DIAMOND
- [lncBoost](#lncboost)
  -  XGBoost model - iterative ML 'gradient boosting' model
- [lncFinder](#lncfinder)
  - support vector machine (SVM) geometric Ml tool - good for finding patterns in nonlinear data
  - image, r script, and sbatch script for this step
- [Minimap2](#minimap2)
  - mapping predictions against truth sequences for pipeline analysis
---
### Alignment array

---
### Bam indexing

---
### CD-HIT

---
### CPAT plant

---
### DIAMOND

---
### fastp

---
### FEELnc

---
### Intersect IDs

---
### lncBoost

---
### lncFinder

---
### Minimap2