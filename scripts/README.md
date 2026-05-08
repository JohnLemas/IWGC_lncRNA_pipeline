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
  - random forest Ml tool for lncRNA prediction
- [DIAMOND](#diamond)
  - sequence homology
- [fastp template](#fastp-tempmlate)
  - template `sbatch` script called by `fastp_wrapper.sh`
- [fastp wrapper](#fastp-wrapper)
  - wrapper script to submit all raw `fastq` files for trimming
- [FEELnc classifier](#feelnc-classifier)
  - edits gtf entries for predicted lncRNAs to include classifications
- [FEELnc filter](#feelnc-filter)
  - linear model ML tool for lncRNA candidate identification
- 