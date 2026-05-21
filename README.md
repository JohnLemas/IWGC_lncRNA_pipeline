# IWGC lncRNA Prediction
## Expanding annotations for non-model plants
### John Lemas
### 5-6-2025
---
## Citation
 This pipeline was originally developed by Xue-Chan Tian *et al*:
 Original publication:[DOI](https://doi.org/10.1093/hr/uhae041)
 Author GitHub:[Plant-LncPipe-v2](https://github.com/xuechantian/Plant-LncRNA-pipeline-v2/tree/master)
---
### Table of contents
- [General overview](#general-overview)
- [Working Directory](#working-directory)
- [Gathering inputs](#gathering-inputs)
- [Set up](#set-up)
---
## General overview
1. Obtain 'truth' dataset from [PLncDB](https://www.tobaccodb.org/plncdb/) for a given speices
2. Generate a `PMID` list associated with RNA-seq projects under the 'species' tab in PLncDB
3. Use `EDirect` and `SRA-Toolkit` to get raw data from each of these projects
4. Download a reference genome and annotation file
5. Trim raw `fastq` files and align trimmed reads to the reference genome
6. Assemble a reference-guided *de novo* transcriptome using `StringTie`
7. Obtain candidate lncRNAs using `grep` and `gffread`
8. Run these candidates through `FEELnc`, `lncBoost`, `lncFinder`, `cpat`, and `DIAMOND`
9. Intersect IDs that passed every tool to obtain high-confidence lncRNAs
10. Use sequence alignment to compare high-confidence *ab initio* candidates to 'truth' dataset
 ---
## Working Directory
- def_files
  - definition files to create necessary imagesfor each tool
- scripts
  - sbatch script templates for hpcc
  - custom `R` and `python` scripts
- species
  - 01_inputs
    - PMID metadata
    - list of sra IDs saved as `species_sra_IDs.txt`
    - reference genome (`species.chromosomes.fa`) and annotation (`species_mRNA.gtf`)
  - 02_scripts
    - scripts copied from parent script directory (see `scripts/README.md`)
  - 03_outputs
    - pipeline output directories
- truth_data
  - lncRNA datasets gathered from [PLncDB](https://www.tobaccodb.org/plncdb/)
  - species directories
    - individual database lncRNA `fasta` and `gff3` files
    - full concatonated lncRNA files
    - condensed lncRNA `fasta` files (see `scripts/README.md` `cd-hit` step)
 ---
### Gathering inputs
- Navigate to [PLncDB](https://www.tobaccodb.org/plncdb/)
- Gather PMID information from the RNA-Seq table for a given species
  - Use EDirect to gather SRA accession numbers from PMIDs; save as list `species_sra_IDs.txt`
- Gather truth sequences from the downloads>species page on PLncDB
  - Download the lncRNA fasta files and gff3 files
- Create working directory
  - Clone this repository. Clone [Plant-lncPipe](https://github.com/xuechantian/Plant-LncRNA-pipeline-v2) and [FEELnc](https://github.com/tderrien/feelnc) into the `scripts` directory.
  - use `singularity build` to build images from `def_files/*`
  - Create the working directory for the species
  - Copy scripts into `species/02_scripts` and run the `sed`  command (see `scripts/README.md`)
---
## Set up
### Prepare the reference genome assembly
The NCBI tool `datasets-cli` is an easy way to download a reference genome via command line; however the assembly usually contains messy fasta headers and gff files. The following steps can be used to obtain the latest reference assembly for a species and clean it for use in this pipeline.

```bash
# Create an environment with datasets-cli
ml purge
ml load Miniforge3
conda create -n datasets -c conda-forge ncbi-datasets-cli
conda activate datasets

# General datasets command
datasets download genome taxon "[Scientific Name]" --reference --include genome,gff3
unzip ncbi_dataset.zip # This will result in a directory, an md5sum file, and a README.md file
```

The script `ncbi_datasets_cleanup.py` pulls full chromosomes from the genome, filters the gff file for mRNA entries, renames chromosome entries for both files so that they match, then outputs the two required downstream files: `species.chromosomes.fa` & `species.gene.gff`.
```bash
# Clean up the raw download
cd ncbi_dataset/ # move into the unzipped directory
python3 /path/to/ncbi_datasets_cleanup.py [species] # run the cleanup script with the chosen species designation as the first argument

# Move into the GCF folder inside ncbi_dataset/data/ and use gffread to convert the gff to gtf
ml purge
ml load gffread
gffread -T -o species_mRNA.gtf species.gene.gff
```
Move the new `.fa` and `.gtf` files into the `01_inputs` directory for indexing, and remove the ncbi_dataset directory, README file, and md5sum file from the datasets download. Optional: keep the compressed datasets download for later reference.

### Create the UniProt database for DIAMOND:
```bash
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz

ml purge
ml load DIAMOND
diamond makedb --in uniprot_sprot.fasta -d uniprot_out
```
---