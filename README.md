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
- [Working Directory](#working-directory)
- [General overview](#general-overview)
- [Gathering inputs](#gathering-inputs)
- [Set-up commands](#set-up-commands)
---
## Working Directory
- def_files
  - def files for the pipeline. Use these to create necessary images.
- scripts
  - sbatch script templates for hpcc
  - custom scripts for the pipeline
- species
  - 01_inputs
    - pmid metadata
    - list of sra IDs saved as `species_sra_IDs.txt`
    - reference genome (`species.fa`) and annotation (`species_mRNA.gtf)
  - 02_scripts
    - all necessary scripts copied from parent script directory (edited using `sed`; see scripts dir)
  - 03_outputs
    - outputs from all 10 steps
- truth_data
  - lncRNA datasets gathered from [PLncDB](https://www.tobaccodb.org/plncdb/)
 ---
## General overview
1. Obtain 'truth' dataset from [PLncDB](https://www.tobaccodb.org/plncdb/) for a given speices
2. get a list of pmids associated with the rna-seq projects under the species tab in PLncDB
3. Use EDirect and SRA-Toolkit to get raw data from each of these projects
4. Download a reference genome and annotation file
5. trim raw fastq files and align trimmed reads to reference
6. assemble a reference-guided de novo transcriptome using stringtie
7. obtain candidate lncRNAs (novel transcripts called by stringtie) using grep and gffread
8. run these candidates through each ML tool
9. intersect IDs that passed each tool to pbtain high-confidence lncRNAs
10. use sequence alignment to compare high-confidence candidates to 'truth' dataset
 ---
### Gathering inputs
- Navigate to [PLncDB](https://www.tobaccodb.org/plncdb/)
- Gather PMID information from the RNA-Seq table for a given species
  - Use EDirect to gather SRA accession numbers from PMIDs; save as list `species_sra_IDs.txt`
- Gather truth and reference genome data from the downloads tab on PLncDB
  - Download the lncRNA fasta files, the genome, and all of the gff3 files
  - rename the genome `species.fa` and the annotation `species_mRNA.gff3`
- Create working directory
  - Clone this repository. Clone [Plant-lncPipe](https://github.com/xuechantian/Plant-LncRNA-pipeline-v2) and [FEELnc](https://github.com/tderrien/feelnc) into the scripts directory.
  - use `singularity build` to build images from `def_files`
  - Create the working directory for the species
  - Copy the scripts into `species/02_scripts` and run the `sed`  command
  - make sure `species_sra_IDs.txt`, `species.fa`, and `species_mRNA.gff3` are in the `species/01_inputs` directory
---
### Set-up
1. Prep the reference genome
```bash
# First, double check chromosome designations in the reference genome
grep ">" species.fa
cut -f 1 species_mRNA.gffs | sort | uniq
# If needed, use sed to ensure consistent naming
```

2. Annotation must be in gtf format for stringtie
```bash
ml purge
ml load gffread
gffread -T -o species_mRNA.gtf species_mRNA.gff3
```

3. Create the UniProt database for DIAMOND
```bash
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz
diamond makedb --in uniprot_sprot.fasta -d uniprot_out
```

4. Double check paths for `star_genomegenerate.sb` then index the genome into `species/01_inputs/species_indexed`

