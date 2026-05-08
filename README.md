# IWGC lncRNA Prediction
## Expanding annotations for non-model plants
### John Lemas
### 5-6-2025
---
### Table of contents
- [Working Directory](#working-directory)
- [General overview](#general-overview)
- [Citation](*citation)
---
## Working Directory
- def_files
  - def files for the pipeline. Use these to create necessary images.
- scripts
  - sbatch script templates for hpcc
  - custom scripts for the pipeline
- corn
  - 01_inputs
  - 02_scripts
  - 03_outputs
- A.thaliana
  - 01_inputs
  - 02_scripts
  - 03_outputs
- soy
  - 01_inputs
  - 02_scripts
  - 03_outputs
 ---
## General overview
1. Obtain 'truth' dataset from PLncDB for a given speices
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
 ## Citation
 This pipeline was originally developed by Xue-Chan Tian *et al*:
 Original publication:[DOI](https://doi.org/10.1093/hr/uhae041)
 Author GitHub:[Plant-LncPipe-v2](https://github.com/xuechantian/Plant-LncRNA-pipeline-v2/tree/master)
