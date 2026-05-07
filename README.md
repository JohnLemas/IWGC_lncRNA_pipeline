# IWGC_lncRNA_pipeline
## Expanding lncRNA prediction to nonmodel weedy species
### Pipeline workflow and resources
---
### Table of Contents
- [General overview](#general-overview)
- [Directories](*directories)
- [Citation](#citation)
---
### General overview
1. Obtain 'truth' dataset from PLncDB for a given speices
2. Get a list of pmids associated with the rna-seq projects under the species tab in PLncDB
3. Use EDirect and SRA-Toolkit to get raw data from each of these projects
4. Download a reference genome and annotation file
5. Trim raw fastq files and align trimmed reads to reference
6. Assemble a reference-guided de novo transcriptome using stringtie
7. Obtain candidate lncRNAs (novel transcripts called by stringtie) using grep and gffread
8. Run these candidates through each ML tool
9. Intersect IDs that passed each tool to pbtain high-confidence lncRNAs
10. Use sequence alignment to compare high-confidence candidates to 'truth' dataset
---
### Directories
- def_files
  - `.def` files used to generate required images
- scripts
  - MSU HPCC sbatch scripts used for each step of the pipeline
---
### Citation
This workflow was adapted from Plant-LncPipe for use in *ab initio* lncRNA prediction in weedy species.
- Original pipeline doi: [Plant-LncPipe](https://doi.org/10.1093/hr/uhae041) 
- Pipeline GitHub: [Plant-LncRNA-pipeline-v2](https://github.com/xuechantian/Plant-LncRNA-pipeline-v2)
