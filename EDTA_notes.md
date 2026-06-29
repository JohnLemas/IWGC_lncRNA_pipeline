# IWGC lncRNA Prediction
## Expanding annotations for non-model plants
### John Lemas
### 6-29-2026
---
## TE content for *Z. mays* and *G. max*
lncRNAs are further classified via their context against annotated transposable elements in the genome. This is partly due to the regulatory role that lncRNAs present, and thier role in epigenetic DNA modifications. TEs are typically heavily silenced in the genome, but locating lncRNAs in proximinity to known TEs suggests potiential *cis* regulation between the element and the TE. I had a difficult time locating trustworthy TE annotaitons from these two reference assemblies, so I decided to annotate them myself. The perk of this process is that I know that the same reference assemblies for both genomes will bo used for both lncRNA annotation and TE annotation. 

---
## Running EDTA on the watership server
I wanted to vet the EDTA singularity image that would eventually be put into the IWGC annotation pipeline, and after completing the analysis I'm glad that I did. Versions of EDTA 2+ gave me lots of trouble, specifically at the TIR candidate search phase. After Claude and I wrestled with resource allocation (and learning about swap space) I discovered that the newest versions of EDTA may be fundamentally broken. I pulled EDTA version 1.9.6 and completed the *G. max* TE annotation in a couple of days. Set up for the runs was fairly simple. I ran the call through a bash script via `tmux` and used `tee` to capture a log file. Since this portion of the pipeline is not necessary for lncRNA ID the scripts will not be included in this repo, but EDTA will be integrated into the IWGC annotation pipeline so that the lncRNA annotation can reference the results.

EDTA requires a genome assembly and a CDS file. General steps for running EDTA in a pipeline are:
1. Genome assembly
2. (Optional) quick provisional gene prediction → bootstrap CDS
3. EDTA (TE library + annotation), using that bootstrap CDS to reduce false positives
4. Soft-mask the genome using EDTA's TE library
5. Real structural/functional gene annotation (BRAKER/MAKER) on the masked genome
6. (Optional) refine TE annotation again using the now-final CDS, if you want a cleaner pass

We'll have to decide how we want to swap out the standing Repeatmodeler/Repeatmasker calls in the IWGC pipeline. Since we are already pulling protein sequences from related species for gene prediction we can also pull their CDS files for use in EDTA. THEY SHOULD NOT BE BACK TRANSLATED CDS FILES FROM THE PROTEIN FILES they need to be directly pulled along with the proten files. EDTA is built to use potential coding sequences to sift out any false positives. Since coding space is under selective pressure, using a CDS file from a closely related species will work fine for this step. Using `datasets` to pull protein and CDS files will make this simple.

---
## General EDTA run notes
```bash
# Pull the EDTA version 1.9.6 image
module purge
module load singularity
singularity pull EDTA_v1.9.6.sif docker://quay.io/biocontainers/edta:1.9.6--hdfd78af_2
# General run command
taskset -c 42-102 singularity exec -B "$INDIR" -B "$OUTDIR" "$SIF" \
  EDTA.pl --genome "$GENOME" \
    --species maize \
    --step all \
    --cds "$CDS" \
    --anno 1 \
    --sensitive 1 \
    --threads 50

# taskset -c 42-102: run on cores 43 through 103
# --species maize: specific run library for maize. Otherwise, use 'others'
# --step all: run the entire pipeline
# --cds: point to the CDS file
# --anno 1: include the gff3 output file
# --sensitive 1: Runs repeatmodeler after EDTA is completed to remove redunancies and clean up annotations
# Note: this is by far the longest, most resource-heavy step. 
# --threads: Number of cores to use.
# Note: a few steps are limited on core count, so core count does not need to be excessive even for large genomes.
```

## Taskset:
Cool tool to subset specific cores to keep the run organized. I started using taskset because version 2 was blowing up cores past the --threads count. The -c flag corresponds to the core number that comes up in htop- note that htop core count starts at 1 and taskset starts at 0.