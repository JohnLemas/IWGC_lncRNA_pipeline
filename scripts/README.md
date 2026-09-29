# Template script directory
### John Lemas
### 5-8-2026
---
This directory contains all scripts used to predict and verify lncRNAs in a given species. Scripts can be copied into a species' sub-directory for editing.
```bash
# 1. Clone the main repo
git clone https://github.com/JohnLemas/IWGC_lncRNA_pipeline.git

# 2. Create the base working directory
mkdir A.thaliana

# 3. Create the pipeline directories
mkdir A.thaliana/01_inputs A.thaliana/02_scripts A.thaliana/03_outputs

# 4. Copy the basic scritps into the working directory
cp IWGC_lncRNA_pipeline/scripts/* A.thaliana/02_scripts

# 5. Insert your species prefix identifier into all of the SLURM scripts
sed -i 's/species/A.thaliana/g' A.thaliana/02_scripts/*sb

# 6. Update the $BASE file path for each script and verify species prefix
```
---
### Contents
- [Alignment array](#alignment-array)
  - Array `SLURM` job that uses an sra list to align trimmed `fastq` files
- [BARRNAP](#barrnap)
  - Builds a custom rRNA database from the input genome
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
- [Intersect-classify](#intersect-classify)
  - Final step in the pipeline
  - Runs R scripts and FEELnc_classify through the rbio_tools and FEELnc images to annotate predictions.
- [lncBoost](#lncboost)
  -  XGBoost model - iterative ML 'gradient boosting' model
- [lncFinder](#lncfinder)
  - support vector machine (SVM) geometric Ml tool - good for finding patterns in nonlinear data
  - image, r script, and sbatch script for this step
- [QC](#qc)
  - Runs `FastQC` and `MultiQC` on raw input files, `Fastp` trimmed files, and `SortMeRNA` cleaned files
- [SortMeRNA](#sortmerna)
  - Filters out rRNA contamination in input raw reads
- [Genome generate](#genome-generate)
  - Indexes the reference genome for STAR read alignment
- [Stringtie](#stringtie)
  - Creates gtf for each fastq file and merges individual gtf files into one *de novo* transcriptome for the species
- [Track Predictions](#track-predictions)
  - Uses metadata to assign tissue type specificity to each predicted lncRNA
  - Can be run directly after StringTie on the entire candidate set
---
### Alignment array
INPUTS
- STAR indexed genome
- trimmed reads
- SRA accession list

COMMAND
```bash
STAR \
    --runThreadN $SLURM_CPUS_PER_TASK \
    --genomeDir $INDEX_DIR \
    --readFilesIn $READS \
    --readFilesCommand zcat \
    --outFileNamePrefix ${OUT_DIR}/${ID}_ \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattrIHstart 0 \
    --outSAMattributes NH HI AS nM NM \
    --outSAMstrandField intronMotif
```
FLAGS
 Flags
- outSAMattrIHstart 0
  - required for compatibility with certain older tools or specific pipelines that expect 0-based indexing for multi-mappers
- outSAMattributes
  - NH: Number of reported alignments for the query
  - HI: Query hit index (distinguishes multiple alignments of the same read)
  - AS: Alignment score
  - nM: Number of mismatches
  - NM: Edit distance to the reference (including insertions and deletions).
- outSAMstrandField
  - intronMotif is used for strand-specific RNA-seq data, which is common in eukaryotic transcriptomics
  - It helps STAR to correctly assign reads to the appropriate strand based on the presence of canonical splice motifs (GT/AG, GC/AG, AT/AC) at intron-exon boundaries
  - This is crucial for accurate gene expression quantification and downstream analyses.

SLURM NOTES
- cpus-per-task=8
- mem=64G
- time=07:00:00
  - based on 160 sra accessions
---
### BARRNAP
INPUTS
- reference genome assembly

COMMAND
```bash
barrnap \
    --kingdom "${kingdom}" \
    --threads "${THREADS}" \
    --outseq "${out_fasta}" \
    "${GENOME}" \
    > "${OUTDIR}/${name}_${kingdom}_rrna.gff"
```
FLAGS
- --kingdom
  - barrnap database selection (euk, bac, mito)
- --threads
  - barrnap runs pretty efficiantly. 8 is enough for most genomes
- --outseq
  - rRNA fasta file with matches to the kingdom database

NOTES
- very fast
  - 8 cpu, 16 G and 30 min for G. max
---
### CPAT plant
Model: Logistic regression
- maps features to a sigmoid curve to output a probability between 0 and 1
- simplest model while still considered ML

INPUTS:
1. Image
2. fasta file of candidate sequences
3. hexamer table from Plant-LncPipe github: 
  - `/base/dir/Plant-LncRNA-pipeline-v2/Model/Plant_Hexamer.tsv`
4. logistic regression R data from Plant-LncPipe github
  - `/base/dir/Plant-LncRNA-pipeline-v2/Model/Plant.logit.RData`
5. output file

COMMAND:
```bash
singularity exec $IMAGE cpat.py -x $HEX -d $LOGIT -g $SEQS -o $OUT_ALL
```

FLAGS:
plain and simple. figure it out.

SLURM NOTES:
- 32 G
- 16 cpus-per-task
- 2 hours for soy
---
### DIAMOND
Double checks that any candidate transcripts do not match any known protein coding transcripts. See main `README.md` for database set-up.

Database path on hpcc: \
`/mnt/scratch/lemasjoh/BASF/uniprot_db/uniprot_out.dmnd`

INPUTS
1. candidate fasta sequences
2. output file path
3. database path

COMMAND
```bash
diamond blastx --threads $SLURM_CPUS_PER_TASK -d $DB \
	-q $INPUT -o $OUTPUT --outfmt 6 --sensitive \
	--max-target-seqs 1 --evalue 1e-5
```
FLAGS
- outfmt 6
  - creates a tsv file instead of a long human readable text file
- sensitive
  - diamond uses a 'seed' method to find matches
  - it is fast but might miss distant homologs
  - this flags tells it to use a more complex seeding strategy to find matches that are more divergent
  - if a transcript is acutally a degraded protein coding gene or other type of artefact this will hopefully catch it
- match-target-seqs 1
  - only reports the single best hit from matched in the db
	reduces output file size
- evalue 1e-5
  - e value: expectation value- is this hit just a random alignment or an actual hit? Lower expectation, lower value
  - currently set to a moderately strict threshold 

SLURM NOTES
- 32 cpus-per-task
- 64 G
- 2 hours for soy (works surprisingly fast- like, minutes)
---
### fastp
Cleans up raw `fastq` files from `prefetch_array.sb`. Trimmed reads are placed in the `species/03_outputs/01_fastp` directory. Double check that the array field in the SLURM parameters correctly reflects the number of accessions, and that the file paths for the species in question are correct.

INPUTS
- Raw `fastq` files from `prefetch_array.sb`
- SRA accession list
- Output directory for trimmed reads

COMMAND
```bash
sbatch fastp.sb
```
FLAGS

NOTES
Updated the original wrapper and template scripts to use an array instead. Easier on SLURM to organize the 25 jobs on its own without me submitting all of them at once. This script identified paired-end or singl end reads and selects fastp parameters accordingly. Works well when pulling lots of random SRA accessions from NCBI without having to filter for paired-end reads. This will be overkill in future runs where I am more selective of input data.

---
### FEELnc
Model: Decision tree ensemble
- builds hundreds of decision trees. each tree looks at a random subset of data and votes on whether a transcript is coding or noncoding
- the final result is the majority vote of all trees

INPUTS:
1. Image path (if using image)
2. Input gtf file: novel transcript gtfs from alignment/stringtie
3. reference annotation file
4. output file designation
5. gitcloned FEELnc directory

COMMAND:
```bash
singularity exec -B /mnt/scratch/lemasjoh $IMAGE \
    env PERL5LIB="${FEELNC_DIR}/lib" \
    perl ${FEELNC_DIR}/scripts/FEELnc_filter.pl \
    -i $INPUT_GTF \
    -a $REF_ANN \
    --monoex=-1 -s 200 -p 20 > $OUTPUT_GTF
```
FLAGS:
- B
  - bind the working directory into the singularity environment; avoids file not found errors
- env
  - makes sure to call the correct perl library; essential for hpcc
- perl
  - the perl script from the gitcloned repo
- monoex=-1
  - monoexonic transcripts are usually noise in rnaseq data. most pipelines just throw them away
  - -1 means keep any monoexonic transcripts that are antisense to an mRNA
  - monoexonic antisense transcripts are much more likely to be a real noncoding transcript than a random monoexonic transcript out in intergenic space
- s 200
  - sets a transcripts size threashold to 200; a common value in lncRNA papers
- p 20
  - parallel threads to use

SLURM NOTES
- 64 G memory
- 20 cpus-per-task
- 4 hours for the 98,984 soy candidates
---
### Intersect-classify
Multiple steps for this script:
1. Intersect IDs
  - Intersects predicted lncRNA IDs from the four tools
  - Calls `intersect_ids.r` in the `rbio_tools` image
2. Pull seqs and gtf
  - Outputs final predicted lncRNAs in `fasta` and `gtf` format
3. Classify lncRNA predictions
  - runs `FEELnc_classifier.pl` through the `FEELnc` image
  - Calls `classifications.r` in the `rbio_tools` image
4. Intersect TEs
  - Intersects predictions with annotated TEs
  - Calls `intersect_lncRNA_TE.r` in the `rbio_tools` image
5. Finalize annotation
  - Builds the tsv file with classifications and TE information
    - Calls `build_final_annotation.r` in the `rbio_tools` image
6. Add tissue type information
  - Adds tissue specificity column to annotation file
  - Calls `add_tissue_annotation.r` in the `rbio_tools` image

INPUTS \
Designed to accept at least two outputs from any of the four tools. Fewer pipeline tools from the entire ensemble can be run if desired, or if a tool/model is not working. 

- FEELnc output
  - accepts either the raw `gtf` file or an indexed list of IDs
- `CPAT` output in `tsv` format
- `lncBoost` output in `tsv` format
- `lncFinder` output in `csv` format
- `DIAMOND` output in `tsv` format
- output file prefix
  - default is `concensus`
- TE annotations
- Tracking summary from `track_preds.sb`


COMMAND
```bash
singularity run $IMAGE \
-f $FEELnc \
-c $CPAT \
-l $lncFinder \
-b $lncBoost \
-d $DIAMOND \
-o $OUTPFX
```
FLAGS

NOTES
4 G mem, 1 hour, 4 threads

---
### lncBoost
Model: XGBoost
- the iterative ML
- 'Gradient boosting' model
  - like random forest, but instead of making them all at once it makes them one by one
  - each new tree is designed to fix the 'errors' of the previous trees
  - currently one of the most powerful 'tabular' Ml tools in existence

INPUTS:
1. Image
2. path to script directory from Plant-LncPipe github
  - `/base/dir/Plant-LncRNA-pipeline-v2/Script`
3. path to model from Plant-LncPipe github
  - `/base/dir/Plant-LncRNA-pipeline-v2/Model/PlantLncBoost_model.cb`
4. fasta file of candidate sequences
5. output directory
  - makes two files. where do you want them
6. output features.csv file
7. output predictions.txt file
  - uses features to make predictions

COMMANDS:
1. feature extraction
```bash
singularity exec -B /mnt/scratch/lemasjoh $IMAGE \
    python3 ${SCRIPT_DIR}/Feature_extraction.py -i $INPUT_FASTA -o $FEATURE_CSV
```

FLAGS:
`${SCRIPT_DIR}` is the path to script directory from Plant-LncPipe github

2. prediction
```bash
singularity exec -B /mnt/scratch/lemasjoh $IMAGE \
    python3 ${SCRIPT_DIR}/PlantLncBoost_prediction.py \
    -i $FEATURE_CSV \
    -m $MODEL \
    -o $FINAL_RESULTS \
    -t 0.5
```
FLAGS:
`${SCRIPT_DIR}`
	path to script directory from Plant-LncPipe github
- i
  - input fasta
- m 
  - path to model from Plant-LncPipe github
- o
  - path to predictions.txt
- t 0.5
  - sets the threshold
  - lncBoost doesnt just 1 or 0 predictions, it assigns probability score between 0 and 1
	  - 0.5 means if 50% probability or higher, call it a lncRNA
  - standard, neutral decision for binary classification
	ok in this case because of the ensemble of Ml tools used in this pipeline
  - dont want one tool to filter out possibilities that other tools thought were lncRNAs

SLURM NOTES
- 20 cpus-per-task
- 64 G
- 4 hours for soy
---
### lncFinder
Model: support vector machine (SVM)
- the geometric ML
- plots transcripts as dots in high-dimensional space
- attempts to find the hyperplane (perfect boundary line) that creates the widest possible gap between the coding dots and noncoding dots
- good at finding patterns in nonlinear data

SIDE NOTE! different than other models. Run it with the rbio_tools image
- The authors of Plant-LncPipe have an R script in the Scripts directory that runs lncFinder
- I didn't like it so Gemini and I rewrote it to work with the custom rbio_tools image I wanted to use
- also, as I eventually discovered through much agony, the model that the Plant-LncPipe authors provide for lncFinder is broken/disfunctional
- luckily, lncFinder has a built in model that was trained on wheat data. use this model instead
- This script call viennaRNA to fold each candidate. All of the dependancies are built into rbio_tools
	
INPUTS:
1. candidate fasta sequences
2. path to output csv file

COMMAND (in the R script):
```bash
singularity run \
    --bind /mnt/scratch/lemasjoh \
    "$IMAGE" lnc_predict \
    --fasta "$INPUT_FASTA" \
    --out "$OUTPUT_CSV" \
    --cores "$SLURM_CPUS_PER_TASK" \
    --model "$MODEL_TYPE"

```

FLAGS:
- --bind
  - bind the path to the project root
- --fasta
  - candidate sequences in fasta format
- --out
  - output directory path with file name
- --cores
  - resource allocation
- --model
  - Passes the selected model type to `lncFinder()`
  - use 'wheat' for all plant species- the only others are human and mouse

NOTES
- Baking the R scripts into rbio_tools makes these runs much easier. Singularity passes these variables directly to the R script eliminating the need to change file paths within the R script itself.
- 32 cpus-per-task
- 64 G
- 12 hours
  - HUGE NOTE: folding cnadidate sequences takes a lot of time
  - tried separating these tasks but `lncFinder()` is really picky about folded input format

---
### QC
INPUTS
- Raw fastq samples from sequencing facility
- Trimmed fastq samples from `fastp.sb`
- Clean fastq output from `sortmerna.sb`

COMMAND
```bash
# General fasp command (not needed for trimmed reads)
ml load FastQC
fastqc -t 8 -o "${OUT_RAW}" "${IN_RAW}"/*.gz

# General multiqc command
ml load MultiQC
multiqc "${OUT_RAW}" -n "raw_report" -o "${OUT_RAW}" --no-data-dir
```
FLAGS
- -n
  - prefix for the report file
- -o
  - output directory
- --no-data-dir
  - skips the default data directory output
NOTES
This step is used to track not only sequence quality but also library prep via ribosomal depletion. SortMeRNA has strict matching thresholds, so read libraries in the `clean_reads` directory will generall drop significantly if ribosomal depletion was not effective.
### SortMeRNA
INPUTS
- `barrnap` custom rDNA library and default SortMeRNA database
- trimmed reads from pastp

COMMAND
```bash
sortmerna \
     --ref ${DEFAULT} \
     --ref ${CUSTOM} \
     --reads "${IN_DIR}/${NEW_ID}_1_trim.fastq.gz" \
     --reads "${IN_DIR}/${NEW_ID}_2_trim.fastq.gz" \
     --workdir "$WORKDIR" \
     --aligned "${RRNA_DIR}/${NEW_ID}_pe_rRNA" \
     --other "${CLEAN_DIR}/${NEW_ID}_pe_non_rRNA" \
     --paired_in \
     --out2 \
     --fastx \
     --threads "$THREADS"
```

FLAGS
- --ref
  - Path to the database for the run. Can use multiple databases via duplicate calls
- --reads
  - Path to trimmed fastq files
- --workdir
  - Path for SortMeRNA temporary working directories. One directory is generated for each sample.
- --aligned
  - File prefix for reads that match the rRNA database(s)
- --other
  - File prefix for reads that did not match the database(s)
- --paired_in
  - Include for paired end reads
- --out2
  - maintain paired-end read format in output files instead of collapsing them
- --fastx
  - output the same file format as input (fasta or fastq) instead of alignment files
- --threads
  - parallel cores to use
NOTES
- This is the longest step in the pipeline next to `lncFinder`
- rRNA contamination is the biggest hurdle we have to navigate with lncRNA prediction
- The combination of barrnap and sortmeRNA should handle most of the rRNA filtering, but downstream predictions should be scrutinized accordingly
### Genome generate
INPUTS
- Reference genome
  - `species/01_inputs/species.chromosomes.fa`
- Reference annotation 
  - `species/01_inputs/species_mRNA.gtf
- Index directory path

COMMAND
```bash
STAR \
     --runMode genomeGenerate \
     --runThreadN $SLURM_CPUS_PER_TASK \
     --genomeDir $INDEX \
     --genomeFastaFiles $GENOME \
     --sjdbGTFfile $ANNO \
     --sjdbGTFtagExonParentTranscript Parent \
     --sjdbOverhang 99 \
     --genomeSAindexNbases 13
```
FLAGS
- -runMode genomeGenerate
- -sjdbGTFtagExonParentTranscript Parent
- -sjdbOverhang 99
- -genomeSAindexNbases 13
NOTES

---
### Stringtie
stringtie transcriptome assembly and lncRNA candidate discovery

1. `stringtie_assemble.sb` \
INPUTS
- metadata file
- STAR bam files
- reference `GTF` 

COMMAND
```bash
stringtie \
    ${BAM_DIR}/${ID}_Aligned.sortedByCoord.out.bam \
    -G $REF_GTF \
    -p $SLURM_CPUS_PER_TASK \
    -m 200 \
    -j 2 \
    -o ${OUT_DIR}/${ID}_assembled.gtf"
```
FLAGS
- -G
  - reference GTF file
- -p
  - threads
- -m
  - minimum length threshold
- -j
  - Minimum junction coverage (higher = fewer false positives)
- -o
  - output file
2. `stringtie_merge.sb` \
INPUTS
- species metadata file
- assembled gtf files
- reference `GTF` and `FASTA` files
- output file paths
COMMAND
```bash
# Run the merge command
stringtie --merge -G $REF_GTF -m 200 -o $OUT_GTF ${gtf_filePath[*]}

# Extract the transcript sequences from the master GTF file
grep 'transcript_id "MSTRG' $OUT_GTF > $CAN_GTF
gffread -w $SEQ_OUT -g $GENOME $CAN_GTF
```
NOTES
- assembly runs as an array job using the metadata to build the array
- merge combines *de novo* transcriptome assembly and candidate feature extraction
- outputs candidate lncRNA files used for all prediction tools
---