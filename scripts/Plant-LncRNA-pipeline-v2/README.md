# Plant-LncRNA-pipeline-v2
A pipeline for identifying and characterizing lncRNAs in plants.

![LncRNA pipline-V2](https://github.com/xuechantian/Plant-LncRNA-pipeline-v2/blob/master/Plant-LncPipeV2.workflow.png) 



**Contact Information:**  
- Email: xuechan.tian@bjfu.edu.cn  
- Email: jianfeng.mao@umu.se


## Table of Contents

- [1. Introduction](#1-Introduction)
- [2. Dependencies](#2-Dependencies)
- [3. Installation](#3-Installation)
  - [3.1. Install Plant-LncRNA-pipeline-v2](#31-install-Plant-LncRNA-pipeline-v2)
  - [3.2. Install HISAT2](#32-install-hisat2)
  - [3.3. Install StringTie](#33-install-stringtie)
  - [3.4. Install CPAT and LncFinder](#34-install-cpat-and-lncfinder)
  - [3.5. Install Diamond](#35-install-diamond)
  - [3.6. Install FEELnc](#36-install-feelnc)
  - [3.7. Install fastp](#37-install-fastp)
- [4. Download Example Data](#4-download-example-data)
- [5. Data Processing](#5-data-processing)
  - [5.1. Convert SRA to FASTQ](#51-convert-sra-to-fastq)
  - [5.2. Data Filtering and Quality Control](#52-data-filtering-and-quality-control)
- [6. Run HISAT2 for RNA-seq Mapping](#6-run-hisat2-for-rna-seq-mapping)
- [7. Transcript Assembly with StringTie](#7-transcript-assembly-with-stringtie)
- [8. LncRNA Identification](#8-lncrna-identification)
  - [8.1. Remove transcripts shorter than 200 bp and overlapping with known mRNAs](#81-remove-short-transcripts)
  - [8.2. Identification Using PlantLncBoost](#82-identification-using-plantlncboost)
  - [8.3. Identification Using LncFinder-plant](#83-identification-using-lncfinder-plant)
  - [8.4. Identification Using CPAT-plant](#84-identification-using-cpat-plant)
  - [8.5. Alignment with UniProt using Diamond](#85-alignment-with-uniprot-using-diamond)
- [9. Classifying LncRNAs Based on Genomic Location](#9-classifying-lncrnas-based-on-genomic-location)
- [10. TE-derived LncRNAs](#10-te-derived-lncrnas)
- [11. Citations](#11-citations)



---

## 1. Introduction

Plant-LncPipe-v2 is a comprehensive pipeline designed to identify and characterize long non-coding RNAs (lncRNAs) in plants. It integrates multiple tools for RNA-seq alignment, transcript assembly, and lncRNA prediction, providing an efficient and reproducible workflow for plant lncRNA studies.

---

### Updates in v2

Plant-LncPipe-v2 incorporates the following major improvements over the previous version:

- Integration of PlantLncBoost, our new CatBoost-based machine learning model specifically optimized for plant lncRNA identification
- Enhanced prediction accuracy across diverse plant species

### Workflow

Plant-LncPipe-v2 consists of four main modules:

1. **Raw Data Preprocessing**: Quality control, adapter trimming, and sequence filtering
2. **Transcript Assembly**: Reference-guided assembly of transcripts from RNA-seq data
3. **lncRNA Identification**: Prediction using PlantLncBoost, CPAT-plant and LncFinder-plant
4. **lncRNA Classification and Characterization**: Categorization and analysis of identified lncRNAs


## 2. Dependencies



    HISAT2
    StringTie
    LncFinder
    CPAT
    Diamond




## 3. Installation

### 3.1. Install Plant-LncRNA-pipeline-v2
    git clone https://github.com/xuechantian/Plant-LncRNA-pipeline-v2.git
    
### 3.2. Install HISAT2
    git clone https://github.com/DaehwanKimLab/hisat2.git
    cd hisat2
    make

### 3.3. Install StringTie
    wget https://github.com/gpertea/stringtie/releases/download/v2.1.4/stringtie-2.1.4.Linux_x86_64.tar.gz
    tar xzf stringtie-2.1.4.Linux_x86_64.tar.gz
    export PATH=$PATH:/path/stringtie-2.1.6.Linux_x86_64/
	
### 3.4. Install CPAT and LncFinder

CPAT (version 1.2.4)


    conda create -n py27 python=2.7 -y
    source activate py27
    pip2 install CPAT


 
LncFinder-plant R Package

 
    install.packages("LncFinder")
    install.package("seqinr")
	
	
### 3.5. Install Diamond
     conda install -c bioconda diamond
          

### 3.6. Install FEELnc
    git clone https://github.com/tderrien/FEELnc.git
    export FEELNCPATH=/path/FEELnc/bin/
    export PERL5LIB=$PERL5LIB:/path/FEELnc/lib/
    export PATH=$PATH:/path/FEELnc/scripts/
	
	
### 3.7. Install fastp
    conda install -c bioconda fastp	
	
	
	
## 4. Download Example Data



### Genome Assembly and Annotation File
    Download the genome assembly fasta file and gff file from Phytozome 13 for Glycine max Wm82.a6.v1


### Download Transcriptome Sequencing Data
Example content of sra.list file


    SRR1174214
    SRR1174217
    SRR1174218
    SRR1174218
    SRR1174232
    ......


 
Download


    prefetch --option-file sra.list
	
	

 
## 5. Data Processing

### 5.1. Convert SRA to FASTQ
    fastq-dump SRR1174214.sra
    fastq-dump SRR1174217.sra
    fastq-dump SRR1174218.sra
    fastq-dump SRR1174228.sra
    fastq-dump SRR1174232.sra

	
	
### 5.2. Data Filtering and Quality Control
    fastp -i SRR1174214.fastq -o SRR1174214_clean.fastq 
    fastp -i SRR1174217.fastq -o SRR1174217_clean.fastq 
    fastp -i SRR1174218.fastq -o SRR1174218_clean.fastq 
    fastp -i SRR1174228.fastq -o SRR1174228_clean.fastq 
    fastp -i SRR1174232.fastq -o SRR1174232_clean.fastq 


## 6. Run HISAT2 for RNA-seq Mapping



### Construct reference genome

    hisat2-build -p 8 genome.fasta genome.index 

	

### Genome alignment with hisat2
#### Single-End RNA-seq  (e.g., Soybean)


If the RNA-seq library is strand-specific, add the parameter "--rna-strandness RF".


    for i in `cat sra.list`; do hisat2 --new-summary --rna-strandness RF -p 10 -x genome.index ${i}_clean.fastq -S ${i}.sam; done

 
If the RNA-seq library is not strand-specific, remove  the parameter "--rna-strandness RF".


    for i in `cat sra.list`; do hisat2 --new-summary -p 10 -x genome.index ${i}_clean.fastq -S ${i}.sam; done

	
#### Paired-End RNA-seq


 If the RNA-seq library is strand-specific, add the parameter "--rna-strandness RF".

 
    for i in `cat sra.list`; do hisat2 --new-summary --rna-strandness RF -p 10 -x genome.index -1 ${i}_1_clean.fastq -2 ${i}_2_clean.fastq -S${i}.sam; done
	
 If the RNA-seq library is not strand-specific, remove  the parameter "--rna-strandness RF".

 
    for i in `cat sra.list`; do hisat2 --new-summary -p 10 -x genome.index -1 ${i}_1_clean.fastq -2 ${i}_2_clean.fastq -S ${i}.sam; done

	
	
### Sort and compress sam files with samtools

    for i in `cat sra.list`; do samtools view -S -b ${i}.sam | samtools sort -o ${i}.bam; done
	

	

## 7. Transcript Assembly with StringTie

### Format of Glycine_max_longest.gtf
    Gm01	phytozomev13	exon 	103572	103594	transcript_id "Glyma.01G000100.2"; gene_id "Glyma.01G000100.2" 
    Gm01	phytozomev13	CDS 	103572	103594	transcript_id "Glyma.01G000100.2"; gene_id "Glyma.01G000100.2"
    Gm01	phytozomev13	exon 	103222	103288	transcript_id "Glyma.01G000100.2"; gene_id "Glyma.01G000100.2"	
    Gm01	phytozomev13	CDS 	103222	103288	transcript_id "Glyma.01G000100.2"; gene_id "Glyma.01G000100.2"	
    Gm01	phytozomev13	exon	102790	102857	transcript_id "Glyma.01G000100.2"; gene_id "Glyma.01G000100.2"	
    Gm01	phytozomev13	CDS	102790	102857	transcript_id "Glyma.01G000100.2"; gene_id "Glyma.01G000100.2";	
    Gm01	phytozomev13	exon	78986	79111	transcript_id "Glyma.01G000100.2"; gene_id "Glyma.01G000100.2";	
    Gm01	phytozomev13	CDS	78986	79111	transcript_id "Glyma.01G000100.2"; gene_id "Glyma.01G000100.2";	

	
	
### Transcription reconstruction of a single sample

If strand-specific, add the parameter "--rf".


	for i in `cat sra.list`; do stringtie -p 10 --rf -G Glycine_max_longest.gtf  -o ${i}.gtf  ${i}.bam; done

 
If not strand-specific, remove the parameter "--rf".


    for i in `cat sra.list`; do stringtie  -p 10 -G Glycine_max_longest.gtf -o ${i}.gtf  ${i}.bam; done
	
	
### Merge transcripts from multiple samples
	stringtie --merge -o merge.gtf  -G Glycine_max_longest.gtf  SRR*.gtf
	grep 'transcript_id "MSTRG' merge.gtf > candidate_transcript.gtf
	gffread -w candidate_transcript.fasta -g genome.fasta candidate_transcript.gtf
	grep '>' candidate_transcript.fasta | awk '{print \$1}' | sed 's/>//g' | sort -u > candidate_transcript.txt

	
	
## 8. LncRNA Identification



### 8.1. Remove transcripts shorter than 200 bp and overlapping with known mRNAs

    FEELnc_filter.pl -i candidate_transcript.gtf -a Glycine_max_longest.gtf --monoex=-1 -s 200 -p 20 > candidate_lncRNA.gtf
    cut -d ";" -f 2 candidate_lncRNA.gtf |sed 's/ transcript_id //g' | sed 's/"//g' | sort -u > candidate_lncRNA.txt


	
### 8.2. Identification Using PlantLncBoost
### Dependency 
    Python (>=3.7)
    Biopython
    NumPy
    Pandas
    SciPy
    CatBoost

### Install via conda
#### Create and activate a conda environment
    conda create -n lncrna_env python=3.7
    conda activate lncrna_env
#### Install core dependencies
    conda install -c conda-forge biopython numpy pandas scipy catboost


### Run PlantLncBoost

#### Feature extraction
    python PlantLncBoost/Script/Feature_extraction.py -i candidate_transcript.fasta -o PlantLncBoost_feature.csv
	
#### LncRNA prediction
In the second column (Predicted_label) of the result file, 1 represents lncRNA and 0 represents mRNA.

    python PlantLncBoost/Script/PlantLncBoost_prediction.py -i PlantLncBoost_feature.csv -m PlantLncBoost/Model/PlantLncBoost_model.cb -t 0.5 -o PlantLncBoost_prediction.csv
    
	
### 8.3.  Identification Using LncFinder-plant

R Package

    R
    library(LncFinder)
    library(seqinr)
	
	
import training data

    mRNA <- seqinr::read.fasta(file ="./data/training/mRNA.fasta")
    lncRNA <- seqinr::read.fasta(file ="./data/training/lncRNA.fasta")
    
    
Use "make_frequencies" function to generate the feature file

    frequencies <- make_frequencies(cds.seq = mRNA, lncRNA.seq = lncRNA, SS.features = FALSE, cds.format = "DNA", lnc.format = "DNA", check.cds = TRUE, ignore.illegal = TRUE)	
	
	
loading the model

    plant = readRDS("./Model/Plant_model.rda")
	
	
Identification of lncRNA 

    Seqs <- seqinr::read.fasta(file ="candidate_transcript.fasta")
    Plant_results <- LncFinder::lnc_finder(Seqs, SS.features = FALSE, format = "DNA", frequencies.file = frequencies, svm.model = plant, parallel.cores = 2)
	
	
Export results

    write.table(Plant_results, file ="plant-lncFinder.txt", sep ="\t",row.names =TRUE, col.names =TRUE,quote =FALSE)




### 8.4. Identification Using CPAT-plant
The coding probability (CP) cutoff: 0.46 (CP >=0.46 indicates coding sequence, CP < 0.46 indicates noncoding sequence).

    source activate py27
    cpat.py -x ./Model/Plant_Hexamer.tsv -d ./Model/Plant.logit.RData -g candidate_transcript.fasta -o CPAT_plant.output





### 8.5. Alignment with UniProt using Diamond
    wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
    gunzip uniprot_sprot.fasta.gz
    diamond makedb --in uniprot_sprot.fasta -d uniprot_out
    diamond blastx -d uniprot_out -q candidate_transcript.fasta -o uniprotoutput.txt




### 8.6. By intersecting the results obtained from the aforementioned steps, a set of high-confidence lncRNAs were obtained

The id of the lncRNA

    Rscript prediction_insersection.sh candidate_lncRNA.txt PlantLncBoost_prediction.csv CPAT_plant.output plant-lncFinder.txt uniprotoutput.txt
    
The lncRNA gtf file

    grep -Fwf Final_lncRNA_results.txt candidate_transcript.gtf > lncRNA.gtf

	

  
	
## 9. Classifying LncRNAs Based on Genomic Location
Classification result file

	FEELnc_classifier.pl -i lncRNA.gtf -a Glycine_max_longest.gtf > lncRNA_classes.txt
	
	
Antisense_exonic-lncRNA

	awk -F '\t' '{if($1==1 && $6 == "antisense" && $10 == "exonic") {print $0}}' lncRNA_classes.txt > LncRNA_antisense_exonic.txt


Intronic-lncRNA

	awk -F '\t' '{if($1==1 && $6 == "sense" && $10 == "intronic") {print $0}}' lncRNA_classes.txt > LncRNA_intronic.txt


Upstream-lncRNA

	awk -F '\t' '{if($1==1 && $6 == "sense" && $7 == "intergenic" && $8 <= 2000 && $10 == "upstream") {print $0}}' lncRNA_classes.txt > LncRNA_upstream.txt
	
	
Downstream-lncRNA

	awk -F '\t' '{if($1==1 && $6 == "sense" && $7 == "intergenic" && $8 <= 2000 && $10 == "downstream") {print $0}}' lncRNA_classes.txt > LncRNA_downstream.txt


Intergenic-lncRNA

	awk -F '\t' '{if($1==1 && $7 == "intergenic" && $8 > 2000) {print $0}}' lncRNA_classes.txt > LncRNA_intergenic.txt


Bidirectional-lncRNA

	awk -F '\t' '{if($1==1 && $6 == "antisense" && $7 == "intergenic" && $8 <= 2000 && $10 == "upstream") {print $0}}' lncRNA_classes.txt > LncRNA_Bidirectional.txt
	
	
	
	
## 10. TE-derived LncRNAs
cat TE.bed

	Chr1    15827287        15838845        LTR_Gypsy
	Chr1    13085455        13085593        LTR_Copia
	Chr1    11181821        11181959        LTR_Copia
	Chr1    20699111        20699248        LTR_Copia
	...
	
	
cat LncRNA.bed

	Chr1    11171031        11171031        MSTRG.2781.1
	Chr1    12199350        12199350        MSTRG.2973.1
	Chr1    13466928        13466928        MSTRG.3115.1
	Chr1    13838536        13838536        MSTRG.3127.1
	...
	
	
TE-lncRNA

	bedtools intersect -a lncRNA.bed -b TE.bed -wo | sort -u | awk '{print $1,$2,$3,$4,$6,$7,$8,$9}' | sed 's/ /\t/g' | sed '1iChr\tLncRNA_start\tLncRNA_end\tLncRNA_ID\tTE_start\tTE_end\tTE_ID\tOverlap' > TE_lncRNA_intersect.txt 
	
	
## 11. Citations
If you use Plant-LncPipe-V2, please cite:

Tian, X.-C., Nie, S., Domingues, D., Rossi Paschoal, A., Jiang, L.-B. and Mao, J.-F. (2025), PlantLncBoost: key features for plant lncRNA identification and significant improvement in accuracy and generalization. New Phytol. https://doi.org/10.1111/nph.70211


	
	
	
	
