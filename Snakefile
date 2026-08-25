"""
lncRNA annotation pipeline (Snakemake port)
============================================

All species/project-specific settings live in config.yaml. This Snakefile
should not need editing between runs/species.

Typical usage
-------------
Dry run (check the DAG without executing anything):
    snakemake -n --use-envmodules

Local run:
    snakemake --use-envmodules -j <cores>

Cluster (SLURM) run, using the Snakemake SLURM executor plugin:
    pip install snakemake-executor-plugin-slurm --break-system-packages
    snakemake --use-envmodules --executor slurm -j <max_concurrent_jobs>

Steps that use `envmodules:` need --use-envmodules on the command line, or
they silently run with whatever is already on PATH.
"""

import os
import re

configfile: "config.yaml"

# ---------------------------------------------------------------------------
# Resolve paths from config
# ---------------------------------------------------------------------------
SPECIES = config["species"]
BASE = config["base_dir"].rstrip("/")
PROJECT = f"{BASE}/{SPECIES}"


def resolve_project(path_template):
    """Resolve a path relative to {base_dir}/{species}, substituting {species}."""
    p = path_template.format(species=SPECIES)
    return p if os.path.isabs(p) else f"{PROJECT}/{p}"


def resolve_base(path_template):
    """Resolve a path relative to {base_dir} (used for shared tools/containers)."""
    p = path_template.format(species=SPECIES)
    return p if os.path.isabs(p) else f"{BASE}/{p}"


SAMPLES_FILE = resolve_project(config["samples_file"])
GENOME_FASTA = resolve_project(config["reference"]["genome_fasta"])
ANNOT_FORMAT = config["reference"]["annotation_format"]
ANNOT_RAW = resolve_project(config["reference"]["annotation_raw"])
ANNOT_GTF = resolve_project(config["reference"]["annotation_gtf"])
PROTEIN_FASTA = resolve_project(config["reference"]["protein_fasta"])
RAW_FASTQ_DIR = resolve_project(config["raw_fastq_dir"])
UNIPROT_DMND = config["uniprot_diamond_db"]
BIND_DIR = config["bind_dir"]
LOGS_DIR = config["logs_dir"]

DIRS = {k: resolve_project(v) for k, v in config["dirs"].items()}

CONTAINERS = {k: resolve_base(v) for k, v in config["containers"].items()}
TOOLS = {k: resolve_base(v) for k, v in config["tool_paths"].items()}

MOD = config["envmodules"]
P = config["params"]

# The reference GTF actually used downstream: produced by the conversion
# rule if the input annotation is GFF3, otherwise the GTF the user supplied.
REF_GTF = ANNOT_GTF

STAR_INDEX_DIR = resolve_project("01_inputs/{species}_indexed".format(species=SPECIES))
SINGULARITY_BIND = f"-B {BIND_DIR}"


# ---------------------------------------------------------------------------
# Sample sheet + paired/single-end detection
# ---------------------------------------------------------------------------
with open(SAMPLES_FILE) as fh:
    SAMPLES = [line.strip() for line in fh if line.strip()]

if not SAMPLES:
    raise ValueError(f"No samples found in {SAMPLES_FILE}")


def _detect_layout(sample):
    pe1 = f"{RAW_FASTQ_DIR}/{sample}_1.fastq.gz"
    pe2 = f"{RAW_FASTQ_DIR}/{sample}_2.fastq.gz"
    se_plain = f"{RAW_FASTQ_DIR}/{sample}.fastq.gz"
    if os.path.exists(pe1) and os.path.exists(pe2):
        return "pe"
    elif os.path.exists(pe1):
        return "se1"  # single-end, but named with a _1 suffix
    elif os.path.exists(se_plain):
        return "se"
    else:
        raise ValueError(
            f"ERROR: no raw fastq files found for sample '{sample}' in {RAW_FASTQ_DIR} "
            f"(expected {sample}_1.fastq.gz [+ _2], or {sample}.fastq.gz)"
        )


SAMPLE_LAYOUT = {s: _detect_layout(s) for s in SAMPLES}
PE_SAMPLES = [s for s, l in SAMPLE_LAYOUT.items() if l == "pe"]
SE_SAMPLES = [s for s, l in SAMPLE_LAYOUT.items() if l in ("se", "se1")]

wildcard_constraints:
    sample = "|".join(re.escape(s) for s in SAMPLES) if SAMPLES else "NEVERMATCH_NO_SAMPLES"


def se_raw_fastq(wildcards):
    layout = SAMPLE_LAYOUT[wildcards.sample]
    if layout == "se1":
        return f"{RAW_FASTQ_DIR}/{wildcards.sample}_1.fastq.gz"
    return f"{RAW_FASTQ_DIR}/{wildcards.sample}.fastq.gz"


def star_reads(wildcards):
    if SAMPLE_LAYOUT[wildcards.sample] == "pe":
        return [
            f"{DIRS['fastp']}/{wildcards.sample}_1_trim.fastq.gz",
            f"{DIRS['fastp']}/{wildcards.sample}_2_trim.fastq.gz",
        ]
    return [f"{DIRS['fastp']}/{wildcards.sample}_trim.fastq.gz"]


def diamond_hits_to_combine(wildcards):
    files = [f"{DIRS['diamond']}/{SPECIES}_diamond_uniprot.tsv"]
    if P["diamond"]["use_species_db"]:
        files.append(f"{DIRS['diamond']}/{SPECIES}_diamond_species.tsv")
    return files


# ---------------------------------------------------------------------------
# Final targets
# ---------------------------------------------------------------------------
rule all:
    input:
        f"{DIRS['intersect']}/{SPECIES}_final_ids.txt",
        f"{DIRS['intersect']}/{SPECIES}_intersection_summary.csv",
        f"{DIRS['intersect']}/{SPECIES}_lncRNA.fasta",
        f"{DIRS['intersect']}/{SPECIES}_lncRNA.gtf",
        f"{DIRS['classify']}/{SPECIES}_lncRNA_classes.txt",
        f"{DIRS['classify']}/.{SPECIES}_classification_parse.done",


# ---------------------------------------------------------------------------
# 0. Reference prep: GFF3 -> GTF (only needed if annotation is GFF3)
# ---------------------------------------------------------------------------
if ANNOT_FORMAT == "gff3":

    rule convert_annotation_to_gtf:
        input:
            gff3=ANNOT_RAW,
        output:
            gtf=ANNOT_GTF,
        envmodules:
            MOD["gffread"],
        log:
            f"{LOGS_DIR}/convert_annotation/{SPECIES}.log",
        shell:
            "gffread {input.gff3} -T -o {output.gtf} > {log} 2>&1"


# ---------------------------------------------------------------------------
# 1. fastp trimming (per sample, paired- or single-end)
# ---------------------------------------------------------------------------
rule fastp_pe:
    input:
        r1=f"{RAW_FASTQ_DIR}/{{sample}}_1.fastq.gz",
        r2=f"{RAW_FASTQ_DIR}/{{sample}}_2.fastq.gz",
    output:
        r1=f"{DIRS['fastp']}/{{sample}}_1_trim.fastq.gz",
        r2=f"{DIRS['fastp']}/{{sample}}_2_trim.fastq.gz",
        html=f"{DIRS['fastp']}/{{sample}}_fastp.html",
        json=f"{DIRS['fastp']}/{{sample}}_fastp.json",
    wildcard_constraints:
        sample="|".join(re.escape(s) for s in PE_SAMPLES) if PE_SAMPLES else "NEVERMATCH_NO_PE",
    threads: P["fastp"]["threads"]
    envmodules:
        MOD["fastp"],
    log:
        f"{LOGS_DIR}/fastp/{{sample}}.log",
    shell:
        "conda run -n {MOD[fastp_env]} fastp "
        "-i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} "
        "--detect_adapter_for_pe --thread {threads} -x -g "
        "--html {output.html} --json {output.json} > {log} 2>&1"


rule fastp_se:
    input:
        r1=se_raw_fastq,
    output:
        r1=f"{DIRS['fastp']}/{{sample}}_trim.fastq.gz",
        html=f"{DIRS['fastp']}/{{sample}}_fastp.html",
        json=f"{DIRS['fastp']}/{{sample}}_fastp.json",
    wildcard_constraints:
        sample="|".join(re.escape(s) for s in SE_SAMPLES) if SE_SAMPLES else "NEVERMATCH_NO_SE",
    threads: P["fastp"]["threads"]
    envmodules:
        MOD["fastp"],
    log:
        f"{LOGS_DIR}/fastp/{{sample}}.log",
    shell:
        "conda run -n {MOD[fastp_env]} fastp "
        "-i {input.r1} -o {output.r1} --thread {threads} -x -g "
        "--html {output.html} --json {output.json} > {log} 2>&1"


# ---------------------------------------------------------------------------
# 2. STAR genome index (once)
# ---------------------------------------------------------------------------
rule star_index:
    input:
        genome=GENOME_FASTA,
        gtf=REF_GTF,
    output:
        directory(STAR_INDEX_DIR),
    threads: P["star_index"]["threads"]
    envmodules:
        MOD["star"],
    log:
        f"{LOGS_DIR}/star_index/{SPECIES}.log",
    params:
        overhang=P["star_index"]["sjdbOverhang"],
        sa_index=P["star_index"]["genomeSAindexNbases"],
    shell:
        "mkdir -p {output} && "
        "STAR --runMode genomeGenerate --runThreadN {threads} "
        "--genomeDir {output} --genomeFastaFiles {input.genome} "
        "--sjdbGTFfile {input.gtf} --sjdbGTFtagExonParentTranscript Parent "
        "--sjdbOverhang {params.overhang} --genomeSAindexNbases {params.sa_index} "
        "> {log} 2>&1"


# ---------------------------------------------------------------------------
# 3. STAR alignment (per sample)
# ---------------------------------------------------------------------------
rule star_align:
    input:
        reads=star_reads,
        index=STAR_INDEX_DIR,
    output:
        bam=f"{DIRS['star']}/{{sample}}_Aligned.sortedByCoord.out.bam",
    params:
        reads_str=lambda wc, input: " ".join(input.reads),
        prefix=f"{DIRS['star']}/{{sample}}_",
    threads: P["star_align"]["threads"]
    envmodules:
        MOD["star"],
    log:
        f"{LOGS_DIR}/star_align/{{sample}}.log",
    shell:
        "STAR --runThreadN {threads} --genomeDir {input.index} "
        "--readFilesIn {params.reads_str} --readFilesCommand zcat "
        "--outFileNamePrefix {params.prefix} --outSAMtype BAM SortedByCoordinate "
        "--outSAMattrIHstart 0 --outSAMattributes NH HI AS nM NM "
        "--outSAMstrandField intronMotif > {log} 2>&1"


# ---------------------------------------------------------------------------
# 4. BAM indexing (per sample)
# ---------------------------------------------------------------------------
rule bam_index:
    input:
        bam=rules.star_align.output.bam,
    output:
        bai=f"{DIRS['star']}/{{sample}}_Aligned.sortedByCoord.out.bam.bai",
    threads: P["bam_index"]["threads"]
    envmodules:
        MOD["samtools"],
    log:
        f"{LOGS_DIR}/bam_index/{{sample}}.log",
    shell:
        "samtools index {input.bam} > {log} 2>&1"


# ---------------------------------------------------------------------------
# 5. StringTie per-sample assembly
# ---------------------------------------------------------------------------
rule stringtie_assemble:
    input:
        bam=rules.star_align.output.bam,
        bai=rules.bam_index.output.bai,
        ref_gtf=REF_GTF,
    output:
        gtf=f"{DIRS['stringtie_asm']}/{{sample}}_assembled.gtf",
    threads: P["stringtie_assemble"]["threads"]
    envmodules:
        MOD["stringtie"],
    log:
        f"{LOGS_DIR}/stringtie_assemble/{{sample}}.log",
    params:
        min_len=P["stringtie_assemble"]["min_transcript_len"],
        min_junction_cov=P["stringtie_assemble"]["min_junction_cov"],
    shell:
        "stringtie {input.bam} -G {input.ref_gtf} -p {threads} "
        "-m {params.min_len} -j {params.min_junction_cov} "
        "-o {output.gtf} > {log} 2>&1"


# ---------------------------------------------------------------------------
# 6. StringTie merge -> candidate GTF -> candidate FASTA
#    (split into three rules for better caching/reruns)
# ---------------------------------------------------------------------------
rule stringtie_merge:
    input:
        gtfs=expand(f"{DIRS['stringtie_asm']}/{{sample}}_assembled.gtf", sample=SAMPLES),
        ref_gtf=REF_GTF,
    output:
        merged=f"{DIRS['stringtie']}/{SPECIES}_denovo_trans.gtf",
    threads: P["stringtie_merge"]["threads"]
    envmodules:
        MOD["stringtie"],
    log:
        f"{LOGS_DIR}/stringtie_merge/{SPECIES}.log",
    params:
        min_len=P["stringtie_merge"]["min_transcript_len"],
    shell:
        "stringtie --merge -G {input.ref_gtf} -m {params.min_len} "
        "-o {output.merged} {input.gtfs} > {log} 2>&1"


rule extract_candidates_gtf:
    input:
        merged=rules.stringtie_merge.output.merged,
    output:
        candidates_gtf=f"{DIRS['stringtie']}/{SPECIES}_candidates.gtf",
    shell:
        """grep 'transcript_id "MSTRG' {input.merged} > {output.candidates_gtf}"""


rule gffread_extract_fasta:
    input:
        candidates_gtf=rules.extract_candidates_gtf.output.candidates_gtf,
        genome=GENOME_FASTA,
    output:
        candidates_fasta=f"{DIRS['gffread']}/{SPECIES}_candidates.fasta",
    envmodules:
        MOD["gffread"],
    log:
        f"{LOGS_DIR}/gffread/{SPECIES}.log",
    shell:
        "gffread -w {output.candidates_fasta} -g {input.genome} "
        "{input.candidates_gtf} > {log} 2>&1"


# ---------------------------------------------------------------------------
# 7. Coding-potential / homology classifiers (all consume the candidate
#    FASTA/GTF and run independently of each other)
# ---------------------------------------------------------------------------
rule feelnc_filter:
    input:
        candidates_gtf=rules.extract_candidates_gtf.output.candidates_gtf,
        ref_gtf=REF_GTF,
    output:
        gtf=f"{DIRS['feelnc']}/{SPECIES}_FEELnc.gtf",
    threads: P["feelnc_filter"]["threads"]
    log:
        f"{LOGS_DIR}/feelnc_filter/{SPECIES}.log",
    params:
        image=CONTAINERS["feelnc"],
        bind=SINGULARITY_BIND,
        feelnc_dir=TOOLS["feelnc_dir"],
        monoex=P["feelnc_filter"]["monoex"],
        min_size=P["feelnc_filter"]["min_size"],
    shell:
        "singularity exec {params.bind} {params.image} "
        "env PERL5LIB={params.feelnc_dir}/lib "
        "perl {params.feelnc_dir}/scripts/FEELnc_filter.pl "
        "-i {input.candidates_gtf} -a {input.ref_gtf} "
        "--monoex={params.monoex} -s {params.min_size} -p {threads} "
        "> {output.gtf} 2> {log}"


rule cpat:
    input:
        seqs=rules.gffread_extract_fasta.output.candidates_fasta,
    output:
        out=f"{DIRS['cpat']}/{SPECIES}_cpat.output",
    threads: P["cpat"]["threads"]
    log:
        f"{LOGS_DIR}/cpat/{SPECIES}.log",
    params:
        image=CONTAINERS["cpat"],
        hexamer=resolve_base(P["cpat"]["hexamer_table"]),
        logit=resolve_base(P["cpat"]["logit_model"]),
    shell:
        "singularity exec {params.image} cpat.py "
        "-x {params.hexamer} -d {params.logit} -g {input.seqs} "
        "-o {output.out} > {log} 2>&1"


rule lncboost_features:
    input:
        seqs=rules.gffread_extract_fasta.output.candidates_fasta,
    output:
        features=f"{DIRS['lncboost']}/{SPECIES}_features.csv",
    threads: P["lncboost"]["threads"]
    log:
        f"{LOGS_DIR}/lncboost/{SPECIES}_features.log",
    params:
        image=CONTAINERS["lncboost"],
        bind=SINGULARITY_BIND,
        script_dir=resolve_base(P["lncboost"]["script_dir"]),
    shell:
        "singularity exec {params.bind} {params.image} "
        "python3 {params.script_dir}/Feature_extraction.py "
        "-i {input.seqs} -o {output.features} > {log} 2>&1"


rule lncboost_predict:
    input:
        features=rules.lncboost_features.output.features,
    output:
        predictions=f"{DIRS['lncboost']}/{SPECIES}_predictions.txt",
    threads: P["lncboost"]["threads"]
    log:
        f"{LOGS_DIR}/lncboost/{SPECIES}_predict.log",
    params:
        image=CONTAINERS["lncboost"],
        bind=SINGULARITY_BIND,
        script_dir=resolve_base(P["lncboost"]["script_dir"]),
        model=resolve_base(P["lncboost"]["model"]),
        threshold=P["lncboost"]["threshold"],
    shell:
        "singularity exec {params.bind} {params.image} "
        "python3 {params.script_dir}/PlantLncBoost_prediction.py "
        "-i {input.features} -m {params.model} -o {output.predictions} "
        "-t {params.threshold} > {log} 2>&1"


rule lncfinder:
    input:
        seqs=rules.gffread_extract_fasta.output.candidates_fasta,
    output:
        csv=f"{DIRS['lncfinder']}/{SPECIES}_lncFinder_results.csv",
    threads: P["lncfinder"]["threads"]
    log:
        f"{LOGS_DIR}/lncfinder/{SPECIES}.log",
    params:
        image=CONTAINERS["rbio_tools"],
        bind=SINGULARITY_BIND,
        script=TOOLS["lncfinder_rscript"],
    shell:
        "singularity exec {params.bind} {params.image} "
        "Rscript {params.script} -i {input.seqs} -o {output.csv} "
        "-p {threads} > {log} 2>&1"


# ---------------------------------------------------------------------------
# 8. DIAMOND homology filter: UniProt (always) + species-specific (optional)
# ---------------------------------------------------------------------------
rule diamond_uniprot:
    input:
        seqs=rules.gffread_extract_fasta.output.candidates_fasta,
    output:
        out=f"{DIRS['diamond']}/{SPECIES}_diamond_uniprot.tsv",
    threads: P["diamond"]["threads"]
    envmodules:
        MOD["diamond"],
    log:
        f"{LOGS_DIR}/diamond/{SPECIES}_uniprot.log",
    params:
        db=UNIPROT_DMND,
        sensitivity=P["diamond"]["sensitivity"],
        max_target_seqs=P["diamond"]["max_target_seqs"],
        evalue=P["diamond"]["evalue"],
    shell:
        "diamond blastx --threads {threads} -d {params.db} -q {input.seqs} "
        "-o {output.out} --outfmt 6 {params.sensitivity} "
        "--max-target-seqs {params.max_target_seqs} --evalue {params.evalue} "
        "> {log} 2>&1"


rule diamond_makedb_species:
    input:
        protein=PROTEIN_FASTA,
    output:
        db=f"{DIRS['diamond']}/{SPECIES}_protein_db.dmnd",
    threads: P["diamond"]["threads"]
    envmodules:
        MOD["diamond"],
    log:
        f"{LOGS_DIR}/diamond/{SPECIES}_makedb.log",
    params:
        db_prefix=f"{DIRS['diamond']}/{SPECIES}_protein_db",
    shell:
        "diamond makedb --in {input.protein} "
        "-d {params.db_prefix} --threads {threads} > {log} 2>&1"


rule diamond_species:
    input:
        seqs=rules.gffread_extract_fasta.output.candidates_fasta,
        db=rules.diamond_makedb_species.output.db,
    output:
        out=f"{DIRS['diamond']}/{SPECIES}_diamond_species.tsv",
    threads: P["diamond"]["threads"]
    envmodules:
        MOD["diamond"],
    log:
        f"{LOGS_DIR}/diamond/{SPECIES}_species.log",
    params:
        sensitivity=P["diamond"]["sensitivity"],
        max_target_seqs=P["diamond"]["max_target_seqs"],
        evalue=P["diamond"]["evalue"],
    shell:
        "diamond blastx --threads {threads} -d {input.db} -q {input.seqs} "
        "-o {output.out} --outfmt 6 {params.sensitivity} "
        "--max-target-seqs {params.max_target_seqs} --evalue {params.evalue} "
        "> {log} 2>&1"


rule diamond_combine:
    # Concatenating outfmt-6 tables is safe: intersect_ids.R only needs
    # unique hit IDs + evalues from column 1 / column 11.
    input:
        diamond_hits_to_combine,
    output:
        combined=f"{DIRS['diamond']}/{SPECIES}_diamond_combined.tsv",
    shell:
        "cat {input} > {output.combined}"


# ---------------------------------------------------------------------------
# 9. Consensus intersection across all classifiers
# ---------------------------------------------------------------------------
rule intersect_ids:
    input:
        feelnc=rules.feelnc_filter.output.gtf,
        cpat=rules.cpat.output.out,
        lncboost=rules.lncboost_predict.output.predictions,
        lncfinder=rules.lncfinder.output.csv,
        diamond=rules.diamond_combine.output.combined,
    output:
        final_ids=f"{DIRS['intersect']}/{SPECIES}_final_ids.txt",
        summary=f"{DIRS['intersect']}/{SPECIES}_intersection_summary.csv",
    threads: P["intersect"]["threads"]
    log:
        f"{LOGS_DIR}/intersect/{SPECIES}.log",
    params:
        image=CONTAINERS["rbio_tools"],
        script=TOOLS["intersect_ids_rscript"],
        prefix=f"{DIRS['intersect']}/{SPECIES}",
    shell:
        "singularity exec {params.image} Rscript {params.script} "
        "-f {input.feelnc} -c {input.cpat} -l {input.lncfinder} "
        "-b {input.lncboost} -d {input.diamond} -o {params.prefix} "
        "> {log} 2>&1"


rule pull_final_fasta:
    input:
        ids=rules.intersect_ids.output.final_ids,
        candidates=rules.gffread_extract_fasta.output.candidates_fasta,
    output:
        fasta=f"{DIRS['intersect']}/{SPECIES}_lncRNA.fasta",
    envmodules:
        MOD["seqkit"],
    log:
        f"{LOGS_DIR}/intersect/{SPECIES}_pull_fasta.log",
    shell:
        "seqkit grep -f {input.ids} {input.candidates} > {output.fasta} 2> {log}"


rule pull_final_gtf:
    input:
        ids=rules.intersect_ids.output.final_ids,
        candidates_gtf=rules.extract_candidates_gtf.output.candidates_gtf,
    output:
        gtf=f"{DIRS['intersect']}/{SPECIES}_lncRNA.gtf",
    shell:
        "grep -Fwf {input.ids} {input.candidates_gtf} > {output.gtf}"


# ---------------------------------------------------------------------------
# 10. Classify the final lncRNA set relative to the reference annotation
# ---------------------------------------------------------------------------
rule feelnc_classify:
    input:
        lnc_gtf=rules.pull_final_gtf.output.gtf,
        ref_gtf=REF_GTF,
    output:
        classes=f"{DIRS['classify']}/{SPECIES}_lncRNA_classes.txt",
    threads: P["feelnc_classifier"]["threads"]
    log:
        f"{LOGS_DIR}/classify/{SPECIES}_classifier.log",
    params:
        image=CONTAINERS["feelnc"],
        bind=SINGULARITY_BIND,
        feelnc_dir=TOOLS["feelnc_dir"],
        window=P["feelnc_classifier"]["window"],
        max_dist=P["feelnc_classifier"]["max_dist"],
    shell:
        "singularity exec {params.bind} "
        "--env PERL5LIB={params.feelnc_dir}/lib "
        "{params.image} perl {params.feelnc_dir}/scripts/FEELnc_classifier.pl "
        "-i {input.lnc_gtf} -a {input.ref_gtf} -w {params.window} "
        "-m {params.max_dist} -v 1 > {output.classes} 2> {log}"


rule parse_classifications:
    input:
        classes=rules.feelnc_classify.output.classes,
    output:
        done=touch(f"{DIRS['classify']}/.{SPECIES}_classification_parse.done"),
    threads: P["classify"]["threads"]
    log:
        f"{LOGS_DIR}/classify/{SPECIES}_parse.log",
    params:
        image=CONTAINERS["rbio_tools"],
        script=TOOLS["classifications_rscript"],
        outdir=DIRS["classify"],
    shell:
        "singularity exec {params.image} Rscript {params.script} "
        "-i {input.classes} -o {params.outdir} > {log} 2>&1"