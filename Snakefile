"""
Author: Sam Nooij
Organisation: Leiden University Medical Center (LUMC)
Department: Medical Microbiology (MM) Research,
    Center for Microbiome Analyses and Therapeutics (CMAT)
Date: 2021-05-18

Workflow for general preprocessing of metagenomic reads,
followed by Mash sketching and Jaccard distance calculation.
Steps included:
 - Host genome removal (i.e. human)
 - Quality trimming, low-complexity filtering, and adapter removal
 - Statistics reporting (number of reads removed/trimmed, etc.)
 - Mash sketch
 - Mash pairwise distance calculation

Input: Raw metagenomics reads (as fastq or fastq.gz)
Output: Several gzipped fastq files and HTML reports with statistics,
 and a tab-separated table with Jaccard distances, p-values and number of shared sketches.

Example use:
    $ snakemake --profile config

N.B. Environment-specific variables may be set in the configuration files
    in the 'config' folder.
"""


### Step 1: Import configuration file ###

configfile: Path("config/parameters.yaml")

SAMPLES = config["samples"]

#Scan for reference files
REFERENCE_DIR = config["host_database"]

#Make sure a slash is appended to the end of the directory name
if not REFERENCE_DIR.endswith("/"):
    REFERENCE_DIR = REFERENCE_DIR + "/"
else:
    pass


### Step 2: Specify the final output files ###

rule all:
   input:
        expand(config["filtered_directory"] + "{sample}-filtered.fastq.gz",
        sample = SAMPLES),
        # Host-filtered reads

        expand(config["trimmed_directory"] + "{sample}-trimmed.fastq.gz",
        sample = SAMPLES),
        # Quality-trimmed reads

        "doc/Preprocessing_report.html",
        # Report with preprocessing statistics

        expand("data/mash/{sample}.msh",
        sample = SAMPLES),
        # Mash sketches

        "data/processed/Mash_distances.tsv",
        # Mash distance table


### Step 3: include processing steps ###

rule host_removal:
    input:
        reads=config["raw_directory"] + "{sample}.fastq.gz",
        reference=REFERENCE_DIR + "reference.1.bt2"
    output:
        filtered=config["filtered_directory"] + "{sample}-filtered.fastq.gz"
    conda:
        "envs/host_removal.yaml"
    threads: config["host_removal"]["threads"]
    resources:
        mem=int(config["host_removal"]["memory"]),
        time=int(config["host_removal"]["time"])
    log:
        "log/preprocessing/host_removal-{sample}.txt"
    benchmark:
        "log/preprocessing/benchmark/host_removal-{sample}.txt"
    params:
        reference=REFERENCE_DIR + "reference",
        mode=config["bowtie2"]["mode"],
        compression_level=config["compression_level"]
    shell:
        """
bowtie2 --threads {threads} {params.mode} -x {params.reference}\
  -U {input.reads} 2> {log} |\
 samtools sort -@ {threads} -n - 2>> {log} |\
 samtools fastq -@ {threads} -c {params.compression_level} -f 4 \
  > {output.filtered} 2>> {log}
        """


rule quality_filtering:
    input:
        filtered=config["filtered_directory"] + "{sample}-filtered.fastq.gz",
    output:
        trimmed=config["trimmed_directory"] + "{sample}-trimmed.fastq.gz",
        json="data/fastp/{sample}-fastp.json",
        html="data/fastp/{sample}-fastp.html"
    conda:
        "envs/quality_filtering.yaml"
    threads: config["quality_filtering"]["threads"]
    resources:
        mem=int(config["quality_filtering"]["memory"]),
        time=int(config["quality_filtering"]["time"])
    log:
        "log/preprocessing/quality_filtering-{sample}.txt"
    benchmark:
        "log/preprocessing/benchmark/quality_filtering-{sample}.txt"
    params:
        trim=config["fastp"]["trim"],
        minlen=config["fastp"]["minlen"],
        adapter=config["fastp"]["adapter"],
        complexity=config["fastp"]["complexity"]
    shell:
        """
fastp -w {threads} --in1 {input.filtered}\
 --out1 {output.trimmed}\
 {params.trim} {params.minlen} {params.adapter} {params.complexity}\
 -j {output.json} -h {output.html}
        """


rule generate_report:
    input:
        filter=expand("log/preprocessing/host_removal-{sample}.txt",
          sample = SAMPLES),
        trim=expand("data/fastp/{sample}-fastp.json",
          sample = SAMPLES)
    output:
        "doc/Preprocessing_report.html"
    conda:
        "envs/report.yaml"
    threads: config["generate_report"]["threads"]
    resources:
        mem=int(config["generate_report"]["memory"]),
        time=int(config["generate_report"]["time"])
    log:
        "log/preprocessing/generate_report.txt"
    benchmark:
        "log/preprocessing/benchmark/generate_report.txt"
    params:
        config=config["multiqc"]["config"]
    shell:
        """
multiqc log/preprocessing/ data/fastp/ --config {params.config}\
 -o doc/ -n Preprocessing_report.html
        """


rule create_sketch:
    input:
        config["trimmed_directory"] + "{sample}-trimmed.fastq.gz",
    output:
        "data/mash/{sample}.msh"
    conda:
        "envs/mash.yaml"
    threads: config["create_sketch"]["threads"]
    resources:
        mem=int(config["create_sketch"]["memory"]),
        time=int(config["create_sketch"]["time"])
    log:
        "log/create_sketch-{sample}.txt"
    benchmark:
        "log/benchmark/create_sketch-{sample}.txt"
    params:
        sketch=config["sketch"]
    shell:
        """
mash sketch -I {wildcards.sample} {params.sketch} -o data/mash/{wildcards.sample} {input}
        """


rule calculate_distances:
    input:
        expand("data/mash/{sample}.msh",
        sample = SAMPLES)
    output:
        "data/processed/Mash_distances.tsv"
    conda:
        "envs/mash.yaml"
    threads: config["calculate_distances"]["threads"]
    resources:
        mem=int(config["calculate_distances"]["memory"]),
        time=int(config["calculate_distances"]["time"])
    log:
        "log/calculate_distances.txt"
    benchmark:
        "log/benchmark/calculate_distances.txt"
    shell:
        """
mash triangle -E {input} > {output}
        """
