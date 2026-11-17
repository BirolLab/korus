#!/usr/bin/env snakemake -s

#
# Snakefile for korusEncode
# Pipeline:
#   ntCard -> ntStat -> korusMatrix.py
#

import os
import shutil
import glob

onsuccess:
    shutil.rmtree(".snakemake", ignore_errors=True)

# ---------------------------------------------------------
# CONFIG
# ---------------------------------------------------------

reads_prefix_full = config.get("reads", "")
reads_prefix = os.path.basename(reads_prefix_full)

genome = config.get("genome", "")
genome_base = os.path.basename(os.path.realpath(genome))

k = config["kmer"]
window_size = config["window_size"]
n_threads = config["threads"]

kmer_file = config["kmer_file"]

# platform-specific time command
mac_time_command = "command time -l -o"
linux_time_command = "command time -v -o"

time_command = (
    mac_time_command
    if os.uname().sysname == "Darwin"
    else linux_time_command
)

reads_files = sorted(
    glob.glob(f"{reads_prefix}*.fastq.gz") +
    glob.glob(f"{reads_prefix}*.fq.gz")
)

# ---------------------------------------------------------
# FINAL TARGETS
# ---------------------------------------------------------

rule all:
    input:
        f"{reads_prefix}_k{k}_t{window_size}_korus-matrix.tsv"

rule genome:
    input:
        f"{genome_base}_k{k}_t{window_size}_korus-matrix.tsv"

# ---------------------------------------------------------
# ntCard
# ---------------------------------------------------------

rule ntcard:
    input:
        reads_files = reads_files 
    output:
        hist = f"{reads_prefix}_k{k}.hist",

    params:
        benchmark = f"{time_command} ntcard_k{k}.time"

    shell:
        """
        {params.benchmark} \
        ntcard -k {k} \
       	-t {n_threads} \
        -p {reads_prefix} \
        {input.reads_files} 
        """

# ---------------------------------------------------------
# ntStat / btllib bloom filter generation
# ---------------------------------------------------------

rule ntstat_bloom:
    input:
        hist = f"{reads_prefix}_k{k}.hist",
        reads_files = reads_files

    output:
        bf = temp(f"{reads_prefix}_k{k}.bf")

    params:
        benchmark = f"{time_command} ntstat_k{k}.time",
        min_cutoff = f"-cmin 0"

    shell:
        """
        {params.benchmark} \
        ntstat filter \
        -k {k} \
        {params.min_cutoff} \
        -t {n_threads} \
        -f {input.hist} \
        -o {output.bf} \
        {input.reads_files} 
        """


# ---------------------------------------------------------
# Make Bloom filter from genome assembly
# ---------------------------------------------------------
rule genome_ntcard:
    input: 
        genome = genome
    output:
        hist=f"{genome_base}.k{k}.hist",
        f_stats=f"{genome_base}.k{k}.tsv"
    params:
        benchmark = f"{time_command} ntcard_{genome_base}_k{k}.time",
        options = f"-t {n_threads} -k {k}"
    shell:
        "{params.benchmark} ntcard {params.options} -o {output.hist} {input.genome} &> {output.f_stats}"

rule ntedit_genome_bf:
    input:
        ntcard_f_stats = f"{genome_base}.k{k}.tsv",
        genome = genome
    output:
        f"{genome_base}_k{k}.bf"
    params:
        options = f"-k {k} -t {n_threads}",
        benchmark = f"{time_command} make_bf_{genome_base}_k{k}.time"
    run:
        with open(input.ntcard_f_stats, 'r', encoding="utf8") as fin:
            for line in fin:
                _, f, num = line.strip().split("\t")
                if f == "F0":
                    num_elements = num
                    break
        shell("{params.benchmark} ntedit_make_genome_bf --genome {input.genome} {params.options} --num_elements {num_elements} -o {output}")


# ---------------------------------------------------------
# Build KORUS matrix
# ---------------------------------------------------------

rule korus_matrix:
    input:
        bf = f"{reads_prefix}_k{k}.bf",
        kmer_file = kmer_file

    output:
        matrix = f"{reads_prefix}_k{k}_t{window_size}_korus-matrix.tsv"

    params:
        benchmark = f"{time_command} korusMatrix_k{k}.time"

    shell:
        """
        {params.benchmark} \
        korusMatrix.py \
            --bf_dir . \
            --sample {reads_prefix:q} \
            --kmer_file {input.kmer_file} \
            --output_matrix {output.matrix} \
            --window_size {window_size} \
            --k {k}
        """

rule korus_matrix_genome:
    input:
        bf = f"{genome_base}_k{k}.bf",
        kmer_file = kmer_file

    output:
        matrix = f"{genome_base}_k{k}_t{window_size}_korus-matrix.tsv"

    params:
        benchmark = f"{time_command} korusMatrix_k{k}.time"

    shell:
        """
        {params.benchmark} \
        korusMatrix.py \
            --bf_dir . \
            --sample {genome_base:q} \
            --kmer_file {input.kmer_file} \
            --output_matrix {output.matrix} \
            --window_size {window_size} \
            --k {k}
        """