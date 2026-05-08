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

reads_prefix_full = config["reads"]
reads_prefix = os.path.basename(reads_prefix_full)

k = config["kmer"]
tile_size = config["tile_size"]
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

reads_files = sorted(glob.glob(f"{reads_prefix}_*.fastq.gz"))

# ---------------------------------------------------------
# FINAL TARGET
# ---------------------------------------------------------

rule all:
    input:
        f"{reads_prefix}_k{k}_t{tile_size}_korus-matrix.tsv"

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
        ntcard -k {k} -t {n_threads} \
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
        -k {k} {params.min_cutoff} \
        -t {n_threads} \
	-f {input.hist} \
        -o {output.bf} \
        {input.reads_files} 
        """

# ---------------------------------------------------------
# Build KORUS matrix
# ---------------------------------------------------------

rule korus_matrix:
    input:
        bf = f"{reads_prefix}_k{k}.bf",
        kmer_file = kmer_file

    output:
        matrix = f"{reads_prefix}_k{k}_t{tile_size}_korus-matrix.tsv"

    params:
        benchmark = f"{time_command} korusMatrix_k{k}.time"

    shell:
        """
        {params.benchmark} \
        korusMatrix.py \
            --bf_dir . \
            --kmer_file {input.kmer_file} \
            --output_matrix {output.matrix} \
            --window_size {tile_size} \
            --k {k}
        """
