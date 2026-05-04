#!/usr/bin/env python3
import os
import re
import argparse
from collections import defaultdict
import btllib
import gzip

# ----------------------
# OPEN KMER FILE
# ----------------------
def open_kmer_file(path):
    """
    Supports plain text, .gz/.bgz (gzip/bgzip) archives.
    """
    if path.endswith((".gz", ".bgz", ".gzip")):
        return gzip.open(path, "rt")

    return open(path, "r")

# ----------------------
# ARGUMENT PARSING
# ----------------------
parser = argparse.ArgumentParser(description="Build SNV window-sample matrix from Bloom filters")
parser.add_argument("--bf_dir", required=True, help="Directory containing Bloom filter files")
parser.add_argument("--kmer_file", required=True, help="File listing SNV kmers")
parser.add_argument("--output_matrix", required=True, help="Output TSV matrix file")
parser.add_argument("--window_size", type=int, default=500_000, help="Window size for aggregation")
parser.add_argument("--k", type=int, default=55, help="Bloom filter(s) k-mer length")

args = parser.parse_args()

BF_DIR = args.bf_dir
KMER_FILE = args.kmer_file
OUTPUT_MATRIX = args.output_matrix
WINDOW_SIZE = args.window_size
KMER = args.k

# ----------------------
# LIST SAMPLES / BF FILES
# ----------------------
bf_files = sorted([
    f for f in os.listdir(BF_DIR)
    if re.search(rf"_k{KMER}\.bf$", f)
])

if not bf_files:
    raise ValueError(f"No .bf files found for k={KMER} in {BF_DIR}")

samples = [re.sub(rf"_k{KMER}\.bf$", "", f) for f in bf_files]

print(f"Found {len(samples)} samples (k={KMER})", flush=True)

# ----------------------
# BUILD SNV → kmer mapping & window assignment
# ----------------------
snv_alleles = {}
snv_window = {}
windows_set = set()

with open_kmer_file(KMER_FILE) as f:
    for line in f:
        chr_, pos, allele, kmer = line.strip().split()
        key = (chr_, pos)
        if key not in snv_alleles:
            snv_alleles[key] = {}
            window_id = (chr_, int(pos)//WINDOW_SIZE)
            snv_window[key] = window_id
            windows_set.add(window_id)
        snv_alleles[key][allele] = kmer

windows_list = sorted(list(windows_set))
print(f"Total windows: {len(windows_list)}", flush=True)

# ----------------------
# Initialize dictionary for matrix
# ----------------------
matrix = {win: {s: 0.0 for s in samples} for win in windows_list}

# ----------------------
# STREAM SAMPLE BY SAMPLE
# ----------------------
for sample, bf_file in zip(samples, bf_files):
    print(f"Processing sample: {sample}", flush=True)
    bf_path = os.path.join(BF_DIR, bf_file)
    bf = btllib.KmerBloomFilter(bf_path)

    window_counts = defaultdict(list)

    for key, alleles in snv_alleles.items():
        counts = {"REF":0, "ALT":0}
        for allele, kmer in alleles.items():
            if bf.contains(kmer):
                counts[allele] += 1
        total = counts["REF"] + counts["ALT"]
        signal = (counts["ALT"] - counts["REF"]) / total if total > 0 else 0.0
        window = snv_window[key]
        window_counts[window].append(signal)

    for win, vals in window_counts.items():
        matrix[win][sample] = sum(vals)/len(vals)

# ----------------------
# WRITE OUTPUT MATRIX
# ----------------------
with open(OUTPUT_MATRIX, "w") as out:
    out.write("chr\twindow_start:" + str(WINDOW_SIZE) + "\t" + "\t".join(samples) + "\n")
    for win in windows_list:
        chr_, win_start = win
        row = [f"{matrix[win][s]:.4f}" for s in samples]
        out.write(f"{chr_}\t{win_start}\t" + "\t".join(row) + "\n")

print(f"Done! Matrix written to {OUTPUT_MATRIX}")

