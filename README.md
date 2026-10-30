[![Release](https://img.shields.io/github/release/birollab/korus.svg)](https://github.com/birollab/korus/releases)
[![ntRoot dependency](https://img.shields.io/badge/depends-ntRoot%20framework-brightgreen)](https://github.com/birollab/ntRoot)
[![Zenodo](https://img.shields.io/badge/korus-zenodo-red)](https://doi.org/10.5281/zenodo.23738447)
[![Python](https://img.shields.io/badge/python-3.9%2B-blue)]()
[![License](https://img.shields.io/badge/license-GPLv3-lightgrey)]()

![Logo](https://github.com/warrenlr/korus/blob/main/korus-logo.png)

# KORUS

Probabilistic ancestry inference beyond variant calls

KORUS (K-mer Orthogonal Representation for Unified Stratification of ancestry) is a probabilistic framework for genomic ancestry inference using k-mer signals encoded in Bloom filters. It supports local and genome-wide ancestry inference without requiring explicit variant calling or sequence alignment of query genomes.

Feature Matrix → Probabilistic Decoding Framework for Genomic Ancestry Inference

---

## Contents

1. [Overview](#overview)
2. [Key Methodological Features](#features)
3. [Credit](#credit)  
4. [Methodological Framework](#methodological-framework)  
5. [KORUS Modes](#korus-modes)
6. [Installation](#installation)  
7. [Tested Environment](#tested-environment)
8. [Requirements](#req)  
9. [Input Files](#inputs)  
10. [Running KORUS using one query sample with the provided reference matrix](#running-korus-using-one-query-sample-with-the-provided-reference-matrix)
11. [Advanced usage: Generating a custom reference matrix with KORUS](#advanced-usage-generating-a-custom-reference-matrix-with-korus)
12. [Advanced usage: Running KORUS using multiple samples](#advanced-usage-running-korus-using-multiple-samples)
13. [Interactive Visualizations](#viz)
14. [Citing](#citing)  
15. [License](#license)  

---

## Overview <a name="overview"></a><a name="synopsis"></a>

KORUS is a two-stage computational framework for ancestry inference from sequencing data using k-mer–derived genomic features encoded in Bloom filters.

Rather than relying on alignment or variant calling, KORUS transforms sequencing data into a structured genomic feature space and performs probabilistic decoding to infer ancestry structure across the genome.

The framework consists of:

- **MATRIX stage** — construction of k-mer–derived genomic feature matrices. KORUS converts k-mer presence in Bloom filters into normalized, per-window ALT–REF scalar signals and aggregates these into a sample-by-window feature matrix.
- **DECODE stage** — probabilistic inference of ancestry from feature space. The feature space is treated as an orthogonalized ancestry embedding in PCA-reduced coordinates, enabling covariance-aware distance classification. In other words, the feature space is treated as a decorrelated representation in PCA space prior to Mahalanobis distance estimation. Ancestry scores are computed using Mahalanobis distances in this PCA-reduced feature space and transformed via softmax normalization to obtain probabilistic ancestry estimates.

An optional hidden Markov model (HMM), decoded using the Viterbi algorithm, smooths window-wise ancestry paths before genome-wide ancestry fractions are calculated.

### Key Methodological Features <a name="features"></a>

- k-mer/Bloom-filter ancestry representation
- Reference-free signal extraction
- Alignment-free local ancestry inference
- Soft-supervised ancestry manifold modeling

---

## Credit <a name="credit"></a>

Developed by René L Warren

---

## Methodological Framework <a name="methodological-framework"></a>

### MATRIX stage (korusMatrix.py)

Transforms Bloom filter–encoded sequencing data into a structured genomic feature representation.

**Key steps:**
- Build per-sample k-mer Bloom filters (ntEdit/ntStat compatible)
- Query per-sample Bloom filters for ancestry-informative SNV k-mer membership 
- Compute allele-specific support (REF vs ALT k-mers)
- Derive normalized directional signal per locus
- Aggregate signals into fixed genomic windows
- Construct samples × window feature matrix

**Output:**  
* Continuous genomic feature matrix across samples and windows


### DECODE stage (korusDecode.py)

Performs probabilistic inference of ancestry from the feature matrix.

**Key Steps:**
- Standardization of the feature matrix
- Dimensionality reduction (PCA)
- Distance-based classification (Mahalanobis)
- Softmax normalization of ancestry likelihoods (Distance-to-probability transformation)
- Optional unsupervised reference clustering (K-means) 
- Spatial smoothing (HMM / Viterbi)

**Output:**
- Global ancestry fractions  
- Local ancestry assignments  
- Distance-derived ancestry probability matrices  
- Smoothed ancestry trajectories  

---

### KORUS Modes

#### Hard supervised (`--mode hard` in `korus`)
- Uses discrete reference labels (`--labels`)
- Performs supervised ancestry classification
- KORUS learns ancestry cluster structure from categorical population labels and performs window-level classification by minimizing Mahalanobis distance to class-specific distributions in PCA space
- This framework assumes each reference genome belongs predominantly to one ancestry class

#### Soft supervised (`--mode soft` in `korus`)
- Uses fractional ancestry supervision (`--fractions`)
- Supports admixed reference modeling
- Instead of learning discrete ancestry clusters from hard assignments, KORUS estimates ancestry distributions using weighted contributions from admixed reference genomes

#### Unsupervised (`--mode unsupervised` in `korus`)
- Partitions the reference panel into `--c` K-means clusters and uses the resulting cluster assignments as pseudo-labels for ancestry inference 
- Operates without explicit ancestry labels
- The inferred clusters represent data-driven ancestry groups rather than predefined ancestry categories

---

## Installation <a name="installation"></a>

```bash
git clone https://github.com/birollab/korus.git
cd korus
# Ensure the KORUS repository is on your PATH
export PATH=$(PWD):$PATH

conda install -c conda-forge -c bioconda numpy pandas scikit-learn matplotlib plotly btllib ntroot snakemake ntcard ntstat
```

---

## Tested Environment

KORUS was developed and tested on:
- Rocky Linux 9.6 (Blue Onyx)
- Python 3.12.3

---

## Requirements <a name="req"></a>

### Core software
- Python >= 3.9
- btllib >= 1.7
- snakemake
- ntCard >= 1.0.1
- ntStat

### Python packages
- numpy
- pandas
- scikit-learn
- matplotlib
- plotly
- packaging

### Optional (needed for `korus`)
- ntRoot >= 1.2.2

---
## Input Files <a name="inputs"></a>

### Reads

Input reference/query reads in FASTQ format

### SNV k-mer File
See `tools/README` for instructions for generating this file from a VCF file.
```
chr  pos  allele  kmer
```
e.g.,
```
10	11501	REF	CCAGCGGTGTCCGGGAGCCTAGCGGCGCGTGTGCAGCGGCCAGTGCACCTGCTCT
10	11501	ALT	CCAGCGGTGTCCGGGAGCCTAGCGGCGAGTGTGCAGCGGCCAGTGCACCTGCTCT
```

### Reference Hard Labels (Hard supervised mode only, `--labels`)
```
sample   ancestry
ERR1     EUR
ERR2     AFR
```

### Reference Soft Labels (Soft supervised mode only, `--fractions`)
```
sample	 EAS	  SAS	   AMR	   AFR	   EUR
ERR1	 0.1961	0.0105	0.7233	0.0595	0.0105
ERR2	 0.0718	0.0298	0.5814	0.2347	0.0823
```
---

## Running KORUS using one query sample with the provided reference matrix

### End-to-End KORUS (`korus`) Pipeline <a name="pipeline"></a>

KORUS provides a single-entry driver script to run all steps in the pipeline for a single query sample:

```
usage: korus [-h] --reads READS [--report_mode {iframe,embed}] --kmer_file KMER_FILE -k {45,55}
             [--window WINDOW] [-t T] [--mode {hard,soft,unsupervised}] [--c C] [--global_pcs GLOBAL_PCS]
             [--local_pcs LOCAL_PCS] [--window_radius WINDOW_RADIUS] [--stay_prob STAY_PROB] [--no_hmm]
             --outprefix OUTPREFIX [--report_title REPORT_TITLE] [-n] [-V] [--ref REF] [--labels LABELS]
             [--fractions FRACTIONS]
End-to-end KORUS ancestry inference pipeline

options:
  -h, --help            show this help message and exit
  --reads READS         Input read prefix
  --report_mode {iframe,embed}
                        Report layout mode [default=iframe]
  --kmer_file KMER_FILE
                        SNV k-mer definition file
  -k {45,55}            Model k-mer size
  --window WINDOW       Window size [default 500000bp]
  -t T                  Number of threads [default 48]
  --mode {hard,soft,unsupervised}
                        KORUS mode (hard, soft or unsupervised)
  --c C                 KMeans n_clusters (unsupervised mode)
  --global_pcs GLOBAL_PCS
                        Number of PCs used for unsupervised K-means clustering
  --local_pcs LOCAL_PCS
                        Maximum number of PCs used for local Mahalanobis ancestry inference
  --window_radius WINDOW_RADIUS
                        Window radius [default 3]
  --stay_prob STAY_PROB
                        HMM self-transition probability [default 0.98]
  --no_hmm              Disable HMM smoothing
  --outprefix OUTPREFIX
                        Prefix for output files
  --report_title REPORT_TITLE
                        Title for output report [KORUS Ancestry Predictions]
  -n, --dry-run         Print commands without executing them
  -V, --version         show program's version number and exit
  --ref REF             Override bundled reference matrix
  --labels LABELS       Override bundled hard-supervised labels
  --fractions FRACTIONS
                        Override bundled soft-supervised fractions
```
Example command:
```bash
korus --reads SAMPLE_PREFIX --kmer_file snv_kmers.tsv -k 55 --outprefix run1 --mode hard
```

---
## Advanced usage: Generating a custom reference matrix with KORUS

1. Run `korusEncode` for all input reference samples - runs `ntCard`, `ntStat` and `korusMatrix.py`
```bash
usage: korusEncode [-h] --reads READS --kmer_file KMER_FILE -k K [--window_size WINDOW_SIZE] [-t T] [-n] [-f] [-V]

korusEncode: build KORUS feature matrices from sequencing reads

options:
  -h, --help            show this help message and exit
  --reads READS         Prefix of input read files (FASTQ/FASTA/gz/bz2 compatible)
  --kmer_file KMER_FILE
                        SNV k-mer definition file used by korusMatrix.py
  -k K                  k-mer size
  --window_size WINDOW_SIZE
                        Window size for matrix aggregation [default=500000]
  -t T                  Number of threads [default=4]
  -n, --dry-run         Print commands without executing
  -f, --force           Force rerun of all rules
  -V, --version         show program's version number and exit

Input reads are identified using a shared prefix.
```

Example command:
```
korusEncode \
  --reads SAMPLE_PREFIX (Prefix of input read files (FASTQ/FASTA/gz/bz2 compatible)) \
  --kmer_file snv_kmers.tsv -k 55 -t 48
```

2. Combine all reference matrices to one single reference matrix
```
usage: korusAggregateMatrices.py [-h] [--matrix-list FILE] [-o FILE]
                            [matrices [matrices ...]]

Merge single-sample matrices with matching window sizes. Preserve input sample order and first-seen coordinate order; fill missing values with NA.

positional arguments:
  matrices              Input matrix paths

optional arguments:
  -h, --help            show this help message and exit
  --matrix-list FILE    File with input matrices (one per line)
  -o FILE, --output FILE
                        Output TSV path (default: standard output; '-' also
                        means stdout)
                        
Examples:
python3 korusAggregateMatrices.py sample1.tsv sample2.tsv -o merged.tsv
python3 korusAggregateMatrices.py --matrix-list matrices.txt -o merged.tsv

Each matrix must have a header: chr window_start:<size> <sample>.
```

This matrix can be input to `korus` using `--ref`. For usage with multi-sample KORUS runs, see below.

---

## Advanced usage: Running KORUS using multiple samples

While the `korus` script can be used to run all steps using a single sample, these steps can be followed to generate a single sample matrix for analysis using the KORUS pipeline.

1. Follow the same steps as for generating a custom reference matrix (see above for details):

    a. Run `korusEncode` for each input sample data, using the same parameters
  
    b. Combine all sample matrices to a combined multi-sample matrix

2. Run `korusDecode.py` to generate ancestry fraction predictions for each input sample
```bash
usage: korusDecode.py [-h] --ref REF --sample SAMPLE [--labels LABELS] [--fractions FRACTIONS] [--c C]
                      --outprefix OUTPREFIX [--window_radius WINDOW_RADIUS] [--require_full_window]
                      [--stay_prob STAY_PROB] [--no_hmm] [--global_pcs GLOBAL_PCS]
                      [--local_pcs LOCAL_PCS] [-V]
```

For hard supervision, use `--labels`, and for soft supervision use `--fractions`. Omitting both of these parameters will run KORUS in unsupervised mode.

Use `--ref` to input the reference matrix (available in directory `model` or use custom-built matrix).

Example command:
```
korusDecode.py --ref model/korusModelMatrixk55.tsv.gz --sample sampleMatrix.tsv \
  --fractions korusModelLabelsFractions.tsv.gz --global_pcs 10 --local_pcs 10 \
  --outprefix korus_run --stay_prob 0.98
```

PCA parameters:
 `--global_pcs` specifies the number of principal components used during unsupervised K-means clustering of the reference panel. `--local_pcs` specifies the maximum number of principal components retained for the window-level Mahalanobis ancestry inference performed during local ancestry assignment. These parameters control distinct PCA operations and may be adjusted independently.


---

### Output Files

- `.lai_summary.tsv` — ancestry fractions and probabilities for each query sample
- `*.lai.tsv` — local ancestry inference results per sample
- `*.pca3d.html` — population context (PCA embedding)
- `.reference_clusters.tsv` — reference sample IDs and their inferred K-means cluster assignments (unsupervised mode only)

- `korus` pipeline (single query sample) only:
  - `*.lai.html` — interactive genome view
  - `*.report.html` — combined PCA + LAI interactive visualization report (see below examples)

---


## Interactive Visualizations <a name="viz"></a>

As part of `korusDecode.py`, KORUS generates interactive PCA visualizations that position the query sample in relation to reference populations in ancestry space. The visualization examples below show sample ERR3988858 from Puerto Rico.

### PCA Visualizations

These HTML files (`*pca3d.html`) allow exploration of population structure under different inference modes:

- **Hard mode**
  [Interactive PCA (hard)](https://www.birollab.ca/viz/korus7Ehard.pca3d.html)

- **Soft mode**
  [Interactive PCA (soft / fractions)](https://www.birollab.ca/viz/korus7Esoft.fractions.pca3d.html)

- **Unsupervised mode**
  [Interactive PCA (unsupervised K=5)](https://www.birollab.ca/viz/korus7Eunsup.pca3d.html)
  [Interactive PCA (unsupervised K=7)](https://www.birollab.ca/viz/korus7EunsupC7.pca3d.html)
---

### Local Ancestry (LAI) Visualization

The end-to-end `korus` pipeline also produces interactive chromosomal local ancestry inference (LAI) visualizations.

These allow users to explore ancestry assignments along chromosomes or across the genome, with linked interaction to population structure views.

Reference implementation and details:
https://github.com/BirolLab/ntRoot/blob/main/README.md#interactive-local-ancestry-lai-visualization-

- **LAI visualization (hard mode example)**
  [Interactive LAI (hard)](https://www.birollab.ca/viz/ERR3988858_korus7Ehard.lai.html)

---

### Combined HTML Report

KORUS generates a standalone report that integrates PCA and LAI views into a unified interface.

- **KORUS report (hard mode example)**
  [Interactive report](https://www.birollab.ca/viz/korus7Ehard.report.html)

---

### Notes

- PCA views are produced directly by `korusDecode.py`
- LAI and report visualizations are generated when running the full `korus` pipeline
- The report provides synchronized navigation between genome view and population context

---

## Citing <a name="citing"></a>

TBD 

---

## License <a name="license"></a>

KORUS is released under the GNU General Public License v3 (GPLv3).

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, version 3.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

A copy of the GPLv3 license (LICENSE) can be found here: https://github.com/BirolLab/korus/LICENSE

See http://www.gnu.org/licenses/ for details.
