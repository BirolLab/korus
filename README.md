[![ntRoot dependency](https://img.shields.io/badge/depends-ntRoot%20framework-brightgreen)](https://github.com/birollab/ntRoot)
[![Python](https://img.shields.io/badge/python-3.9%2B-blue)]()
[![License](https://img.shields.io/badge/license-GPLv3-lightgrey)]()

![Logo](https://github.com/warrenlr/korus/blob/main/korus-logo.png)

# KORUS

K-mer Orthogonal Representation for Unified Stratification (of ancestry) — an orthogonalized ancestry feature representation in PCA-reduced space

KORUS derives window-level ancestry by embedding k-mer-based ALT/REF signals into a feature matrix, assigning ancestry likelihoods via Mahalanobis distance in PCA-reduced feature space followed by softmax normalization, and enforcing genomic continuity using HMM smoothing before computing ancestry fractions

## Feature Matrix → Probabilistic Decode Framework for Genomic Ancestry Inference

---

## Contents

1. [Synopsis](#synopsis)
2. [Features](#features)
3. [Credit](#credit)  
4. [Overview](#overview)  
5. [Methodological Framework](#methodological-framework)  
6. [Installation](#installation)  
7. [Requirements](#req)  
8. [Usage](#usage)  
9. [Inputs](#inputs)  
10. [Outputs](#outputs)  
11. [Workflow Overview](#workflow-overview)  
12. [End-to-End Pipeline (korus)](#pipeline)
13. [Interactive Visualizations](#viz)
14. [Demo](#demo)  
15. [Citing](#citing)  
16. [License](#license)  

---

## Synopsis <a name="synopsis"></a>

KORUS converts k-mer presence in Bloom filters into per-window scalar signals (ALT–REF normalized), aggregates these into a sample-by-window feature matrix, assigning ancestry likelihoods via Mahalanobis distance in PCA-reduced feature space followed by softmax normalization, and optionally applies an HMM (Viterbi) to smooth window-wise ancestry paths before aggregating them into genome-wide fractions


### Key Methodological Features <a name="features"></a>

- k-mer/Bloom-filter ancestry representation
- Reference-free signal extraction
- Local ancestry inference without alignment
- Soft-supervised ancestry manifold modeling

---

## Credit <a name="credit"></a>

Developed by René L Warren


---

## Overview <a name="overview"></a>

KORUS is a two-stage computational framework for ancestry inference from sequencing data using k-mer–derived genomic features encoded in Bloom filters

Rather than relying on alignment or variant calling, KORUS transforms sequencing data into a structured genomic feature space and performs probabilistic decoding to infer ancestry structure across the genome

The feature space is treated as an orthogonalized ancestry embedding in PCA-reduced coordinates, enabling covariance-aware distance classification. In other words, the feature space is treated as a decorrelated representation in PCA space prior to Mahalanobis distance estimation

The framework consists of:

- **MATRIX stage** — construction of k-mer–derived genomic feature matrices  
- **DECODE stage** — probabilistic inference of ancestry from feature space  

---

## Methodological Framework <a name="methodological-framework"></a>

### MATRIX (korusMatrix.py)

Transforms Bloom filter–encoded sequencing data into a structured genomic feature representation.

**Key steps:**
- Build per-sample k-mer Bloom filters (ntEdit/ntStat compatible)
- Query per-sample Bloom filters for ancestry-informative SNV k-mer membership 
- Compute allele-specific support (REF vs ALT k-mers)
- Derive normalized directional signal per locus
- Aggregate signals into fixed genomic windows
- Construct samples × window feature matrix

**Output:**  
Continuous genomic feature matrix across samples and windows

---

### DECODE (korusDecode.py)

Performs probabilistic inference of ancestry from the feature matrix.

#### Hard Supervision (--labels)

Each reference sample is assigned a single discrete ancestry class

KORUS learns ancestry cluster structure from categorical population labels and performs window-level classification by minimizing Mahalanobis distance to class-specific distributions in PCA space

This framework assumes each reference genome belongs predominantly to one ancestry class

In practice this behaves similarly to:
- Supervised metric learning
- Centroid/covariance-based discriminant assignment
- Gaussian distance classification in reduced-dimensional space

with optional HMM smoothing across genomic windows.


#### Soft Supervision (--fractions)

Each reference sample contributes fractional ancestry memberships rather than a single categorical label

Instead of learning discrete ancestry clusters from hard assignments, KORUS estimates ancestry distributions using weighted contributions from admixed reference genomes

Mahalanobis distances are therefore computed relative to ancestry-specific distributions informed by continuous ancestry weights

This allows:
- Admixed references
- Continuous ancestry structure
- Smoother ancestry manifolds
- Reduced dependence on discrete continental boundaries

Conceptually, this behaves like:
- Weighted Gaussian supervision
- Probabilistic prototype learning
- Soft-label metric embedding
Rather than categorical classification

#### Unsupervised (--mode unsupervised)

No ancestry labels or ancestry fractions are required

KORUS first standardizes the reference feature matrix, performs PCA, and partitions the reference samples into `k` clusters using K-means (see below, --c in korusDecode/korus). The resulting cluster assignments are then used as pseudo-labels for the same Gaussian/Mahalanobis classification framework used in hard-supervised mode

This allows ancestry structure to be inferred directly from the reference panel without predefined ancestry labels

Conceptually, this behaves like:
- Unsupervised population structure discovery
- K-means clustering in PCA space
- Pseudo-label Gaussian classification

The inferred clusters represent data-driven ancestry groups rather than predefined ancestry categories


#### Hard vs Soft Supervision Behavior

When using --labels, KORUS can model any discrete ancestry category represented in the reference labels

When using --fractions, the ancestry space is constrained to the ancestry columns provided in the fraction matrix. Ancestries absent from the fraction table are not explicitly modeled during probabilistic decoding


**Core Steps:**
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

## Conceptual Model

Sequencing Data  
→ k-mer Bloom Filter Space  
→ **MATRIX**  
→ **DECODE**

---

## Installation <a name="installation"></a>

```bash
git clone https://github.com/birollab/korus.git
cd korus

pip install numpy pandas scikit-learn matplotlib plotly btllib 

chmod +x korusMatrix.py korusDecode.py
```
Alternatively, installation can be performed within a Conda environment.
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

### Python packages
- numpy
- pandas
- scikit-learn
- matplotlib
- plotly

### Input Data Requirements
KORUS expects Bloom filters generated using a compatible ntHash/ntCard/ntEdit-style workflow.

Compatible `.bf` files can be generated using:
- ntEdit
- ntStat
- btllib utilities

### Python standard library modules used
The following standard-library modules are used and do not require installation:
- argparse
- collections
- gzip
- io
- math
- os
- re
- sys

---

## Usage <a name="usage"></a>

### 1. MATRIX Stage

The SNV-to-signal encode runs via the `korusEncode` snakemake (runs `ntCard`, `ntStat` and `korusMatrix.py`):
```bash
usage: korusEncode [-h] --reads READS --kmer_file KMER_FILE -k K [--tile TILE] [-t T] [-n] [-f] [-V]
korusEncode: the following arguments are required: --reads, --kmer_file, -k 
```

```
./korusEncode
  --reads Prefix of input read files (FASTQ/FASTA/gz/bz2 compatible)
  --kmer_file snv_kmers.tsv
  -k 55
  -t 48
```

```
usage: korusMatrix.py [-h] --bf_dir BF_DIR --kmer_file KMER_FILE --output_matrix OUTPUT_MATRIX
                      [--window_size WINDOW_SIZE] [--k K] [-V]
korusMatrix.py: the following arguments are required: --bf_dir, --kmer_file, --output_matrix
```

```
./korusMatrix.py 
  --bf_dir /path/to/bloom_filters 
  --snv_kmers snv_kmers.tsv 
  --output_matrix sampleMatrix.tsv 
  --window_size 500000 
  --k 55
```

### 2. DECODE Stage

```bash
usage: korusDecode.py [-h] --ref REF --sample SAMPLE [--labels LABELS] [--fractions FRACTIONS] [--c C]
                      --outprefix OUTPREFIX [--window_radius WINDOW_RADIUS] [--require_full_window]
                      [--stay_prob STAY_PROB] [--no_hmm] [--global_pcs GLOBAL_PCS]
                      [--local_pcs LOCAL_PCS] [-V]
korusDecode.py: the following arguments are required: --ref, --sample, --outprefix
```

```
./korusDecode.py
  --ref korusModelMatrixk55.tsv
  --sample sampleMatrix.tsv
  --labels korusModelLabels.tsv
  --fractions korusModelLabelsFractions.tsv.gz
  --c 5
  --global_pcs 10
  --local_pcs 10
  --outprefix korus_run
  --stay_prob 0.98
```

```
PCA parameters. --global_pcs specifies the number of principal components used during unsupervised K-means clustering of the reference panel. --local_pcs specifies the maximum number of principal components retained for the window-level Mahalanobis ancestry inference performed during local ancestry assignment. These parameters control distinct PCA operations and may be adjusted independently.
```

---

## Input Files <a name="inputs"></a>

### Bloom filters
Per-sample k-mer Bloom filters compatible with ntEdit/ntStat.

Example:
```
SAMPLE_k55.bf
```

### SNV k-mer File
```
chr  pos  allele  kmer
```
e.g.,
```
10	11501	REF	CCAGCGGTGTCCGGGAGCCTAGCGGCGCGTGTGCAGCGGCCAGTGCACCTGCTCT
10	11501	ALT	CCAGCGGTGTCCGGGAGCCTAGCGGCGAGTGTGCAGCGGCCAGTGCACCTGCTCT
```

### Reference Hard Labels (Optional, --labels)
```
sample   ancestry
ERR1     EUR
ERR2     AFR
```

### Reference Soft Labels (Optional, --fractions)
```
sample	 EAS	  SAS	   AMR	   AFR	   EUR
ERR1	 0.1961	0.0105	0.7233	0.0595	0.0105
ERR2	 0.0718	0.0298	0.5814	0.2347	0.0823
```
---

## Output Files <a name="outputs"></a>

### Global Outputs
- .lai_summary.tsv
- .ancestry_fractions.tsv  
- .posterior_probabilities.tsv  

### Local Ancestry Outputs
- .SAMPLE.lai.tsv   

### Unsupervised Mode (additional output)
- .reference_clusters.tsv — reference sample IDs and their inferred K-means cluster assignments

### Visualizations
- PCA plots  
- ancestry trajectories  

### Logs
- runtime logs per stage  

---

## Workflow Overview <a name="workflow-overview"></a>

1. Bloom filter interrogation
2. SNV-level signal extraction
3. Window aggregation (MATRIX)
4. PCA embedding
5. Probabilistic classification
6. Optional smoothing (DECODE)
7. Report generation (korus pipeline, below)

---

## End-to-End KORUS (korus) Pipeline <a name="pipeline"></a>

In addition to the modular workflow (`korusMatrix.py` → `korusDecode.py`), KORUS provides a single-entry driver script:

```bash
./korus
  --reads SAMPLE_PREFIX
  --kmer_file snv_kmers.tsv
  -k 55
  --outprefix run1
  --mode hard
```

```
usage: korus [-h] --reads READS [--report_mode {iframe,embed}] --kmer_file KMER_FILE -k {45,55}
             [--tile TILE] [-t T] [--mode {hard,soft,unsupervised}] [--c C] [--global_pcs GLOBAL_PCS]
             [--local_pcs LOCAL_PCS] [--window_radius WINDOW_RADIUS] [--stay_prob STAY_PROB] [--no_hmm]
             --outprefix OUTPREFIX [--report_title REPORT_TITLE] [-n] [-V] [--ref REF] [--labels LABELS]
             [--fractions FRACTIONS]
korus: the following arguments are required: --reads, --kmer_file, -k, --outprefix
```

### What It Does

The `korus` driver orchestrates the full analysis chain for single-sample ancestry inference.

### Pipeline Stages

#### 1. korusEncode (internal call)
- Builds a k-mer Bloom filter–derived feature matrix from input reads
- Generates the sample × genomic window representation used for inference

#### 2. korusDecode.py
- Performs ancestry inference using one of:
  - hard supervision
  - soft supervision
  - unsupervised (KMeans-based initialization)

#### 3. plot_ntroot_lai.py (ntRoot dependency)
- Generates interactive local ancestry inference (LAI) visualizations

#### 4. HTML report generation
- Produces a combined PCA + LAI report:
  - `*.report.html`
- Supported layouts:
  - iframe mode (default)
  - embedded mode (`--report_mode embed`)

---

### Modes

#### `--mode hard`
- Uses discrete reference labels (`--labels`)
- Performs supervised ancestry classification

#### `--mode soft`
- Uses fractional ancestry supervision (`--fractions`)
- Supports admixed reference modeling

#### `--mode unsupervised`
- Partitions the reference panel into --c K-means clusters and uses the resulting cluster assignments as pseudo-labels for ancestry inference 
- Operates without explicit ancestry labels

---

### Output Files

For a single-sample run, the pipeline produces:

- `*.lai.tsv` — local ancestry inference results
- `*.pca3d.html` — population context (PCA embedding)
- `*.lai.html` — interactive genome view
- `*.report.html` — combined PCA + LAI interactive visualization report (see below examples)


---

### Dependency Note

The pipeline requires:

- `plot_ntroot_lai.py` (from ntRoot)

Available from:
- https://github.com/BirolLab/ntroot
- Conda:
  ```bash
  conda install -c bioconda -c conda-forge ntroot
  ```

---

## Interactive Visualizations <a name="viz"></a>

As part of `korusDecode.py`, KORUS generates interactive PCA visualizations that position the query sample in relation to reference populations in ancestry space. Below visualization examples show sample ERR3988858 from Puerto Rico.

### PCA Visualizations

These HTML files allow exploration of population structure under different inference modes:

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

## Demo <a name="demo"></a>

```bash
./korusMatrix.py --help
./korusDecode.py --help
```

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
