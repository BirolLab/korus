[![ntRoot dependency](https://img.shields.io/badge/depends-ntRoot%20framework-brightgreen)](https://github.com/birollab/ntRoot)
[![Python](https://img.shields.io/badge/python-3.9%2B-blue)]()
[![License](https://img.shields.io/badge/license-GPLv3-lightgrey)]()

![Logo](https://github.com/warrenlr/korus/blob/main/korus-logo.png)

# KORUS

K-mer Orthogonal Representation for Unified Stratification (of ancestry) — an orthogonalized ancestry feature representation in PCA-reduced space

KORUS derives window-level ancestry by embedding k-mer-based ALT/REF signals into a feature matrix, assigning ancestry likelihoods via Mahalanobis distance in PCA-reduced feature space followed by softmax normalization, and enforcing genomic continuity using HMM smoothing before computing ancestry fractions.

## Feature Matrix → Probabilistic Decode Framework for Genomic Ancestry Inference

---

## Contents

1. [Synopsis](#synopsis)
2. [Features](#features)
3. [Credit](#credit)  
4. [Overview](#overview)  
5. [Methodological Framework](#methodological-framework)  
6. [Installation](#installation)  
7. [Dependencies](#dependencies)  
8. [Usage](#usage)  
9. [Inputs](#inputs)  
10. [Outputs](#outputs)  
11. [Workflow Overview](#workflow-overview)  
12. [Demo](#demo)  
13. [Citing](#citing)  
14. [License](#license)  

---

## Synopsis <a name="synopsis"></a>

KORUS converts k-mer presence in Bloom filters into per-window scalar signals (ALT–REF normalized), aggregates these into a sample-by-window feature matrix, assigning ancestry likelihoods via Mahalanobis distance in PCA-reduced feature space followed by softmax normalization, and optionally applies an HMM (Viterbi) to smooth window-wise ancestry paths before aggregating them into genome-wide fractions.


### Key methodological features <a name="features"></a>

- k-mer/Bloom-filter ancestry representation
- Reference-free signal extraction
- Local ancestry inference without alignment
- Soft-supervised ancestry manifold modeling

---

## Credit <a name="credit"></a>

Developed by XX, XX and René L Warren.

---

## Overview <a name="overview"></a>

KORUS is a two-stage computational framework for ancestry inference from sequencing data using k-mer–derived genomic features encoded in Bloom filters.

Rather than relying on alignment or variant calling, KORUS transforms sequencing data into a structured genomic feature space and performs probabilistic decoding to infer ancestry structure across the genome. The feature space is treated as an orthogonalized ancestry embedding in PCA-reduced coordinates, enabling covariance-aware distance classification. In other words, the feature space is treated as a decorrelated representation in PCA space prior to Mahalanobis distance estimation.

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

#### Hard supervision (--labels)

Each reference sample is assigned a single discrete ancestry class.

KORUS learns ancestry cluster structure from categorical population labels and performs window-level classification by minimizing Mahalanobis distance to class-specific distributions in PCA space.

This framework assumes each reference genome belongs predominantly to one ancestry class.

In practice this behaves similarly to:
- supervised metric learning
- centroid/covariance-based discriminant assignment
- Gaussian distance classification in reduced-dimensional space

with optional HMM smoothing across genomic windows.


#### Soft supervision (--fractions)

Each reference sample contributes fractional ancestry memberships rather than a single categorical label.

Instead of learning discrete ancestry clusters from hard assignments, KORUS estimates ancestry distributions using weighted contributions from admixed reference genomes.

Mahalanobis distances are therefore computed relative to ancestry-specific distributions informed by continuous ancestry weights.

This allows:
- admixed references
- continuous ancestry structure
- smoother ancestry manifolds
- reduced dependence on discrete continental boundaries

Conceptually, this behaves like:
- weighted Gaussian supervision
- probabilistic prototype learning
- soft-label metric embedding

rather than categorical classification.

#### Hard vs soft supervision behavior

When using --labels, KORUS can model any discrete ancestry category represented in the reference labels.

When using --fractions, the ancestry space is constrained to the ancestry columns provided in the fraction matrix. Ancestries absent from the fraction table are not explicitly modeled during probabilistic decoding.


**Core steps:**
- Standardization of the feature matrix
- Dimensionality reduction (PCA)
- Distance-based classification (Mahalanobis)
- Softmax normalization of ancestry likelihoods (Distance-to-probability transformation)
- Optional unsupervised cluster initialization (KMeans)
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
chmod +x korusMatrix.py korusDecode.py
```

---

## Dependencies

- Python ≥ 3.9  
- numpy  
- pandas  
- scikit-learn  
- matplotlib
- plotly
- btllib
- gzip 
- ntEdit / ntStat compatible Bloom filter ecosystem  

---

## Usage <a name="usage"></a>



### 1. MATRIX stage

```bash
The SNV-to-signal encode runs via the korusEncode snakemake (runs ntCard, ntStat and korusMatrix.py):
usage: korusEncode [-h] --reads READS --kmer_file KMER_FILE -k K [--tile TILE] [-t T] [-n] [-f] [-V]
korusEncode: the following arguments are required: --reads, --kmer_file, -k 


usage: korusMatrix.py [-h] --bf_dir BF_DIR --kmer_file KMER_FILE --output_matrix OUTPUT_MATRIX
                      [--window_size WINDOW_SIZE] [--k K] [-V]
korusMatrix.py: the following arguments are required: --bf_dir, --kmer_file, --output_matrix

./korusMatrix.py 
  --bf_dir /path/to/bloom_filters 
  --snv_kmers snv_kmers.tsv 
  --output_matrix sampleMatrix.tsv 
  --window_size 500000 
  --k 55
```

### 2. DECODE stage

```bash
usage: korusDecode.py [-h] --ref REF --sample SAMPLE [--labels LABELS] [--fractions FRACTIONS]
                          [--k K] --outprefix OUTPREFIX [--window_radius WINDOW_RADIUS]
                          [--require_full_window] [--stay_prob STAY_PROB] [--no_hmm] [-V]
korusDecode.py: the following arguments are required: --ref, --sample, --outprefix

./korusDecode.py 
  [-h]
  --ref korusModelMatrixk55.tsv 
  --sample sampleMatrix.tsv 
  --labels korusModelLabels.tsv
  --fractions korusModelLabelsFractions.tsv.gz
  --k 5 
  --outprefix korus_run
  --stay_prob 0.98
```

---

## Inputs <a name="inputs"></a>

### Bloom filters
Per-sample k-mer Bloom filters compatible with ntEdit/ntStat.

Example:
```
SAMPLE_k55.bf
```

### SNV k-mer file
```
chr  pos  allele  kmer
```

### Reference hard labels (optional, --labels)
```
sample   ancestry
ERR1     EUR
ERR2     AFR
```

### Reference soft labels (optional, --fractions)
```
sample	 EAS	  SAS	   AMR	   AFR	   EUR
ERR1	 0.1961	0.0105	0.7233	0.0595	0.0105
ERR2	 0.0718	0.0298	0.5814	0.2347	0.0823
```
---

## Outputs <a name="outputs"></a>

### Global outputs
- .lai_summary.tsv
- .ancestry_fractions.tsv  
- .posterior_probabilities.tsv  

### Local ancestry outputs
- .SAMPLE.lai.tsv   

### Visualizations
- PCA plots  
- ancestry trajectories  

### Logs
- runtime logs per stage  

---

## Workflow overview <a name="workflow-overview"></a>

1. Bloom filter interrogation
2. SNV-level signal extraction
3. Window aggregation (MATRIX)
4. PCA embedding
5. Probabilistic classification
6. Optional smoothing (DECODE)

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
