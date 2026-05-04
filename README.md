[![ntRoot dependency](https://img.shields.io/badge/depends-ntRoot%20framework-brightgreen)](https://github.com/birollab/ntRoot)
[![Python](https://img.shields.io/badge/python-3.9%2B-blue)]()
[![License](https://img.shields.io/badge/license-GPLv3-lightgrey)]()

![Logo](https://github.com/warrenlr/korus/blob/main/korus-logo.png)

# KORUS

K-mer Orthogonal Representation for Unified Stratification (of ancestry)  
K-mer Orthogonal Representation of Unseen Structure  

## Signal > Decode Framework for Genomic Ancestry Inference

---

## Contents

1. Credit  
2. Overview  
3. Methodological Framework  
4. Installation  
5. Dependencies  
6. Usage  
7. Inputs  
8. Outputs  
9. Workflow overview  
10. Demo  
11. Citing  
12. License  

---

## Credit

Developed by René L Warren.

---

## Overview

KORUS is a two-stage computational framework for ancestry inference from sequencing data using k-mer–derived genomic signals encoded in Bloom filters.

Instead of relying directly on alignment or variant calling, KORUS transforms sequencing data into a structured genomic signal field and performs probabilistic decoding to infer ancestry structure across the genome.

The framework is composed of:

- **SIGNAL stage**: construction of k-mer–derived genomic signal matrices  
- **DECODE stage**: probabilistic inference of ancestry from signal space  

---

## Methodological Framework

### SIGNAL (korusSIGNAL.py)

The SIGNAL stage converts raw Bloom-filter–encoded sequencing data into a structured genomic signal representation.

Key steps:
- Query per-sample k-mer Bloom filters (ntEdit/ntStat compatible)
- Map k-mer hits to ancestry-informative SNVs
- Compute allele-specific support (REF vs ALT k-mers)
- Normalize directional signal
- Aggregate signals into fixed genomic windows
- Output sample × window signal matrix

Output: continuous genomic signal matrix per sample.

---

### DECODE (korusDECODE.py)

The DECODE stage performs probabilistic inference of ancestry from the signal matrix.

Core steps:
- Standardization of signal matrix
- PCA embedding into latent space
- Gaussian / Mahalanobis classification
- Softmax probability calibration
- Optional KMeans clustering (unsupervised mode)
- HMM / Viterbi smoothing across windows

Output:
- Global ancestry fractions
- Local ancestry states
- Posterior probability matrices
- Smoothed ancestry trajectories

---

## Conceptual Model

Sequencing Data  
> k-mer Bloom Filter Space  
> SIGNAL (windowed genomic signal field)  
> DECODE (latent ancestry reconstruction)

---

## Installation

```bash
git clone https://github.com/birollab/korus.git
cd korus
chmod +x korusSIGNAL.py korusDECODE.py
```

---

## Dependencies

- Python ≥ 3.9  
- numpy  
- pandas  
- scikit-learn  
- matplotlib  
- btllib  
- ntEdit / ntStat compatible Bloom filter ecosystem  

---

## Usage

### 1. SIGNAL stage

```bash
korusSIGNAL.py   --bf_dir /path/to/bloom_filters   --snv_kmers snv_kmers.tsv   --output_matrix signal_matrix.tsv   --window_size 500000   --k 55
```

### 2. DECODE stage

```bash
korusDECODE.py   --ref reference_matrix.tsv   --sample signal_matrix.tsv   --labels reference_labels.tsv   --k 5   --outprefix korus_run
```

---

## Inputs

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

### Reference labels (optional)
```
sample   ancestry
ERR1     EUR
ERR2     AFR
```

---

## Outputs

### Global outputs
- ancestry_fractions.tsv  
- posterior_probabilities.tsv  

### Local ancestry outputs
- *.lai.tsv  
- *.lai_smoothed.tsv  

### Visualizations
- PCA plots  
- ancestry trajectories  

### Logs
- runtime logs per stage  

---

## Workflow overview

1. Bloom filter interrogation of k-mers  
2. SNV-level allele signal extraction  
3. Window-based aggregation (SIGNAL)  
4. Latent space embedding (PCA)  
5. Probabilistic classification (Gaussian / Mahalanobis)  
6. Optional clustering and HMM smoothing (DECODE)  

---

## Demo

```bash
./korusSIGNAL.py --help
./korusDECODE.py --help
```

---

## Citing

TBD 

---

## License

KORUS is released under the GNU General Public License v3 (GPLv3).

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, version 3.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See http://www.gnu.org/licenses/ for details.
