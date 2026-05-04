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

KORUS is a two-stage computational framework for ancestry inference from sequencing data using k-mer–derived genomic features encoded in Bloom filters.

Rather than relying on alignment or variant calling, KORUS projects sequencing data into a structured genomic feature space and performs probabilistic decoding to infer ancestry structure across the genome.

The framework is composed of:

-MATRIX stage: construction of k-mer–derived genomic feature matrices
-DECODE stage: probabilistic inference of ancestry from feature space
-Methodological Framework
-MATRIX (korusMatrix.py)

### The MATRIX stage transforms Bloom filter–encoded sequencing data into a structured genomic feature representation.

Key steps:

-Query per-sample k-mer Bloom filters (ntEdit/ntStat compatible)
-Map k-mer hits to ancestry-informative SNVs
-Compute allele-specific support (REF vs ALT k-mers)
-Derive normalized directional signal per locus
-Aggregate signals into fixed genomic windows
-Construct a sample × window feature matrix

Output: continuous genomic feature matrix across samples and windows

### DECODE (korusDecode.py)

The DECODE stage performs probabilistic inference of ancestry from the feature matrix.

Core steps:

-Standardization of the feature matrix
-Dimensionality reduction (PCA) into latent space
-Distance-based classification (e.g., Mahalanobis)
-Probability calibration (softmax)
-Optional unsupervised clustering (KMeans)
-Spatial smoothing via HMM / Viterbi decoding

Output:

-Global ancestry fractions
-Local ancestry assignments
-Posterior probability matrices
-Smoothed ancestry trajectories across the genome
-Conceptual Model

Sequencing Data
> k-mer Bloom Filter Space
> MATRIX (windowed genomic feature representation)
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

### 1. MATRIX stage

```bash
korusMatrix.py   --bf_dir /path/to/bloom_filters   --snv_kmers snv_kmers.tsv   --output_matrix signal_matrix.tsv   --window_size 500000   --k 55
```

### 2. DECODE stage

```bash
korusDecode.py   --ref reference_matrix.tsv   --sample signal_matrix.tsv   --labels reference_labels.tsv   --k 5   --outprefix korus_run
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
./korusMatrix.py --help
./korusDecode.py --help
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
