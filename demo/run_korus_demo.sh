#1/bin/bash

# KORUS demo
set -eux -o pipefail

# Prep - download required files using curl if not present in working directory

if [ ! -f "ERR3240246.grch38.aln.chr22.fq.gz" ]; then
    curl -O https://zenodo.org/records/23964981/files/ERR3240246.grch38.aln.chr22.fq.gz
fi

if [ ! -f "ERR3242841.grch38.aln.chr22.fq.gz" ]; then
    curl -O https://zenodo.org/records/23964981/files/ERR3242841.grch38.aln.chr22.fq.gz
fi

if [ ! -f "center_kmers.k55ivc.chr22.txt.gz" ]; then
    curl -O https://zenodo.org/records/23964981/files/center_kmers.k55ivc.chr22.txt.gz
fi

if [ ! -f "korusModelMatrixk55.chr22.tsv.gz" ]; then
    curl -O https://zenodo.org/records/23964981/files/korusModelMatrixk55.chr22.tsv.gz
fi

# 1: Running single sample using korus driver script (--mode hard)

korus --reads ERR3240246.grch38.aln.chr22 --kmer_file center_kmers.k55ivc.chr22.txt.gz \
    -k 55 --outprefix ERR3240246_chr22_korus --mode hard --ref korusModelMatrixk55.chr22.tsv.gz

python3 get_gai_labels.py ERR3240246_chr22_korus.lai_summary.tsv > ERR3240246_chr22_korus.gai_labels.txt

gai=$(cat ERR3240246_chr22_korus.gai_labels.txt |cut -f2)
if [ "$gai" != "EUR" ]; then
    echo "Unexpected ancestry GAI for ERR3240246: $gai (expected: EUR) - please check your KORUS installation"
    exit 1
fi
rm ERR3240246_chr22_korus.gai_labels.txt

# 2: Running KORUS using multiple samples (--mode soft)
 
# a. Run korusEncode
korusEncode --reads ERR3240246.grch38.aln.chr22 --kmer_file center_kmers.k55ivc.chr22.txt.gz -k 55 -t4
korusEncode --reads ERR3242841.grch38.aln.chr22 --kmer_file center_kmers.k55ivc.chr22.txt.gz -k 55 -t4

# b. Combine matrices into a single matrix file
korusAggregateMatrices.py ERR3240246.grch38.aln.chr22_k55_t500000_korus-matrix.tsv \
    ERR3242841.grch38.aln.chr22_k55_t500000_korus-matrix.tsv -o multi-sample.korus-matrix.tsv

# c. Run korusDecode.py 
korusDecode.py --ref korusModelMatrixk55.chr22.tsv.gz \
    --sample multi-sample.korus-matrix.tsv --fractions ../model/korusModelLabelsFractions.tsv.gz \
    --outprefix multi-sample.korus

python3 get_gai_labels.py  multi-sample.korus.lai_summary.tsv > multi-sample.korus.lai_summary.gai_labels.txt

gai1=$(cat multi-sample.korus.lai_summary.gai_labels.txt |grep ERR3240246 |cut -f2)
gai2=$(cat multi-sample.korus.lai_summary.gai_labels.txt |grep ERR3242841 |cut -f2)
if [ "$gai1" != "EUR" ]; then
    echo "Unexpected ancestry GAI for ERR3240246: $gai1 (expected: EUR) - please check your KORUS installation"
    exit 1
fi

if [ "$gai2" != "SAS" ]; then
    echo "Unexpected ancestry GAI for ERR3242841: $gai2 (expected: SAS) - please check your KORUS installation"
    exit 1
fi

rm multi-sample.korus.lai_summary.gai_labels.txt

# 3: Running KORUS using multiple samples (unsupervised mode) - only to check for runtime errors
korusDecode.py --ref korusModelMatrixk55.chr22.tsv.gz \
    --sample multi-sample.korus-matrix.tsv \
    --outprefix multi-sample.korus_unsupervised

rm multi-sample.korus_unsupervised*

echo "KORUS demo completed successfully!"