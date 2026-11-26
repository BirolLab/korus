#1/bin/bash

# KORUS demo
set -eux -o pipefail

# Prep - download required files using curl if not present in working directory

if [ ! -f "ERR3240246.grch38.aln.chr22.fq.gz" ]; then
    curl -O https://zenodo.org/records/24564514/files/ERR3240246.grch38.aln.chr22.fq.gz
fi

if [ ! -f "ERR3242841.grch38.aln.chr22.fq.gz" ]; then
    curl -O https://zenodo.org/records/24564514/files/ERR3242841.grch38.aln.chr22.fq.gz
fi

if [ ! -f "HG02132_mat.chr20.fa.gz" ]; then
    curl -O https://zenodo.org/records/24564514/files/HG02132_mat.chr20.fa.gz
fi

if [ ! -f "HG02132_pat.chr20.fa.gz" ]; then
    curl -O https://zenodo.org/records/24564514/files/HG02132_pat.chr20.fa.gz
fi

if [ ! -f "center_kmers.k55ivc.chr22.txt.gz" ]; then
    curl -O https://zenodo.org/records/24564514/files/center_kmers.k55ivc.chr22.txt.gz
fi

if [ ! -f "center_kmers.k55ivc.chr20.txt.gz" ]; then
    curl -O https://zenodo.org/records/24564514/files/center_kmers.k55ivc.chr20.txt.gz
fi

if [ ! -f "korusModelMatrixk55.chr22.tsv.gz" ]; then
    curl -O https://zenodo.org/records/24564514/files/korusModelMatrixk55.chr22.tsv.gz
fi

if [ ! -f "korusModelMatrixk55.chr20.tsv.gz" ]; then
    curl -O https://zenodo.org/records/24564514/files/korusModelMatrixk55.chr20.tsv.gz
fi

# 1: Running single reads sample using korus driver script (--mode hard)

korus --reads ERR3240246.grch38.aln.chr22 --kmer_file center_kmers.k55ivc.chr22.txt.gz \
    -k 55 --outprefix ERR3240246_chr22_korus --mode hard --ref korusModelMatrixk55.chr22.tsv.gz

python3 get_gai_labels.py ERR3240246_chr22_korus.lai_summary.tsv > ERR3240246_chr22_korus.gai_labels.txt

gai=$(cat ERR3240246_chr22_korus.gai_labels.txt |cut -f2)
if [ "$gai" != "EUR" ]; then
    echo "Unexpected ancestry GAI for ERR3240246: $gai (expected: EUR) - please check your KORUS installation"
    exit 1
fi
rm ERR3240246_chr22_korus.gai_labels.txt

# 2: Running KORUS using multiple reads samples (--mode soft)
 
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

# 4: Running KORUS in genome mode (--hard mode)
korus --genome HG02132_mat.chr20.fa.gz HG02132_pat.chr20.fa.gz \
    --kmer_file center_kmers.k55ivc.chr20.txt.gz -k 55 -t 4 \
    --mode hard --outprefix HG02132_korus --ref korusModelMatrixk55.chr20.tsv.gz

python3 get_gai_labels.py HG02132_korus.lai_summary.tsv > HG02132_korus.gai_labels.txt

gai=$(cat HG02132_korus.gai_labels.txt |cut -f2)
if [ "$gai" != "EAS" ]; then
    echo "Unexpected ancestry GAI for HG02132: $gai (expected: EAS) - please check your KORUS installation"
    exit 1
fi
rm HG02132_korus.gai_labels.txt

# 5: Running KORUS in genome mode (--soft), using previously built matrix
korusDecode.py --ref korusModelMatrixk55.chr20.tsv.gz --sample HG02132_korus_k55_t500000_korus-matrix.tsv \
    --fractions ../model/korusModelLabelsFractions.tsv.gz --outprefix HG02132_korus_soft -t 4

python3 get_gai_labels.py HG02132_korus_soft.lai_summary.tsv > HG02132_korus_soft.gai_labels.txt

gai=$(cat HG02132_korus_soft.gai_labels.txt |cut -f2)
if [ "$gai" != "EAS" ]; then
    echo "Unexpected ancestry GAI for HG02132: $gai (expected: EAS) - please check your KORUS installation"
    exit 1
fi
rm HG02132_korus_soft.gai_labels.txt

# 6: Running KORUS in genome mode (unsupervised) - only to check for runtime errors
korusDecode.py --ref korusModelMatrixk55.chr20.tsv.gz --sample HG02132_korus_k55_t500000_korus-matrix.tsv \
    --outprefix HG02132_korus_unsupervised -t 4

echo "KORUS demo completed successfully!"
