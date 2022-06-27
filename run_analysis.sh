#!/bin/bash

python3 ../code/nanoCT_preprocess.py --fastq ../../bcd_CT/single-cell/fastq/P24004/P24004_1001/02-FASTQ/211221_A00621_0569_BHTNK3DRXY/ \
  --barcodes ATAGAGGC TATAGCCT CCTATCCT --modalities ATAC H3K27ac H3K27me3 \
  --cellranger_ref /data/ref/cellranger-atac/refdata-cellranger-atac-mm10-2020-A-2.0.0/ \
  --threads 40 --tempdir /data/proj/GCB_MB/tmp --snakeargs='--profile htcondor'