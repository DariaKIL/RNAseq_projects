#!/bin/bash

# List of SRR IDs
samples=(
"SRR12649401"
"SRR12649399"
"SRR12649403"
"SRR12649402"
"SRR12649404")

# Downloading using prefetch и converting to FASTQ using fastq-dump
for sample in "${samples[@]}"
do
    echo "Downloading $sample"
    if prefetch "$sample"; then
        echo "Successfully downloaded $sample"
    else
        echo "Failed to download $sample" >&2
        continue
    fi
    
    echo "Converting $sample to FASTQ.gz"
    if fastq-dump --split-files --gzip "$sample"; then
        echo "Successfully converted $sample to FASTQ"
    else
        echo "Failed to convert $sample to FASTQ.gz" >&2
        continue
    fi
done

echo "All samples processed"
