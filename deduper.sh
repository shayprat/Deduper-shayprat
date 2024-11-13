#!/bin/bash

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --mem=16G
#SBATCH --job-name=deduper
#SBATCH --error=log/dedupe_%j.err
#SBATCH --output=log/dedupe_%j.out


/usr/bin/time -v /projects/bgmp/spratap/bioinfo/Bi624/Deduper-shayprat/pratap_deduper.py \
    -f /projects/bgmp/spratap/bioinfo/Bi624/Deduper-shayprat/sorted_C1_SE_uniqAlign.sam \
    -o deduped_sorted_C1_SE_uniqAlign.sam \
    -u /projects/bgmp/spratap/bioinfo/Bi624/Deduper-shayprat/STL96.txt