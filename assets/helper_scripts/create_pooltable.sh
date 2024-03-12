#!/bin/bash

# mads 2024-03-12
# Create pooltable from a directory of fastq files.
#
# Usage: create_pooltable.sh <directory>

set -euo pipefail

TEMP1="create_pooltable_1.TEMP"
TEMP2="create_pooltable_2.TEMP"
OUTFILE="pooltable.tsv"

echo "Creating pooltabletable $(date -Is)"
echo "Searching for files in $*"

# Get all forward read files
find $@ -name '*_1.fastq.gz' -print > ${TEMP1}

echo "Found $(cat ${TEMP1} | wc -l) files ending in _1.fastq.gz"

# Then remove duplicated filenames (basenames). Also in case of different subdirectories.
paste ${TEMP1}  <(cat ${TEMP1} | xargs -n 1 basename)  \
    | sort -k 2 \
    | uniq -f 1 \
    | cut -f 1  \
    > ${TEMP2}

echo "  of which $(cat ${TEMP2} | wc -l) files were not duplicates and thus retained"

# Create full pooltable (poolID, fq1, fq2) based on these unique filenames
paste <(perl -lpe 's|^.+/([^/]+)_1.fastq.gz$|\1|' ${TEMP2}) \
    ${TEMP2} \
    <(sed 's/_1.fastq.gz/_2.fastq.gz/' ${TEMP2}) \
    > ${OUTFILE}

# Clean up temporary files
rm ${TEMP1} ${TEMP2}

echo "Finished with pooltable!"