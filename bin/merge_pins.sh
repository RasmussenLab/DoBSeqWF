#!/bin/bash
#
# mads
# script for merging two-dimensional vcfs
#

mkdir -p .tmp output

# Create annotation header file
cat <<EOL > .tmp/annots.hdr
##FORMAT=<ID=POOL_ID,Number=1,Type=String,Description="Pool ID">
##INFO=<ID=SAMPLE_ID,Number=1,Type=String,Description="Sample ID">
EOL

for i in *.vcf.gz; do 
    echo "Annotating $i..."
    
    # Extract sample IDs and assign row and column IDs
    pool_ids=$(bcftools query -l "$i")
    row_id=$(echo "$pool_ids" | awk 'NR==1')  # First sample ID
    col_id=$(echo "$pool_ids" | awk 'NR==2')  # Second sample ID
    
    # Extract sample_id from the filename
    sample_id=$(basename "${i%_unique_pins.vcf.gz}")
    echo "Sample ID: $sample_id"
    echo "Row ID: $row_id"
    echo "Column ID: $col_id"
    
    # Query basic VCF data
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$i" > .tmp/tmp.tsv

    # Add annotations to TSV
    awk -F'\t' -v sample_id="$sample_id" -v row_id="$row_id" -v column_id="$col_id" \
        '{ print $0 "\t" sample_id "\t" row_id "\t" column_id }' .tmp/tmp.tsv > .tmp/tmp_with_ids.tsv
    
    # Compress and index the annotated TSV
    bgzip -f .tmp/tmp_with_ids.tsv
    tabix -f -b 2 -s 1 -e 2 .tmp/tmp_with_ids.tsv.gz
    
    # Annotate the VCF with new fields
    bcftools annotate \
        -s "${row_id}" \
        -a .tmp/tmp_with_ids.tsv.gz \
        -c CHROM,POS,REF,ALT,-,FMT/POOL_ID,- \
        -h .tmp/annots.hdr \
        "$i" > .tmp/tmp_wrow.vcf.gz
    
    bcftools annotate \
        -s "${col_id}" \
        -a .tmp/tmp_with_ids.tsv.gz \
        -c CHROM,POS,REF,ALT,-,-,FMT/POOL_ID \
        -h .tmp/annots.hdr \
        .tmp/tmp_wrow.vcf.gz > .tmp/tmp_wcol.vcf.gz
    
    bcftools annotate \
        -a .tmp/tmp_with_ids.tsv.gz \
        -c CHROM,POS,REF,ALT,SAMPLE_ID,-,- \
        -h .tmp/annots.hdr \
        .tmp/tmp_wcol.vcf.gz > .tmp/tmp_wsample.vcf
    
    echo -e "$row_id\tROW\n$col_id\tCOLUMN" > .tmp/rename_samples.txt
    bcftools reheader --samples .tmp/rename_samples.txt .tmp/tmp_wsample.vcf > output/$(basename ${i%_unique_pins.vcf.gz}).vcf    
done

for f in output/*.vcf; do bgzip $f; tabix $f.gz; done
bcftools concat -a output/*.vcf.gz -o unsorted.vcf
bcftools sort --temp-dir . unsorted.vcf -o pinpointables.vcf