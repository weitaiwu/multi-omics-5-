for barcode_file in cluster_*_barcodes.txt; do
  cluster=${barcode_file%_*}
  cluster=${cluster##*_}

  subset-bam \
    --bam /home/dell/data/SCAFE/ISSAAC_5p/data/ATAC_5pRNA/ZKBL-20241228-L-01-2025-01-141457/Sample_JZ24252891-241225-ISSAAC-5P-6N-RNA-241225-ISSAAC-5P-6N-RNA/star_50/filtered_6N.bam \
    --bam-tag CB \
    --cell-barcodes "$barcode_file" \
    --out-bam "cluster_${cluster}.bam" \
    --cores 4
done
