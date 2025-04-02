#!/bin/bash


atac_peaks="/home/dell/data/SCAFE/ISSAAC_5p/data/ATAC_5pRNA/ZKBL-20241228-L-01-2025-01-041610/Sample_JZ24252890-241225-ISSAAC-5P-6N-ATAC-241225-ISSAAC-5P-6N-ATAC/macs2_pk/aggregate_peaks.narrowPeak"
output_dir="overlap_results2"


mkdir -p "$output_dir"


sorted_atac_peaks="${output_dir}/sorted_ATAC_peaks.bed"
sort -k1,1 -k2,2n -k3,3n "$atac_peaks" | uniq > "$sorted_atac_peaks"


for bed6_file in *.bed6; do
    echo "Processing $bed6_file ..."
    
   
    sorted_bed6="${output_dir}/${bed6_file%.bed6}_sorted.bed"
    sort -k1,1 -k2,2n -k3,3n "$bed6_file" | uniq > "$sorted_bed6"
    
    # 执行交集操作，使用 -u 确保唯一性
    bedtools intersect \
        -a "$sorted_bed6" \
        -b "$sorted_atac_peaks" \
        -u \
        > "${output_dir}/${bed6_file%.bed6}_overlap.bed"
done

echo "All overlaps completed! Results saved to: $output_dir"
