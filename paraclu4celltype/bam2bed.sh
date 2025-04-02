#!/bin/bash

for bam_file in *.bam; do
    [ ! -f "$bam_file" ] && continue
    output_txt="${bam_file%.bam}_counts.txt"
    echo "Processing $bam_file → $output_txt"
    
    samtools view "$bam_file" | awk '{
        chrom = $3;
        pos = $4;
        flag = $2 + 0;
        
       
        if (int(flag / 16) % 2 == 1) {
            strand = "-";
        } else {
            strand = "+";
        }
        
       
        cigar = $6;
        ref_len = 0;
        remaining = cigar;
        while (remaining != "") {
            if (match(remaining, /^[0-9]+[MDN=X]/)) {
                op_str = substr(remaining, 1, RLENGTH);
                len = substr(op_str, 1, RLENGTH - 1) + 0;
                op = substr(op_str, RLENGTH, 1);
                if (op ~ /[MDN=X]/) {
                    ref_len += len;
                }
                remaining = substr(remaining, RLENGTH + 1);
            } else {
                break;
            }
        }
        
      
        if (strand == "-") {
            start = pos + ref_len - 1;
        } else {
            start = pos;
        }
        
        print chrom "\t" strand "\t" start;
    }' | sort -k1,1 -k2,2 -k3,3n | uniq -c | awk -v OFS="\t" '{print $2, $3, $4, $1}' > "$output_txt"
done

echo "All BAM files processed!"
