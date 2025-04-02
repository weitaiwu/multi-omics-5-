#!/bin/bash

# 遍历所有以 .output 结尾的文件
for output_file in *.output; do
    # 生成对应的 .bed6 文件名（如 paraclu_0.output → paraclu_0.bed6）
    bed6_file="${output_file%.output}.bed6"
    
    # 使用 awk 转换格式
    awk '{print $1 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $2}' "$output_file" > "$bed6_file"
    
    echo "Converted: $output_file → $bed6_file"
done

echo "All conversions completed!"
