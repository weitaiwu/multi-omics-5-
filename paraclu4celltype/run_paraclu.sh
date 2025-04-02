#!/bin/bash

for counts_file in *_counts.txt; do
    # 提取文件名中的唯一标识符（假设格式为 XXX_数字_counts.txt）
    base_name="${counts_file%_counts.txt}"  # 移除末尾的 "_counts.txt"
    identifier="${base_name##*_}"           # 提取最后一个下划线后的部分（如 0）
    
    # 生成输出文件名（如 paraclu_0.output）
    output_file="paraclu_${identifier}.output"
    
    # 执行 paraclu 命令
    echo "Running: paraclu 1 $counts_file > $output_file"
    paraclu 1 "$counts_file" > "$output_file"
done

echo "All paraclu jobs completed!"
