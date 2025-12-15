#!/bin/bash
# 提取小鼠基因体上下游500kb范围
# 输入: gencode.vM38.annotation.gff3.gz
# 输出: mm10_genebody_500kb.bed

GFF3="gencode.vM38.annotation.gff3.gz"
OUTPUT="mm10_genebody_500kb.bed"

echo "正在处理小鼠 GENCODE GFF3 文件..."

# 使用 awk 提取基因信息并计算500kb范围
zcat "$GFF3" | awk '
BEGIN {
    OFS = "\t"
    # 定义扩展范围（500kb）
    extension = 500000
}

# 只处理基因行
$3 == "gene" {
    # 从属性字段提取基因名称
    gene_name = "unknown"
    if (match($9, /gene_name=([^;]+)/)) {
        gene_name = substr($9, RSTART+10, RLENGTH-10)
        # 清理可能的引号
        gsub(/"/, "", gene_name)
    } else if (match($9, /gene_id=([^;]+)/)) {
        gene_name = substr($9, RSTART+8, RLENGTH-8)
        gsub(/"/, "", gene_name)
    }
    
    # 计算扩展后的区域
    start_pos = $4 - extension
    end_pos = $5 + extension
    
    # 确保起始位置不小于1
    if (start_pos < 1) start_pos = 1
    
    # 输出BED格式：染色体、起始、结束、基因名
    # 注意：GFF3是1-based，BED是0-based，所以起始位置要减1
    print $1, start_pos-1, end_pos, gene_name
}' | sort -k1,1 -k2,2n > "$OUTPUT"

echo "完成！输出文件: $OUTPUT"
echo "生成的区域数量: $(wc -l < $OUTPUT)"
