import sys

tso = "CTAACGGG"  # your TSO序列
max_mismatch = 1

complement = str.maketrans("ACGT", "TGCA")
tso_rc = tso.translate(complement)[::-1]  # CCCGTTAG

def hamming(s1, s2):
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

for line in sys.stdin:
    line = line.strip()
    
    if line.startswith('@'):
        print(line)
        continue
    
    fields = line.split('\t')
    if len(fields) < 10:
        continue  
        
    seq = fields[9]
    seq_len = len(seq)
    found = False
    
    # 检查正向25-32bp区域
    if seq_len >= 32:
        forward_sub = seq[24:32]
        if len(forward_sub) == 8:
            if hamming(forward_sub, tso) <= max_mismatch:
                found = True
    
    # 检查倒数25-32bp区域
    if not found and seq_len >= 32:
        rev_sub = seq[seq_len-32 : seq_len-24]
        if len(rev_sub) == 8:
            if hamming(rev_sub, tso_rc) <= max_mismatch:
                found = True
    
    if found:
        print(line)
EOF
