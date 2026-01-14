configfile: 'config.json'

rule all:
    input:
        'Sample_JZ24252890-241225-ISSAAC-5P-6N-ATAC-241225-ISSAAC-5P-6N-ATAC/ATA"
                             1C_filter/split_R1.fastq.gz',
        'Sample_JZ24252890-241225-ISSAAC-5P-6N-ATAC-241225-ISSAAC-5P-6N-ATAC/ATAC_filter/split_R2.fastq.gz',
        'Sample_JZ24252890-241225-ISSAAC-5P-6N-ATAC-241225-ISSAAC-5P-6N-ATAC/ATAC_filter/split_R3.fastq.gz',
        'genome.sizes',
        'outs/fragments.tsv.gz',
        'outs/fragments.tsv.gz.tbi',
        'outs/per_cell_barcode_total_fragment_count.tsv',
        'outs/fragment_size_distribution.tsv',
        'outs/reads.bed.gz',
        'macs2_pk/aggregate_peaks.narrowPeak',
        'macs2_pk/aggregate_peaks_formatted.bed',
        'outs/peak_read_ov.tsv.gz',
        'outs/raw_mtx/barcodes.tsv',
        'outs/raw_mtx/peaks.bed',
        'outs/raw_mtx/matrix.mtx',
        'outs/raw_mtx/features.tsv',
        'outs/filtered_mtx/barcodes.tsv',
        'outs/filtered_mtx/matrix.mtx',
        'outs/filtered_mtx/features.tsv',
        'outs/filtered_mtx/metrics.csv'


rule getGenomeSize:
    input:
        config['chromap_fa']
    output:
        'genome.sizes'
    shell:
        ''' ~/Downloads/faSize -detailed {input} | sort -k1,1 -k2,2n > {output}
        '''

rule chromap:
    input:
        cb='Sample_JZ24252890-241225-ISSAAC-5P-6N-ATAC-241225-ISSAAC-5P-6N-ATAC/ATAC_filter/split_R2.fastq.gz',
        r1='Sample_JZ24252890-241225-ISSAAC-5P-6N-ATAC-241225-ISSAAC-5P-6N-ATAC/ATAC_filter/split_R1.fastq.gz',
        r2='Sample_JZ24252890-241225-ISSAAC-5P-6N-ATAC-241225-ISSAAC-5P-6N-ATAC/ATAC_filter/split_R3.fastq.gz',
        wl='/home/dell/reference/737K-cratac-v1_rc.txt'
    output:
        'outs/fragments.tsv.gz',
        'outs/fragments.tsv.gz.tbi'
    log:
        'logs/chromap.err'
    params:
        idx=config['chromap_idx'],
        fa=config['chromap_fa']
    threads: 40
    shell:
        ''' chromap -t {threads} --preset atac -x {params.idx} -r {params.fa} \
            -1 {input.r1} -2 {input.r2} -b {input.cb} \
            --barcode-whitelist {input.wl} \
            -o outs/fragments.tsv 2> {log}
            bgzip outs/fragments.tsv && \
            tabix -p bed outs/fragments.tsv.gz
        '''

rule fragToReads:
    input:
        frag='outs/fragments.tsv.gz',
        gs='genome.sizes'
    output:
        met='outs/per_cell_barcode_total_fragment_count.tsv',
        bed='outs/reads.bed.gz'
    shell:
        ''' ./convert_frag_to_reads.py | \
            ~/Downloads/bedClip stdin {input.gs} stdout | \
            sort -k1,1 -k2,2n | \
            gzip > {output.bed}
        '''

rule peakCalling:
    input:
        'outs/reads.bed.gz'
    output:
        'macs2_pk/aggregate_peaks.narrowPeak'
    log:
        'logs/macs2.err'
    params:
        gs=config['macs2_gsize'],
        broad=config['macs2_bpk'],
        fmt=config['macs2_format'],
        shift=config['macs2_shift']
    shell:
        ''' macs2 callpeak -t {input} -g {params.gs} {params.broad} \
            -f {params.fmt} {params.shift} --keep-dup all \
            -B --SPMR --outdir macs2_pk -n aggregate \
            2> {log}
        '''

rule formatPeak:
    input:
        npk='macs2_pk/aggregate_peaks.narrowPeak',
        bdt=config['bedtools_bin'],
        bl=config['blacklist']
    output:
        'macs2_pk/aggregate_peaks_formatted.bed'
    shell:
        ''' cut -f 1-4 {input.npk} | \
            sed '/chrM/d' | \
            {input.bdt} intersect -a - -b {input.bl} -v | \
            sort -k1,1 -k2,2n > {output}
        '''

rule prepareReadPeaksOVs:
    input:
        bdt=config['bedtools_bin'],
        reads='outs/reads.bed.gz',
        peaks='macs2_pk/aggregate_peaks_formatted.bed',
        gsize='genome.sizes'
    output:
        'outs/peak_read_ov.tsv.gz'
    shell:
        ''' {input.bdt} intersect -a {input.peaks} -b {input.reads} \
            -wo -sorted -g {input.gsize} | sort -k8,8 | \
            {input.bdt} groupby -g 8 -c 4 -o freqdesc | gzip > {output}
        '''

rule generate10xStyleOutput1:
    input:
         peaks='macs2_pk/aggregate_peaks_formatted.bed',
         ov='outs/peak_read_ov.tsv.gz'
    output:
        bc='outs/raw_mtx/barcodes.tsv',
        pk='outs/raw_mtx/peaks.bed'
    shell:
        ''' cut -f 1-3 {input.peaks} > {output.pk}
            zcat {input.ov} | cut -f 1 > {output.bc}
        '''

rule generate10xStyleOutput2:
    input:
        ov='outs/peak_read_ov.tsv.gz',
        bc='outs/raw_mtx/barcodes.tsv',
        peaks='macs2_pk/aggregate_peaks_formatted.bed'
    output:
        'outs/raw_mtx/matrix.mtx'
    script:
        './generate_csc_mtx.py'

rule prepareSoloFilter:
    input:
        'outs/raw_mtx/peaks.bed'
    output:
        'outs/raw_mtx/features.tsv'
    shell:
        ''' awk 'BEGIN{{OFS="\t"}}{{print $1 "-" $2 "-" $3, $1 "-" $2 "-" $3}}' \
            {input} > {output}
        '''

rule soloFilterCell:
    input:
        pg=config['star_bin'],
        bc='outs/raw_mtx/barcodes.tsv',
        ct='outs/raw_mtx/matrix.mtx',
        ft='outs/raw_mtx/features.tsv'
    output:
        'outs/filtered_mtx/barcodes.tsv',
        'outs/filtered_mtx/matrix.mtx',
        'outs/filtered_mtx/features.tsv'
    shell:
        ''' {input.pg} --runMode soloCellFiltering \
            outs/raw_mtx outs/filtered_mtx/ \
            --soloCellFilter EmptyDrops_CR
        '''

rule getFragmentSizeDistributionChromap:
    input:
        'outs/fragments.tsv.gz'
    output:
        'outs/fragment_size_distribution.tsv'
    shell:
        ''' echo -e "size\tcount" > {output}
            zcat {input} | sed '/^#/d' | awk '{{print $3-$2}}' | \
            sort | uniq -c | sort -b -k2,2n | \
            awk 'BEGIN{{OFS="\t"}}{{print $2, $1}}' >> {output}
        '''

rule getMetricsChromap:
    input:
        bc='outs/filtered_mtx/barcodes.tsv',
        mtx='outs/filtered_mtx/matrix.mtx'
    output:
        'outs/filtered_mtx/metrics.csv'
    shell:
        ''' ./output_qc_mtx {input.mtx} {input.bc} {output}
        '''
