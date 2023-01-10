rule run_mosdepth:
    """
    Run Mosdepth to get coverage depth.
    """
    input:
        bam = "%s/{assembly}.{sra}.bam" % minimap_path,
        csi = "%s/{assembly}.{sra}.bam.csi" % minimap_path,
        bed = "%s/{assembly}.chunk_stats.mask.bed" % chunk_stats_path
    output:
        gz = "{assembly}.{sra}.regions.bed.gz"
    params:
        prefix = lambda wc: "%s.%s" % (wc.assembly, wc.sra)
    threads: 4
    log:
        "logs/{assembly}/run_mosdepth/{sra}.log"
    benchmark:
        "logs/{assembly}/run_mosdepth/{sra}.benchmark.txt"
    shell:
        """mosdepth -n -x -b {input.bed} \
                 -t 4 {params.prefix} {input.bam} 2> {log}"""