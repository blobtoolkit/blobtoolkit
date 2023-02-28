rule run_blobtk_depth:
    """
    Run blobtk depth to get coverage depth.
    """
    input:
        bam = "%s/{assembly}.{sra}.bam" % minimap_path,
        csi = "%s/{assembly}.{sra}.bam.csi" % minimap_path,
    output:
        bed = "{assembly}.{sra}.regions.bed.gz"
    params:
        prefix = lambda wc: "%s.%s" % (wc.assembly, wc.sra)
    threads: 4
    log:
        "logs/{assembly}/run_blobtk_depth/{sra}.log"
    benchmark:
        "logs/{assembly}/run_blobtk_depth/{sra}.benchmark.txt"
    shell:
        """blobtk depth -b {input.bam} \
                 -s 1000 -O {output.bed} 2> {log}"""