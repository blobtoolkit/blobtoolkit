rule bamtools_stats:
    """
    Run bamtools stats to generate summary statistics for each BAM file
    """
    input:
        bam = "%s/{assembly}.{sra}.bam" % minimap_path
    output:
        "{assembly}.{sra}.bam.stats"
    threads: 1
    log:
        "logs/{assembly}/bamtools_stats/{sra}.log"
    benchmark:
        "logs/{assembly}/bamtools_stats/{sra}.benchmark.txt"
    shell:
        """bamtools stats \
            -in {input.bam} \
            -insert \
            > {output} 2> {log}"""