rule combine_stats:
    """
    Combine sequence and coverage statistics.
    """
    input:
        bedgz = expand("{{assembly}}.{sra}.regions.bed.gz", sra=reads_by_prefix(config).keys()),
        bed = "{assembly}.fasta.bed",
    output:
        "{assembly}.stats.bed"
    threads: 4
    log:
        "logs/{assembly}/combine_stats.log"
    benchmark:
        "logs/{assembly}/combine_stats.benchmark.txt"
    shell:
        """(cp {input.bed} {output}; \
        for BEDGZ in {input.bedgz}; do \
            awk '{print $1 $2 $3 $4}' >> {output};
        done) 2> {log}"""
