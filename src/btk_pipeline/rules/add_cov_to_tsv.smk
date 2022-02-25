rule add_cov_to_tsv:
    """
    Add coverage values to tsv.
    """
    input:
        bed = expand("{{assembly}}.{sra}.regions.bed.gz", sra=reads_by_prefix(config).keys()),
        tsv = "%s/{assembly}.chunk_stats.tsv" % busco_path
    output:
        tsv = "{assembly}.chunk_stats.tsv"
    params:
        bedfiles = gzipped_bed_cols(config)
    threads: 4
    log:
        "logs/{assembly}/add_cov_to_tsv.log"
    benchmark:
        "logs/{assembly}/add_cov_to_tsv.benchmark.txt"
    shell:
        """(paste {input.tsv} {params.bedfiles} > {output.tsv}) 2> {log}"""