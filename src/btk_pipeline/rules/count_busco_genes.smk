rule count_busco_genes:
    """
    Count busco genes.
    """
    input:
        busco = expand("{{assembly}}.busco.{lineage}/full_table.tsv.gz", lineage=config["busco"]["lineages"]),
        mask = "%s/{assembly}.chunk_stats.tsv" % chunk_stats_path
    output:
        tsv = "{assembly}.chunk_stats.tsv",
    params:
        chunk = set_stats_chunk(config),
        overlap = 0,
        max_chunks = 1000000000,
        min_length = set_stats_chunk(config),
    threads: 1
    log:
        "logs/{assembly}/count_busco_genes.log"
    benchmark:
        "logs/{assembly}/count_busco_genes.benchmark.txt"
    script:
        "../lib/count_busco_genes.py"
