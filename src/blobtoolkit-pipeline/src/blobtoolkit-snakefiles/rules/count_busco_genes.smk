rule count_busco_genes:
    """
    Count busco genes.
    """
    input:
        busco = expand("{{assembly}}.busco.{lineage}/full_table.tsv.gz", lineage=config["busco"]["lineages"]),
        mask = "%s/{assembly}.chunk_stats.tsv" % chunk_stats_path
    output:
        tsv = "{assembly}.chunk_stats.tsv"
    params:
        busco = lambda wc: " --in ".join(expand("%s.busco.{lineage}/full_table.tsv.gz" % wc.assembly, lineage=config['busco']['lineages'])),
    threads: 1
    log:
        "logs/{assembly}/count_busco_genes.log"
    benchmark:
        "logs/{assembly}/count_busco_genes.benchmark.txt"
    shell:
        """(btk pipeline count-busco-genes \
            --in {params.busco} \
            --mask {input.mask} \
            --out {output.tsv}) 2> {log}"""
