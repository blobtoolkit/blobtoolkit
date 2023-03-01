rule extract_busco_genes:
    """
    Extract busco genes into a single fasta file.
    """
    input:
        busco = expand("%s/%s.busco.{lineage}/full_table.tsv.gz" % (busco_path, config["assembly"]["prefix"]), lineage=get_basal_lineages(config)),
    output:
        fasta = "{assembly}.busco_genes.fasta"
    params:
        busco = lambda wc: " --busco ".join(expand("%s/%s.busco.{lineage}/full_table.tsv.gz" % (busco_path, wc.assembly), lineage=get_basal_lineages(config))),
    threads: 1
    log:
        "logs/{assembly}/extract_busco_genes.log"
    benchmark:
        "logs/{assembly}/extract_busco_genes.benchmark.txt"

    shell:
        """(btk pipeline extract-busco-genes \
            --busco {params.busco} \
            --out {output.fasta}) 2> {log}"""
