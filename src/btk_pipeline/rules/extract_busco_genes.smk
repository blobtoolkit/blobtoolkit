rule extract_busco_genes:
    """
    Extract busco genes into a single fasta file.
    """
    input:
        busco = expand("%s/%s.busco.{lineage}/full_table.tsv.gz" % (busco_path, config["assembly"]["prefix"]), lineage=get_basal_lineages(config)),
    output:
        fasta = "{assembly}.busco_genes.fasta"
    threads: 1
    log:
        "logs/{assembly}/extract_busco_genes.log"
    benchmark:
        "logs/{assembly}/extract_busco_genes.benchmark.txt"
    # script:
    #     "../lib/extract_busco_genes.py"
    shell:
        """(btk pipeline extract-busco-genes \
            --in {input.busco} \
            --out {output.fasta}) 2> {log}"""
