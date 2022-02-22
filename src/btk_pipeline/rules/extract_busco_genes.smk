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
    script:
        "../lib/extract_busco_genes.py"
    # shell:
    #     """(> {output}; \
    #     for TABLE in {input.busco}; do \
    #         if [ -s $TABLE ]; then \
    #             SEQS=${{TABLE/full_table.tsv.gz/busco_sequences.tar.gz}};
    #             tar xf $SEQS \
    #                 --to-command='FILE=$(basename $TAR_FILENAME); \
    #                               TYPE=$(echo $TAR_FILENAME | awk -F '"'"'[/_]'"'"' '"'"'{{print $3}}'"'"');
    #                               awk -v busco=${{FILE%.faa}} -v type=$TYPE '"'"'{{if($1 ~ /^>/){{print $1 "=" busco "=" type}} else {{print $1}}}}'"'"''; \
    #         fi; \
    #     done) >> {output} 2> {log}"""
