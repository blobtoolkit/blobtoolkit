rule run_diamond_blastp:
    """
    Run Diamond blastp to search protein database with busco gene query.
    """
    input:
        fasta = "{assembly}.busco_genes.fasta",
        dmnd = "%s/%s.dmnd" % (similarity_setting(config, "diamond_blastp", "path"), similarity_setting(config, "diamond_blastx", "name"))
    output:
        "{assembly}.diamond.busco_genes.out"
    params:
        evalue = similarity_setting(config, "diamond_blastp", "evalue"),
        max_target_seqs = similarity_setting(config, "diamond_blastp", "max_target_seqs"),
        taxid_flag = taxid_flag(config, "diamond_blastp")
    threads: 32
    log:
        "logs/{assembly}/run_diamond_blastp.log"
    benchmark:
        "logs/{assembly}/run_diamond_blastp.benchmark.txt"
    shell:
        """(if [ -s {input.fasta} ]; then \
            diamond blastp \
                --query {input.fasta} \
                --db {input.dmnd} \
                --outfmt 6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
                --max-target-seqs {params.max_target_seqs} \
                --max-hsps 1 \
                --evalue {params.evalue} \
                --threads {threads} {params.taxid_flag} \
                > {output}; \
        else \
            > {output}; \
        fi) 2> {log}"""