rule run_diamond_blastx:
    """
    Run Diamond blastx to search protein database with assembly query.
    """
    input:
        fasta = "{assembly}.fasta.chunks",
        dmnd = "%s/%s.dmnd" % (similarity_setting(config, "diamond_blastx", "path"), similarity_setting(config, "diamond_blastx", "name"))
    output:
        "{assembly}.diamond.reference_proteomes.out.raw"
    params:
        evalue = similarity_setting(config, "diamond_blastx", "evalue"),
        max_target_seqs = similarity_setting(config, "diamond_blastx", "max_target_seqs"),
        taxid_flag = taxid_flag(config, "diamond_blastx")
    threads: 32
    log:
        "logs/{assembly}/run_diamond_blastx.log"
    benchmark:
        "logs/{assembly}/run_diamond_blastx.benchmark.txt"
    shell:
        """(if [ -s {input.fasta} ]; then \
            diamond blastx \
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