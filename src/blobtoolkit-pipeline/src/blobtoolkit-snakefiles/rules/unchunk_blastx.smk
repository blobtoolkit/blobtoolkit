rule unchunk_blastx:
    """
    Unchunk chunked blastx results.
    """
    input:
        "{assembly}.diamond.reference_proteomes.out.raw"
    output:
        "{assembly}.diamond.reference_proteomes.out"
    params:
        max_target_seqs = similarity_setting(config, "diamond_blastx", "max_target_seqs")
    threads: 1
    log:
        "logs/{assembly}/unchunk_blast.log"
    benchmark:
        "logs/{assembly}/unchunk_blast.benchmark.txt"
    shell:
        """(btk pipeline unchunk-blast \
            --in {input} \
            --count {params.max_target_seqs} \
            --out {output}) 2> {log}"""