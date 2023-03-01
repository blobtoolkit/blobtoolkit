rule unchunk_blastn:
    """
    Unchunk chunked blastn results.
    """
    input:
        "{assembly}.blastn.nt.out.raw"
    output:
        "{assembly}.blastn.nt.out"
    params:
        max_target_seqs = similarity_setting(config, "blastn", "max_target_seqs")
    threads: 1
    log:
        "logs/{assembly}/unchunk_blastn.log"
    benchmark:
        "logs/{assembly}/unchunk_blastn.benchmark.txt"
    shell:
        """(btk pipeline unchunk-blast \
            --in {input} \
            --count {params.max_target_seqs} \
            --out {output}) 2> {log}"""
