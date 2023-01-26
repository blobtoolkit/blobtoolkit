rule run_windowmasker:
    """
    Run windowmasker to mask repeats in the assembly.
    """
    input:
        fasta = "%s.fasta" % config["assembly"]["prefix"]
    output:
        counts = "%s.windowmasker.counts" % config["assembly"]["prefix"],
        masked = "%s.windowmasker.fasta" % config["assembly"]["prefix"]
    threads: 1
    log:
        "logs/%s/run_windowmasker.log" % config["assembly"]["prefix"]
    benchmark:
        "logs/%s/run_windowmasker.benchmark.txt" % config["assembly"]["prefix"]
    shell:
        """windowmasker -in {input.fasta} \
                     -infmt fasta \
                     -mk_counts \
                     -sformat obinary \
                     -out {output.counts} 2> {log} && \
        windowmasker -in {input} \
                     -infmt fasta \
                     -ustat {output.counts} \
                     -dust T \
                     -outfmt fasta \
                     -out {output.masked} 2>> {log}"""
