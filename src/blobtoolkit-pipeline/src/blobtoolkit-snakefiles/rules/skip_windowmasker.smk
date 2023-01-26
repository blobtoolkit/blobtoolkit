rule skip_windowmasker:
    """
    Skip windowmasker step.
    """
    input:
        fasta = "%s.fasta" % config["assembly"]["prefix"]
    output:
        masked = "%s.windowmasker.fasta" % config["assembly"]["prefix"]
    threads: 1
    log:
        "logs/%s/skip_windowmasker.log" % config["assembly"]["prefix"]
    benchmark:
        "logs/%s/skip_windowmasker.benchmark.txt" % config["assembly"]["prefix"]
    shell:
        """cp {input.fasta} {output.masked} 2>> {log}"""
