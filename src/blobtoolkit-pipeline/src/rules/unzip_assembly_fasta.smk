rule unzip_assembly_fasta:
    """
    Unzip an assembly fasta file.
    """
    input:
        config["assembly"]["file"]
    output:
        temp("{assembly}.fasta")
    threads: 1
    log:
        "logs/{assembly}/unzip_assembly_fasta.log"
    benchmark:
        "logs/{assembly}/unzip_assembly_fasta.benchmark.txt"
    shell:
        "(gunzip -c {input} > {output}) &> {log}"
