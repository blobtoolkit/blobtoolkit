rule extract_nohit_fasta:
    """
    Extract sequences with no blastx hits into a separate file.
    """
    input:
        blastx = "%s/%s.diamond.reference_proteomes.out" % (diamond_path, config["assembly"]["prefix"]),
        fasta = "%s/{assembly}.windowmasker.fasta" % windowmasker_path
    output:
        "{assembly}.nohit.fasta"
    params:
        evalue = similarity_setting(config, "diamond_blastx", "import_evalue")
    threads: 4
    log:
        "logs/{assembly}/extract_nohit_fasta.log"
    benchmark:
        "logs/{assembly}/extract_nohit_fasta.benchmark.txt"
    shell:
        """seqtk subseq {input.fasta} <(grep '>' {input.fasta} | \
            grep -v -w -f <(awk '{{if($14<{params.evalue}){{print $1}}}}' {input.blastx} | \
                sort | uniq) \
            | cut -f1 | sed 's/>//') > {output}"""
