rule get_chunked_stats:
    """
    Get chunked sequence stats.
    """
    input:
        fasta = "%s/{assembly}.windowmasker.fasta" % windowmasker_path,
    output:
        bed = "{assembly}.chunk_stats.mask.bed",
        tsv = "{assembly}.chunk_stats.tsv",
    params:
        chunk = set_stats_chunk(config),
        overlap = 0,
        max_chunks = 1000000000,
        min_length = set_stats_chunk(config),
        bed = lambda wc: "%s.chunk_stats.bed" % wc.assembly
    threads: 1
    log:
        "logs/{assembly}/get_seq_stats.log"
    benchmark:
        "logs/{assembly}/get_seq_stats.benchmark.txt"
    shell:
        """(btk pipeline chunk-fasta \
            --in {input.fasta} \
            --chunk {params.chunk} \
            --overlap {params.overlap} \
            --max-chunks {params.max_chunks} \
            --min-length {params.min_length} \
            --busco None \
            --out None \
            --bed {params.bed}) 2> {log}"""