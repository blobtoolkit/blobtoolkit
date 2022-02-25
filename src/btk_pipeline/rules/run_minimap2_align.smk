rule run_minimap2_align:
    """
    Run minimap2 read alignment.
    """
    input:
        fasta = config["assembly"]["file"],
        index = lambda wc: "%s.%s.mmi" % (wc.assembly, minimap_tuning(config, wc.sra)),
        fastq = lambda wc: read_files(config, wc.sra)
    output:
        bam = "{assembly}.{sra}.bam",
        csi = "{assembly}.{sra}.bam.csi"
    params:
        tuning = lambda wc: minimap_tuning(config, wc.sra),
        assembly = lambda wc: wc.assembly,
        subsample = lambda wc: seqtk_sample_input(config, wc.sra)
    threads: 24
    log:
        "logs/{assembly}/run_minimap2_align/{sra}.log"
    benchmark:
        "logs/{assembly}/run_minimap2_align/{sra}.benchmark.txt"
    shell:
        """(minimap2 -t {threads} -ax {params.tuning} {input.index} {params.subsample} | \
        samtools view -h -T {input.fasta} - | \
        samtools sort -@4 -O BAM -o {output.bam} - &&
        samtools index -c {output.bam} {output.csi}) &> {log}"""