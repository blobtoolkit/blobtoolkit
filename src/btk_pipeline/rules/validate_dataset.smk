rule validate_dataset:
    """
    Run BlobToolKit validator on a dataset to check all expected fields are present.
    """
    input:
        copied = "{blobdir}.copied",
        cov = expand("{{blobdir}}/{sra}_cov.json", sra=reads_by_prefix(config).keys()),
        tax = "{blobdir}/buscogenes_phylum_positions.json",
        busco = "{blobdir}/%s_busco.json" % config["busco"]["lineages"][0],
        ids = "{blobdir}/identifiers.json"
    output:
        touch(temp("{blobdir}.valid"))
    params:
        blobdir = lambda wc: wc.blobdir
    threads: 1
    log:
        "logs/{blobdir}/validate_dataset.log"
    benchmark:
        "logs/{blobdir}/validate_dataset.benchmark.txt"
    shell:
        """blobtools validate {params.blobdir} > {log} 2>&1"""
