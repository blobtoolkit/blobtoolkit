rule transfer_dataset:
    """
    Transfer dataset out of working directory.
    """
    input:
        "{assembly}/CHECKSUM"
    output:
        temp("{assembly}.complete"),
    params:
        assembly = lambda wc: wc.assembly,
        destdir = destination_dir,
    threads: 1
    log:
        "logs/{assembly}/transfer_dataset.log"
    benchmark:
        "logs/{assembly}/transfer_dataset.benchmark.txt"
    shell:
        """(rsync -av --remove-source-files {params.assembly}* {params.destdir}/ \
        && rmdir {params.assembly} \
        && touch {params.assembly}.complete) 2> {log}"""
