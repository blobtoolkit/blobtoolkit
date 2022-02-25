rule checksum_files:
    """
    Calculate SHA1 checksum for all files in dataset.
    """
    input:
        "{blobdir}/summary.json"
    output:
        "{blobdir}/CHECKSUM"
    params:
        blobdir = lambda wc: wc.blobdir
    threads: 1
    log:
        "logs/{blobdir}/checksum_files.log"
    benchmark:
        "logs/{blobdir}/checksum_files.benchmark.txt"
    shell:
        """(find {params.blobdir}/ -type f -exec sha1sum {{}} ';' \
        | sort -k 2 \
        | sed 's:{params.blobdir}/::' > {params.blobdir}/CHECKSUM) 2> {log}"""
