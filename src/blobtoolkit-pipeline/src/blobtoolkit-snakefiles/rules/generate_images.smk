rule generate_images:
    """
    Use BlobTools2 view to generate a set of static images.
    """
    input:
        valid = "{blobdir}.valid",
        cov = expand("{{blobdir}}/{sra}_cov.json", sra=reads_by_prefix(config).keys())
    output:
        "{blobdir}/cumulative.png"
    params:
        blobdir = lambda wc: wc.blobdir,
        cov = lambda wc: " --cov " if reads_by_prefix(config) else "",
        host = "http://localhost",
        ports = "8000-8099",
        timeout = set_view_timeout(config)
    threads: 3
    log:
        "logs/{blobdir}/generate_images.log"
    benchmark:
        "logs/{blobdir}/generate_images.benchmark.txt"
    shell:
        """(btk pipeline generate-static-images \
            --blobdir {params.blobdir} {params.cov} \
            --host {params.host} \
            --ports {params.ports} \
            --timeout {params.timeout}) 2> {log}"""
