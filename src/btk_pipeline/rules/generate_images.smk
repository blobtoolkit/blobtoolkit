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
        host = "http://localhost",
        ports = "8000-8099",
        timeout = set_view_timeout(config)
    threads: 3
    log:
        "logs/{blobdir}/generate_images.log"
    benchmark:
        "logs/{blobdir}/generate_images.benchmark.txt"
    script:
        """../lib/generate_static_images.py"""
