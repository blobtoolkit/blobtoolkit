rule add_summary_to_metadata:
    output:
        "{assembly}.meta.yaml"
    params:
        parent_dir = parent_dir,
    threads: 1
    log:
        "logs/{assembly}/add_summary_to_metadata.log"
    benchmark:
        "logs/{assembly}/add_summary_to_metadata.benchmark.txt"
    shell:
        """(btk pipeline add-summary-to-metadata \
            --config {params.parent_dir}/config.yaml \
            --out {output}) 2> {log}"""