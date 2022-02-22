rule add_summary_to_metadata:
    output:
        "{assembly}.meta.yaml"
    threads: 1
    log:
        "logs/{assembly}/add_summary_to_metadata.log"
    benchmark:
        "logs/{assembly}/add_summary_to_metadata.benchmark.txt"
    script:
        "../lib/add_summary_to_metadata.py"
