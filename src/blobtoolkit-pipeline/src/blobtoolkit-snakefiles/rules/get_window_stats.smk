rule get_window_stats:
    """
    Get average sequence stats across windows.
    """
    input:
        tsv = "%s/{assembly}.chunk_stats.tsv" % cov_stats_path
    output:
        tsv = "{assembly}.window_stats.tsv"
    params:
        window = " --window ".join([str(window) for window in set_stats_windows(config)]),
    threads: 1
    log:
        "logs/{assembly}/get_window_stats.log"
    benchmark:
        "logs/{assembly}/get_window_stats.benchmark.txt"
    shell:
        """(btk pipeline window-stats \
            --in {input.tsv} \
            --window {params.window} \
            --out {output.tsv}) 2> {log}"""