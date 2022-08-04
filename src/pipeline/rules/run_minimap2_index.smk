rule run_minimap2_index:
    """
    Run minimap2 reference genome indexing.
    """
    input:
        config["assembly"]["file"]
    output:
        "{assembly}.{tuning}.mmi"
    params:
        tuning = lambda wc: wc.tuning,
        assembly = lambda wc: wc.assembly,
        span = config["assembly"].get("span", 0)
    threads: 3
    log:
        "logs/{assembly}/run_minimap2_index/{tuning}.log"
    benchmark:
        "logs/{assembly}/run_minimap2_index/{tuning}.benchmark.txt"
    shell:
        """(SPAN={params.span}; \
        if [ $SPAN == 0 ]; then \
            SPAN=$(gunzip -c {input} | grep -v ">" | tr -d '\\n' | wc -c); \
        fi; \
        minimap2 -t {threads} -x {params.tuning} -I $SPAN -d {output} {input}) &> {log}"""