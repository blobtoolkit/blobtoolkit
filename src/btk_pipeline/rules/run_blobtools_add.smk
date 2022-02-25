rule run_blobtools_add:
    """
    Run blobtools add.
    """
    input:
        meta = "%s/meta.json" % blobdir_name(config),
        blastx = "%s/%s.diamond.reference_proteomes.out" % (diamond_path, config["assembly"]["prefix"]),
        blastn = "%s/%s.blastn.nt.out" % (blastn_path, config["assembly"]["prefix"]),
        taxdump = config["settings"]["taxdump"],
    output:
        "%s/buscoregions_phylum.json" % blobdir_name(config),
    params:
        blobdir = blobdir_name(config),
        evalue = similarity_setting(config, "diamond_blastx", "import_evalue"),
        max_target_seqs = similarity_setting(config, "diamond_blastx", "import_max_target_seqs"),
        update_plot = set_update_taxrule(config),
        fields = set_fields(config)
    threads: 4
    log:
        "logs/%s/run_blobtools_add.log" % config["assembly"]["prefix"]
    benchmark:
        "logs/%s/run_blobtools_add.benchmark.txt" % config["assembly"]["prefix"]
    shell:
        """blobtools replace \
            --taxdump {input.taxdump} \
            --taxrule bestdistorder=buscoregions {params.update_plot} \
            --hits {input.blastx} \
            --hits {input.blastn} \
            --evalue {params.evalue} \
            --hit-count {params.max_target_seqs} \
            --threads {threads} {params.fields} \
            {params.blobdir} > {log} 2>&1"""
