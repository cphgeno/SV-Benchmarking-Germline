rule summary_analysis:
    input:
        expand(
            rules.truvari.output.summary,
            type=config["truvari_runs"],
            sample=SAMPLES,
            caller=CALLERS,
            output=config["output"],
        ),
        script="workflow/r-scripts/01-truvari-summary.R",
    output:
        config["output"] + "/truvari/stat_summary.png",
    log:
        config["output"] + "/logs/truvari-summary.log",
    resources:
        mem_mb=8000,
        walltime_h=1,
    envmodules:
        "R/" + config["R_version"],
        "R-bundle-Bioconductor/3.20-foss-2024a-R-4.4.2",
    shell:
        """
        Rscript {input.script} > {log} 2>&1
        """


rule vcf_analysis:
    input:
        expand(
            rules.truvari.output.summary,
            type=config["truvari_runs"],
            sample=SAMPLES,
            caller=CALLERS,
            output=config["output"],
        ),
        script="workflow/r-scripts/02-vcf-investigation.R",
        data=rules.summary_analysis.output,
    output:
        config["output"] + "/truvari/type-ignored/pctsizesimilarity_violin.png",
    log:
        config["output"] + "/logs/truvari-summary.log",
    resources:
        mem_mb=16000,
        walltime_h=2,
    envmodules:
        "openssl/" + config["openssl_version"],
        "R/" + config["R_version"],
    shell:
        """

        Rscript --vanilla {input.script} > {log} 2>&1

        """


rule upset_analysis:
    input:
        json=expand(
            rules.truvari_consistency.output.json,
            type="type-ignored",
            sample="NA12878",
            svtype=SVTYPES,
            output=config["output"],
        ),
        script="workflow/r-scripts/upsetr_const.R",
        data=rules.summary_analysis.output,
    output:
        config["output"] + "/truvari/type-ignored/_upset-degreesort.png",
    log:
        config["output"] + "/logs/truvari-upset.log",
    resources:
        mem_mb=8000,
        walltime_h=1,
    envmodules:
        "openssl/" + config["openssl_version"],
        "R/" + config["R_version"],
    shell:
        """

        Rscript --vanilla {input.script} > {log} 2>&1

        """


rule width_analysis:
    input:
        expand(
            rules.truvari.output.summary,
            type=config["truvari_runs"],
            sample=SAMPLES,
            caller=CALLERS,
            output=config["output"],
        ),
        script="workflow/r-scripts/05-tpbase.R",
        data=rules.summary_analysis.output,
    output:
        config["output"] + "/truvari/type-ignored/count_plot.png",
    log:
        config["output"] + "/logs/truvari-upset.log",
    resources:
        mem_mb=8000,
        walltime_h=1,
    envmodules:
        "openssl/" + config["openssl_version"],
        "R/" + config["R_version"],
    shell:
        """

        Rscript --vanilla {input.script} FALSE > {log} 2>&1

        """


rule benchmark_analysis:
    input:
        expand(
            config["output"] + "/benchmarks/{caller}/{sample}.tsv",
            sample=SAMPLES,
            caller=CALLERS,
            output=config["output"],
        ),
        script="workflow/r-scripts/03-benchmarks.R",
    output:
        config["output"] + "/benchmarks/caller-summary.csv",
    log:
        config["output"] + "/logs/benchmark-summary.log",
    resources:
        mem_mb=8000,
        walltime_h=1,
    envmodules:
        "openssl/" + config["openssl_version"],
        "R/" + config["R_version"],
    shell:
        """

        R -f {input.script} > {log} 2>&1

        """
