rule tiddit:
    input:
        ref=config["reference"],
        ref_index=multiext(config["reference"], ".pac", ".fai"),
        bam=config["sample_path"] + "{sample}.bam",
    output:
        vcf=config["output"] + "/callers/tiddit/{sample}.vcf",
        tmp=temp(directory(config["output"] + "/temp/tiddit/{sample}")),
    params:
        outdir=config["output"] + "/callers/tiddit/",
    benchmark:
        config["output"] + "/benchmarks/tiddit/{sample}.tsv"
    log:
        config["output"] + "/logs/tiddit/{sample}.log",
    threads: 40
    resources:
        walltime_h=24,
        mem_mb=80000,
    singularity:
        "/PATH/TO/singularity/tiddit.sif"
    envmodules:
        "tiddit/" + config["tiddit_version"],
    shell:
        """
        (rm -rf {output.tmp}*
        mkdir -p {params.outdir}
        echo 'Starting tiddit on {wildcards.sample}'
        tiddit --sv --bam {input.bam} -o {output.tmp} --ref {input.ref} --threads {threads}
        echo 'Moving tiddit output from {output.tmp} to {params.outdir}'
        cp --verbose --remove-destination {output.tmp}* {params.outdir}) > {log} 2>&1

        """
