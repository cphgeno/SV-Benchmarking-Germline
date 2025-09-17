rule delly:
    input:
        bam=config["sample_path"] + "{sample}.bam",
        ref_index=multiext(config["reference"], ".pac", ".fai"),
        ref=config["reference"],
    output:
        vcf=config["output"] + "/callers/delly/{sample}.vcf",
    resources:
        walltime_h=15,
        mem_mb=80000,
    benchmark:
        config["output"] + "/benchmarks/delly/{sample}.tsv"
    log:
        config["output"] + "/logs/delly/{sample}.log",
    threads: 40
    shadow:
        "shallow"
    singularity:
        "/PATH/TO/singularity/delly.sif"
    envmodules:
        "delly2/" + config["delly_version"],
    shell:
        """

        export OMP_NUM_THREADS={threads}

        delly call -g {input.ref} -o {output.vcf} {input.bam} > {log} 2>&1

        """
