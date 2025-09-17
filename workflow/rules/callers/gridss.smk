rule gridss:
    input:
        ref=config["reference"],
        ref_index=multiext(config["reference"], ".pac", ".fai"),
        bam=config["sample_path"] + "{sample}.bam",
    output:
        vcf=config["output"] + "/callers/gridss/{sample}.vcf",
    log:
        config["output"] + "/logs/gridss/{sample}.log",
    threads: 40
    resources:
        mem_mb=80000,
        walltime_h=10,
    benchmark:
        config["output"] + "/benchmarks/gridss/{sample}.tsv"
    shadow:
        "shallow"
    shell:
        """
        module load gridss/{config[gridss_version]}
        mkdir -p $(dirname {output.vcf})
        ( gridss -r {input.ref} -o {output.vcf} -t {threads} --jvmheap {resources.mem_mb}m {input.bam} ) > {log} 2>&1 

        """
