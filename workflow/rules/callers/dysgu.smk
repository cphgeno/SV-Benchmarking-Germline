rule dysgu:
    input:
        ref=config["reference"],
        ref_index=multiext(config["reference"], ".pac", ".fai"),
        bam=config["sample_path"] + "{sample}.bam",
    output:
        vcf=config["output"] + "/callers/dysgu/{sample}.vcf",
        tmp=temp(directory(config["output"] + "/temp/dysgu/{sample}")),
    benchmark:
        config["output"] + "/benchmarks/dysgu/{sample}.tsv"
    resources:
        walltime_h=24,
        mem_mb=80000,
    log:
        config["output"] + "/logs/dysgu/{sample}.log",
    params:
        reference_dir=os.path.abspath(os.path.dirname(config["reference"])),
        sample_dir=os.path.abspath(config["sample_path"]),
        output_dir=os.path.abspath(config["output"]),
        temp_dir=os.path.abspath(config["output"] + "/temp/dysgu/{sample}"),
    threads: 40
    shell:
        """
        ml dysgu/1.6.2
        mkdir -p {params.temp_dir} > {log} 2>&1 
        export CONTAINER_WORK_DIRS={params.reference_dir},{params.output_dir},{params.sample_dir},{params.temp_dir}
        dysgu run -v 1 -x -p2 -o {output.vcf} {input.ref} {output.tmp} {input.bam} > {log} 2>&1
        """
