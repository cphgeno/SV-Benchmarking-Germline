rule octopus:
    input:
        ref=config["reference"],
        ref_index=multiext(config["reference"], ".pac", ".fai"),
        bam=config["sample_path"] + "{sample}.bam",
    output:
        raw_vcf=config["output"] + "/callers/octopus/{sample}.vcf.gz",
    resources:
        walltime_h=200,
        mem_mb=80000,
    params:
        reference_dir=os.path.dirname(config["reference"]),
        output_dir=config["output"] + "/callers/octopus",
        temp_dir=config["output"] + "/temp/octopus/{sample}",
    threads: 40
    benchmark:
        config["output"] + "/benchmarks/octopus/{sample}.tsv"
    log:
        config["output"] + "/logs/octopus/{sample}.log",
    envmodules:
        "octopus/" + config["octopus_version"],
    shell:
        """
        echo "Running octopus on {wildcards.sample}" > {log}
        octopus --reference {input.ref} --reads {input.bam} --output {output.raw_vcf} --threads {threads}  --debug {log}.debug >> {log} 2>&1  
        """


# TODO move to general
rule bam_coverage:
    input:
        bam=config["sample_path"] + "{sample}.bam",
    output:
        coverage_bed=config["sample_path"] + "{sample}_coverage.bed",
        nonzero_bed=config["sample_path"] + "{sample}_nonzero.bed",
    log:
        profile=config["output"] + "/logs/{sample}-coverage.log",
    resources:
        walltime_h=24,
        mem_mb=80000,
    benchmark:
        config["output"] + "/benchmarks/coverage/{sample}.tsv"
    threads: 40
    envmodules:
        "bedtools/1.30.0",
    shell:
        """
        bedtools genomecov -ibam {input.bam} -bga > {output.bed}
        """
