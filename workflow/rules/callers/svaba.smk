rule svaba:
    input:
        ref=config["reference"],
        ref_index=multiext(config["reference"], ".pac", ".fai"),
        bam=config["sample_path"] + "{sample}.bam",
    output:
        sv=config["output"] + "/callers/svaba/{sample}.svaba.sv.vcf",
        indel=config["output"] + "/callers/svaba/{sample}.svaba.indel.vcf",
    params:
        dir=config["output"] + "/callers/svaba/",
        sample="{sample}",
    threads: 40
    resources:
        walltime_h=48,
        mem_mb=80000,
    benchmark:
        config["output"] + "/benchmarks/svaba/{sample}.tsv"
    log:
        config["output"] + "/logs/svaba/{sample}.log",
    envmodules:
        "svaba/" + config["svaba_version"],
    shell:
        """

        (svaba run -t {input.bam} -p {threads} --germline -a {wildcards.sample} -G {input.ref};

        mkdir -p {params.dir};

        mv {wildcards.sample}* {params.dir};) > {log} 2>&1  

        """


rule svaba_combine:
    input:
        sv=rules.svaba.output.sv,
        indel=rules.svaba.output.indel,
    output:
        vcf=config["output"] + "/callers/svaba/{sample}.vcf.gz",
    threads: 1
    resources:
        walltime_h=1,
        mem_mb=10000,
    benchmark:
        config["output"] + "/benchmarks/svaba-postprocess/{sample}.tsv"
    log:
        config["output"] + "/logs/svaba/postprocess-{sample}.log",
    envmodules:
        "bcftools/" + config["bcftools_version"],
    shell:
        """

        (bcftools view {input.sv} -o {input.sv}.gz -Oz;

        bcftools index -t -f {input.sv}.gz;

        bcftools view {input.indel} -o {input.indel}.gz -Oz;

        bcftools index -t -f {input.indel}.gz;

        bcftools concat -a {input.sv}.gz {input.indel}.gz | bcftools annotate -c 'INFO/SVLEN:=INFO/SPAN' -o {output.vcf} -Oz;

        ) > {log} 2>&1

        """
