rule smoove:
    input:
        ref=config["reference"],
        ref_index=multiext(config["reference"], ".pac", ".fai"),
        cram=config["sample_path"] + "{sample}.cram",
    output:
        vcf=config["output"] + "/callers/smoove/{sample}.vcf.gz",
    params:
        sample="{sample}",
    shadow:
        "shallow"
    threads: 40
    resources:
        walltime_h=24,
        mem_mb=80000,
    benchmark:
        config["output"] + "/benchmarks/smoove/{sample}.tsv"
    log:
        config["output"] + "/logs/smoove/{sample}.log",
    singularity:
        "/PATH/TO/singularity/smoove.sif"
    envmodules:
        "anaconda2/" + config["anaconda2_version"],
        "samblaster/" + config["samblaster_version"],
        "sambamba/" + config["sambamba_version"],
        "bcftools/" + config["bcftools_version"],
        "samtools/" + config["samtools_version"],
        "gsort/" + config["gsort_version"],
        "lumpy/" + config["lumpy_version"],
        "smoove/" + config["smoove_version"],
    shell:
        """
        smoove call --name tmpdir_{params.sample} --fasta {input.ref} --genotype -x {input.cram} -p {threads} > {log} 2>&1
        mv tmpdir_{params.sample}-smoove.genotyped.vcf.gz {output.vcf}
        """


