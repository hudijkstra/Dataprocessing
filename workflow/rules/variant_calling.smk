rule call_variants:
    input:
        genome=config["reference_genome"] + config["reference_genome_ext"],
        fai=config["reference_genome"] + config["reference_genome_ext"] + ".fai",
        dict=config["reference_genome"] + ".dict",
        bam=config["results_dir"] + "post_processing/grouped/{sample}.bam",
        bai=config["results_dir"] + "post_processing/grouped/{sample}.bai"
    output:
        config["results_dir"] + "variants/{sample}.vcf"
    message:
        "VARIANT CALLING: calling variants for: {wildcards.sample}..."
    benchmark:
        "../logs/benchmarks/variant_calling/{sample}.txt"
    resources:
        mem_mb=config["variant_calling_params"]["memory_mb"]
    threads:
        config["variant_calling_params"]["threads"]
    log:
        "../logs/variant_calling/{sample}.log"
    shell:
        "{config[java]} -jar {config[gatk]} HaplotypeCaller -R {input.genome} -I {input.bam} -O {output} 2> {log}"

rule create_vcf_list:
    input:
        vcf_files=expand("{results_dir}variants/{sample}.vcf", sample=SAMPLES, results_dir=config["results_dir"])
    output:
        temp("vcf_files.list")
    shell:
        "echo {input.vcf_files} | tr ' ' '\\n' > {output}"
        
rule merge_vcfs:
    input:
        vcf_files_list="vcf_files.list"
    output:
        config["results_dir"] + "combined_variants/all_samples_combined.vcf"
    message:
        "VARIANT CALLING: Combining variants for all samples..."
    shell:
        "{config[java]} -jar {config[gatk]} MergeVcfs -I {input.vcf_files_list} -O {output}"