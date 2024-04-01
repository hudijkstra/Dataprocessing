
rule subset_vcf_by_chromosome:
    input:
        genome=config["reference_genome"] + config["reference_genome_ext"],
        vcf=config["results_dir"] + "combined_variants/all_samples_combined.vcf"
    output:
        expand(config["results_dir"] + "chromosome_subsets/{chr_id}_variants.vcf", chr_id=config["chromosome_name"])
    benchmark:
        "../logs/benchmarks/variant_analysis/subset_chromosome.txt"
    params:
        chromosome=config["chromosome_name"]
    log:
        "../logs/variant_analysis/subset_chromosome.log"
    shell:
       "{config[java]} -jar {config[gatk]} SelectVariants -R {input.genome} -V {input.vcf} -O {output} -L {params.chromosome} > {log}"

rule analyse_variants:
    input:
        vcf=expand( "{results_dir}chromosome_subsets/{chr_id}_variants.vcf", results_dir=config["results_dir"], chr_id=config["chromosome_name"])
    output:
        expand("{results_dir}variant_evaluation/{chr_id}.pdf",results_dir=config["results_dir"], chr_id=config["chromosome_name"])
    message:
        "VARIANT EVALUATION: Generating plots..."
    benchmark:
        "../logs/benchmarks/variant_analysis/chorm.txt"
    params:
        genome=config["reference_genome"] + config["reference_genome_ext"],
        gff=check_directory("../resources/", ".gff")
    resources:
        mem_mb=config["variant_analysis_params"]["memory_mb"]
    threads:
        config["variant_analysis_params"]["threads"]
    log:
        "../logs/variant_analysis/.log"
    shell:
        "Rscript scripts/variant_analysis.R {config[results_dir]} {input.vcf} {params.genome} {params.gff} > {log} 2>&1"