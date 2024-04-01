rule align_reads:
    input:
        sample="../results/pre_processed/{sample}.fastq.gz",
        genome=config["reference_genome"] + config["reference_genome_ext"],
        reference_files = expand(
            "{reference_genome}{reference_genome_ext}.{ext}",
            reference_genome=config["reference_genome"],
            reference_genome_ext=config["reference_genome_ext"],
            ext=["amb", "ann", "bwt", "pac", "sa"])
    output:
        config["results_dir"] + "aligned/{sample}.bam"
    message:
        "ALIGNMENT: Aligning reads for sample: {wildcards.sample}"
    benchmark:
        "../logs/benchmarks/{sample}.bwa.benchmark.txt"
    resources:
        mem_mb=config["alignment_params"]["memory_mb"]
    threads:
        config["alignment_params"]["threads"]
    log:
        "../logs/alignment/{sample}.log"
    shell:
        "bwa mem -t 4 {input.genome} {input.sample} | samtools sort -o {output}"