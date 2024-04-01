rule trim_reads:
    input:
        config["sample_dir"] + "{sample}" + config["sample_ext"]
    output:
        config["results_dir"] + "pre_processed/{sample}.fastq.gz"
    message:
        "PRE-PROCESSING: Trimming reads of sample: {wildcards.sample}..."
    params:
        quality_cutoff=config["trimming_parameters"]["quality_cutoff"]
    benchmark:
        "../logs/benchmarks/{sample}.trimming.benchmark.txt"
    threads:
        config["trimming_parameters"]["threads"]
    log:
        "../logs/pre_processing/{sample}.log"
    shell:
        "cutadapt -q {params.quality_cutoff} -o {output} {input} > {log} 2>&1"