rule mark_duplicates:
    input:
        genome=config["reference_genome"] + config["reference_genome_ext"],
        bam="../results/aligned/{sample}.bam"
    output:
        config["results_dir"] + "post_processing/marked_duplicates/{sample}.bam"
    message:
        "POST-PROCESSING: Marking duplicates for sample: {wildcards.sample}..."
    log:
        metrics="../logs/post_processing/marked_duplicates/{sample}_metrics.txt",
        stoud_log="../logs/post_processing/marked_duplicates/{sample}.log"
    shell:
        """
        {config[java]} -jar {config[gatk]} MarkDuplicates \
        --REFERENCE_SEQUENCE {input.genome} \
        --INPUT {input.bam} \
        --METRICS_FILE {log.metrics} \
        --OUTPUT {output} \
        --REMOVE_DUPLICATES true \
        2> {log.stoud_log}
        """

rule add_read_group:
    input:
        config["results_dir"] + "post_processing/marked_duplicates/{sample}.bam"
    output:
        bam=config["results_dir"] + "post_processing/grouped/{sample}.bam",
        bai=config["results_dir"] + "post_processing/grouped/{sample}.bai"
    message:
        "POST-PROCESSING: Adding read group for sample: {wildcards.sample}..."
    params:
        SORT_ORDER=config["read_group_params"]["SORT_ORDER"],
        RGID=config["read_group_params"]["RGID"],
        RGLB=config["read_group_params"]["RGLB"],
        RGPL=config["read_group_params"]["RGSM"],
        RGSM=config["read_group_params"]["RGSM"],
        RGPU=config["read_group_params"]["RGPU"]
    log:
        "../logs/post_processing/grouped/{sample}.log"
    shell:
        """
        {config[java]} -jar {config[gatk]} AddOrReplaceReadGroups \
        -I {input} \
        -O {output.bam} \
        -SORT_ORDER {params.SORT_ORDER} \
        -RGID {params.RGID} \
        -RGLB {params.RGLB} \
        -RGPL {params.RGPL} \
        -RGSM {params.RGSM} \
        -RGPU {params.RGPU} \
        -CREATE_INDEX True \
        2> {log}
        """
        