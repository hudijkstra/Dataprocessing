rule index_reference_genome:
    input:
        config["reference_genome"] + config["reference_genome_ext"]
    output:
        expand("{reference_genome}{reference_genome_ext}.{ext}", 
                reference_genome=config["reference_genome"], 
                reference_genome_ext=config["reference_genome_ext"],
                ext=["amb", "ann", "bwt", "pac", "sa"])
    message:
        "ALIGNMENT: Indexing reference genome.."
    log:
        "../logs/alignment/bwa_index.log",
    shell:
        "bwa index {input} 2> {log}"

rule samtools_faidx:
    input:
        config["reference_genome"] + config["reference_genome_ext"]
    output:
        "{genome}.fai"
    message:
        "ALIGNMENT: Generating .fai file..."
    shell:
        "samtools faidx {input}"

rule create_sequence_dict:
    input:
        config["reference_genome"] + config["reference_genome_ext"]
    output:
        config["reference_genome"] + ".dict"
    message:
        "ALIGNMENT: Generating .dict file..."
    log:
        "../logs/alignment/sequence_dict.log",
    shell:
        "{config[java]} -jar {config[gatk]} CreateSequenceDictionary -R {input} 2> {log}"