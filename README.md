**Thema 11 - Dataprocessing**  
Author: H. Dijkstra  
Date: 30th March 2024  
Version: 1.0

This repo contains the snakemake pipeline made as a learning project for the minor "High Throughput / High-performance Biocomputing". The pipeline is inspired by this [article](https://bmcbioinformatics.biomedcentral.com/counter/pdf/10.1186/s12859-016-1431-9.pdf). 

---
## Prerequisites  
The program is dependant on the following tools and packages:  

- GATK     (version 4.5.0.0)
- cutadapt (version 4.6)
- samtools (version >= 1.13)
- BWA      (version 0.7.17)
- vcfR     (version 1.15.0)
- Java     (jdk/jre >= 17)

## Installation  
Clone the repository to your device.  
- Download [Miniconda or Anaconda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/download.html) latest version.
- Create a conda environment from the yml file.  
`conda env create -f conda_env.yml` 
- Install samtools.
`sudo apt install samtools=1.13`
- Install BWA.
`sudo apt install bwa=0.7.17 `
- Download [GATK](http://hgdownload.cse.ucsc.edu/admin/exe/) version 4.5.0.0.  
- Add the path to `gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar` in the `paths.yaml` file.
- Add the path to the java installation to the `paths.yaml` file.

## Data
- The pipeline was tested on data available at NCBI sequence read archive (SRA)# SRP059747 with reference genome [Genome assembly Glycine_max_v4.0](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000004515.6/)

## Set-up
1. Place your sample files in the `data` directory and **configure the extension** in the `paths.yaml` file. The program will automatically detect the sample files.
2. Place your reference genome file in the `resources` directory and once again **configure the extension** in the `paths.yaml` file.
3. (Optional) Place the .gff annotation file in the `resources` directory. Make sure the name matches the reference genome.
4. Add the chromosome name to the `parameters.yaml` file.

The paths.yaml and parameters.yaml files include supplementary settings for user customization.

## Usage
1. Make sure you have created the conda environment (Look at the installation caption.)
2. Open the terminal.
3. Activate the conda environment.  
`conda activate conda_env` 
4. Navigate to the `Snakefile` located in the `workflow` folder. 
5. Run the following command:

```bash
snakemake -j [number of jobs]
```

## Output
The pipeline produces and processes a VCF file, generating a PDF report located in results/variant_evaluation. This report includes chromosome-wise annotation density and read depth/quality metrics, aiding in the analysis of genetic variations.

## Rules

`Rule: trim_reads`
This rule trims reads using Cutadapt based on specified quality cutoff, saving trimmed reads to pre-processed directory in fastq.gz format, with associated logging and benchmarking.

`Rule: align_reads`
This rule aligns reads to the reference genome using BWA mem, producing a sorted BAM file as output, with associated logging and benchmarking.

`Rule: index_reference_genome`
This rule indexes the reference genome using BWA, generating auxiliary files required for alignment.

`Rule: samtools_faidx`
This rule generates a .fai index file for the reference genome using samtools.

`Rule: create_sequence_dict`
This rule creates a .dict file for the reference genome using GATK's CreateSequenceDictionary tool.

`Rule: mark_duplicates`
This rule marks duplicate reads in the aligned BAM file using GATK's MarkDuplicates tool, generating metrics and logging.

`Rule: add_read_group`
This rule adds read group information to the BAM file, essential for downstream analysis, using GATK's AddOrReplaceReadGroups tool, with associated logging.

`Rule: call_variants`
This rule performs variant calling on the grouped BAM file using GATK's HaplotypeCaller, producing a VCF file, with associated benchmarking and logging.

`Rule: create_vcf_list`
This rule creates a list of VCF files for all samples to be merged.

`Rule: merge_vcfs`
This rule merges individual VCF files from multiple samples into a single combined VCF file using GATK's MergeVcfs tool.

`Rule: subset_vcf_by_chromosome`
This rule subsets the combined VCF file by chromosome using GATK's SelectVariants tool, producing chromosome-specific VCF files for further analysis, with associated benchmarking and logging.

`Rule: analyse_variants`
This rule generates plots for variant evaluation using a custom R script, taking chromosome-specific VCF files as input along with the reference genome and GFF annotation file, with associated benchmarking and logging.

## Support 
- e-mail: h.u.dijkstra@st.hanze.nl
