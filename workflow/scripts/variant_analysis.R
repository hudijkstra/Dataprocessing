# Variant Evaluation Script

# Purpose:
# The Variant Evaluation module assesses and visualizes genetic variants stored in a VCF file. 
# It generates a quality control plot for a specified chromosome, considering read depth (DP) information. 
# It supports the inclusion of genomic annotations from a GFF file.

# Usage:
# Run from the command line with user-provided arguments:
#   --input-vcf: Input VCF file
#   --reference-seq: Reference sequence in FASTA format
#   --output-dir: Output directory to save the plot
#   --annotation-file (optional): GFF file providing genomic annotations

# Input:
# 1. VCF file: Contains genetic variant information.
# 2. Reference sequence (FASTA): Represents the genomic reference for alignment.
# 3. Output directory: Specifies the location to save the generated quality control plot.
# 4. Annotation file (optional): GFF file providing genomic annotations for visualizing features on the plot.

# Output:
# A quality control plot in PDF format, highlighting read depth information and genomic annotations if provided. 
# The plot is saved in the specified output directory.

# Example Usage:
# Rscript variant_evaluation.R --input-vcf input.vcf --reference-seq reference.fasta --output-dir output/ --annotation-file annotations.gff

# Load necessary libraries
library("vcfR")
library("ape")
library("stringr")

# Get command line arguments
args <- commandArgs()

# Extract output directory from command line arguments
output_dir <- args[6]

# Extract chromosome name
chromosome <- str_match(args[7], "/chromosome_subsets/([^/]+)_variants\\.vcf")[2]

# Read VCF file
vcf <- read.vcfR(args[7], verbose = TRUE)

# Read reference sequence in FASTA format
ref_seq <- ape::read.dna(args[8], format = "fasta")
ref_seq_subset <- ref_seq[grep(chromosome, names(ref_seq))]
names(ref_seq_subset) <- chromosome
ref_seq_subset <- as.matrix(ref_seq_subset)

# Construct the PDF path
pdf_path <- paste0(output_dir, "variant_evaluation/", chromosome, ".pdf")

# Start PDF device
pdf(pdf_path)

# Check if annotation file is provided
if (args[9] == "None") {
    # If no annotation, create chromosome without annotation
    chrom <- create.chromR(name=chromosome, vcf=vcf)
} else {
    # If annotation is provided, read GFF file and create chromosome with annotation
    gff <- read.table(args[9], sep="\t", quote="")
    gff_subset <- gff[grep(chromosome, gff[,1]),]
    chrom <- create.chromR(name=chromosome, vcf=vcf, seq=ref_seq_subset, ann=gff_subset)
}

# Process the chromosome
chrom <- proc.chromR(chrom, verbose=TRUE)

# Generate chromosome quality control plot
chromoqc(chrom, dp.alpha=20)

# Close the PDF device
dev.off()