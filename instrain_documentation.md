
# Strain Tracking in Gut Metagenomes of Mother-Infants Using InStrain and Complementary Tools

## Description
This script is designed for strain-level analysis in gut metagenomes of mother-infant studies. It integrates InStrain with other bioinformatics tools for insights into strain diversity and transmission. The workflow includes mapping reads, predicting genes, and strain profiling.

## Requirements
- **InStrain**: Microbial strain identification and comparison.
- **bowtie2**: Aligning sequencing reads.
- **samtools**: SAM format manipulation.
- **prodigal**: Gene prediction.
- **GNU parallel**: Parallel processing.

## Usage

### 1. Preparing Input Data
Ensure sequencing data and reference genomes are correctly formatted.

### 2. Generating Scaffold to Bin File
Create a mapping between scaffolds and bins using drep.

### 3. Gene Prediction with Prodigal
Predict genes in metagenomic data using Prodigal.

### 4. Creating Bowtie Database
Build a bowtie2 database from reference genomes.

### 5. Mapping Reads
Map sequencing reads to the database with bowtie2.

### 6. Running InStrain
Execute InStrain for strain-level profiling and comparison.

### 7. Data Analysis and Interpretation
Analyze InStrain output focusing on strain diversity and mother-infant sharing patterns.
