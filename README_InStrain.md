# InStrain Script for Gut Metagenome Analysis

## Description
This script is used for strain-level analysis in gut metagenomes, particularly focusing on mother-infant studies. It integrates the InStrain tool with other bioinformatics tools to provide insights into strain diversity and transmission.

## Script
```bash

### To generate a scaffold to bin file which lists the genome assignment of each scaffold
### In Drep program

#module load drep/3.2.0

#parse_stb.py --reverse -f /data/Food/analysis/R0372_MicrobeMom/Sai/metawrap_vs_metabat/drep/metawrap/swedish/output/dereplicated_genomes/*  -o genomes.stb 

##### Change above -f to the dereplicated genomes folder ####

#module unload drep/3.2.0

#module load instrain/1.8.0

### Running prodigal to get the genes for the representative genomes

#current_dir=$(pwd)
#mkdir genes

#echo "Copying all the drep genomes into one genomes rep file"

#cat /data/Food/analysis/R0372_MicrobeMom/Sai/metawrap_vs_metabat/drep/metawrap/swedish/output/dereplicated_genomes/* > ./genomes_reps.fasta

#cd /data/Food/analysis/R0372_MicrobeMom/Sai/metawrap_vs_metabat/drep/metawrap/swedish/output/dereplicated_genomes    

#### CHange this to the derpelicated genomes directory                                  

#echo "Started prodigal parallely on genomes"

#ls *fa | parallel -k --max-args=1 -j 20 "prodigal -i {1} -o $current_dir/genes/{1}.genes -a $current_dir/genes/{1}.gene.faa -d $current_dir/genes/{1}.gene.fna  -m -p single"

### Concatenate all the genomes,genes

#cd  $current_dir

#cat genes/*.gene.fna > genome_reps.genes.fna

#cat genes/*.gene.faa > genome_reps.genes.faa


###### Creating bowtie-database

#echo "Creating bowtie-database"

#mkdir bowtie-database

#bowtie2-build ./genomes_reps.fasta ./bowtie-database/genomes_reps.fasta.bt2 --threads 20 --large-index

###### Mapping to rep genome database

#mkdir bowtie-mapping
#mkdir instrain-profiles

for x in $(cat mother_samples.txt );             ####### CHANGe the input file names file

do	

echo  $x "Copying the reads from metawrap directory and unzipping"

cp  /data/Food/analysis/R0372_MicrobeMom/Sai/microbemom/public_datasets/swedish_dataset_1year/metawrap_mom_samples/output/"$x"/READ_QC/*gz . ### CHANGE THE PATH TO READ FILES

tar -xzvf qced_host_removed_reads.fastq.tar.gz

if [ -e "$x" ]; then

        echo $x "exits in normal id folder"
	mv "$x"/READ_QC/final_pure_reads_1.fastq ./"$x"_1.fastq
	mv "$x"/READ_QC/final_pure_reads_2.fastq ./"$x"_2.fastq

else
        echo $x "exists in output folder"
	mv ./output/"$x"/READ_QC/final_pure_reads_1.fastq ./"$x"_1.fastq
	mv ./output/"$x"/READ_QC/final_pure_reads_2.fastq ./"$x"_2.fastq
	
fi


###  -1 -2 fastq paths from the below

echo $x "Mapping to bowtie"

# Mapping Reads with bowtie2
module load bowtie2/2.4.4

# Mapping Reads with bowtie2
bowtie2 -p 20 -x bowtie-database/genomes_reps.fasta.bt2 -1 ./"$x"_1.fastq -2 ./"$x"_2.fastq > ./bowtie-mapping/"$x".sam

# Mapping Reads with bowtie2
module unload bowtie2/2.4.4

echo "Deleting the direcoty of metawrap reads and fastq files copied"

rm -r "$x" ./"$x"_1.fastq ./"$x"_2.fastq qced_host_removed_reads.fastq.tar.gz ./output

module load samtools/1.10

samtools view -bS ./bowtie-mapping/"$x".sam > ./bowtie-mapping/"$x".bam -@ 20

samtools sort ./bowtie-mapping/"$x".bam -o ./bowtie-mapping/sorted_"$x".bam -@ 20

module unload samtools/1.10

########## InSTRAIN profile

echo "Instrain profiling started" $x

module load instrain/1.8.0
module load samtools/1.10

# Running InStrain for Microbial Strain Profiling
inStrain profile ./bowtie-mapping/sorted_"$x".bam genomes_reps.fasta -o ./instrain-profiles/"$x"_instrain_profile -p 20 -g genome_reps.genes.fna  -s genomes.stb --database_mode

module unload instrain/1.8.0
module unload samtools/1.10

done

######## InStrain compare

inStrain compare -i ./merged_instrain-profiles/* -s genomes.stb -p 25 -o instrain_compare --database_mode


######## InStrain parse

inStrain parse_annotations -i ./merged_instrain-profiles/* -a ./gene_annotations/eggnog/pfam_annotation_table.csv -p 25 

module unload instrain/1.8.0

```
