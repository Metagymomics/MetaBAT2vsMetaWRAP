# InStrain Script for Gut Metagenome Analysis

## Description
This script is used for strain-level analysis in gut metagenomes, particularly focusing on mother-infant studies. It integrates the InStrain tool with other bioinformatics tools to provide insights into strain diversity and transmission.

## Detailed Steps and Script

###  To generate a scaffold to bin file which lists the genome assignment of each scaffold

###  In Drep program 3.2.0

```bash
parse_stb.py --reverse -f dereplicated_genomes_folder  -o genomes.stb 
```

###  Change above -f to the dereplicated genomes folder 

### instrain/1.8.0

### Running prodigal to get the genes for the representative genomes

```bash
current_dir=$(pwd)
mkdir genes
```

```bash
echo "Copying all the drep genomes into one genomes rep file"
```

```bash
cat dereplicated_genomes_folder/* > ./genomes_reps.fasta
```

```bash
cd dereplicated_genomes_folder    
```

###  CHange this to the derpelicated genomes directory                                  

```bash
echo "Started prodigal parallely on genomes"
```

```bash
ls *fa | parallel -k --max-args=1 -j 20 "prodigal -i {1} -o $current_dir/genes/{1}.genes -a $current_dir/genes/{1}.gene.faa -d $current_dir/genes/{1}.gene.fna  -m -p single"
```

###  Concatenate all the genomes,genes

```bash
cd  $current_dir
```

```bash
cat genes/*.gene.fna > genome_reps.genes.fna
```

```bash
cat genes/*.gene.faa > genome_reps.genes.faa
```

###  Creating bowtie-database

```bash
echo "Creating bowtie-database"
```

```bash
mkdir bowtie-database
```

```bash
bowtie2-build ./genomes_reps.fasta ./bowtie-database/genomes_reps.fasta.bt2 --threads 20 --large-index
```

###  Mapping to rep genome database

```bash
mkdir bowtie-mapping
mkdir instrain-profiles
```

```bash
for x in $(cat input_samples);             ####### CHANGe the input file names file
```

```bash
do	
```

```bash
echo  $x "Copying the reads from metawrap directory and unzipping"
```

### ## CHANGE THE PATH TO FASTQ READS

```bash
cp  FASTQ_READS.gz  .
```

```bash
tar -xzf FASTQ_READS.gz
```

```bash
echo $x "Mapping to bowtie"
```

### ## bowtie2/2.4.4

```bash
bowtie2 -p 20 -x bowtie-database/genomes_reps.fasta.bt2 -1 ./"$x"_1.fastq -2 ./"$x"_2.fastq > ./bowtie-mapping/"$x".sam
```

```bash
echo "Deleting copied fastq files "
```

```bash
rm -r "$x" ./"$x"_1.fastq ./"$x"_2.fastq 
```

###  samtools/1.10

```bash
samtools view -bS ./bowtie-mapping/"$x".sam > ./bowtie-mapping/"$x".bam -@ 20
```

```bash
samtools sort ./bowtie-mapping/"$x".bam -o ./bowtie-mapping/sorted_"$x".bam -@ 20
```

###  InStrain profiling

```bash
echo "InStrain profiling started" $x
```

###  instrain/1.8.0

```bash
inStrain profile ./bowtie-mapping/sorted_"$x".bam genomes_reps.fasta -o ./instrain-profiles/"$x"_instrain_profile -p 20 -g genome_reps.genes.fna  -s genomes.stb --database_mode
```

```bash
done
```

###  InStrain comparing strain across all the profiled samples

```bash
inStrain compare -i ./instrain-profiles/* -s genomes.stb -p 25 -o instrain_compare --database_mode
```

###  InStrain parse

```bash
inStrain parse_annotations -i ./instrain-profiles/* -a ./pfam_annotation_table.csv -p 25 
```
