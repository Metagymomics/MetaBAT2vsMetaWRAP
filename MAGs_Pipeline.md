# Metagenome Assembled Genomes (MAGs) Pipeline

## Description
This script is designed for generating Metagenome Assembled Genomes (MAGs) from gut metagenomes using METAWRAP.

## Pipeline Steps and Script

### metawrap/1.2.1
### samtools/1.10
### READS INPUT FILE NAME SHOULD BE IN THIS FORMAT name_1.fastq
```bash
mkdir output

for x in $(cat input_samples_list);
do 

	cp ../raw_fastq/"$x"_1.fastq.gz ../raw_fastq/"$x"_2.fastq.gz  . 
	gunzip "$x"_1.fastq.gz	"$x"_1.fastq
	gunzip "$x"_2.fastq.gz  "$x"_2.fastq
	mkdir output/"$x"  
	echo "STARTED QUALITY FILTERING"
	metawrap read_qc -1 ./"$x"_1.fastq -2 ./"$x"_2.fastq -t 20 -o ./output/"$x"/READ_QC;
	echo "DONE QC"
	echo "STARTED ASSEMBLING WITH METASPADES"
	metawrap assembly -m 800 -t 20 -1 ./output/"$x"/READ_QC/final_pure_reads_1.fastq -2 ./output/"$x"/READ_QC/final_pure_reads_2.fastq  --metaspades -o ./output/"$x"/ASSEMBLY
	echo "DONE METASPADES"
	echo "STARTED BINNING WITH THREE BINNING TOOLS"
	metawrap binning -o ./output/"$x"/INITIAL_BINNING -t 20 -m 800 -a ./output/"$x"/ASSEMBLY/final_assembly.fasta --metabat2 --maxbin2 --concoct  ./output/"$x"/READ_QC/final_pure*
	echo "DONE BINNING NOW REFINING"


### # checkm/1.0.18
	metawrap bin_refinement -o ./output/"$x"/BIN_REFINEMENT  -m 800 -t 20 -A ./output/"$x"/INITIAL_BINNING/metabat2_bins/ -B ./output/"$x"/INITIAL_BINNING/maxbin2_bins/ -C ./output/"$x"/INITIAL_BINNING/concoct_bins/ -c 50 -x 10
	echo "DONE REFINING NOW REASSMBLING"
	metawrap reassemble_bins -t 20 -m 800 -o ./output/"$x"/BIN_REASSEMBLY -1 ./output/"$x"/READ_QC/final_pure_reads_1.fastq -2 ./output/"$x"/READ_QC/final_pure_reads_2.fastq  -c 50 -x 10 -b ./output/"$x"/BIN_REFINEMENT/metawrap_50_10_bins
	echo "METAWRAP DONE on $x"
	rm "$x"_1.fastq "$x"_2.fastq

### Compressing the QC and host contamination removed reads
tar -cvzf ./output/"$x"/READ_QC/qced_host_removed_reads.fastq.tar.gz ./output/"$x"/READ_QC/final_pure_reads_*

### Cleaning up temp files
rm -r ./output/"$x"/ASSEMBLY/metaspades/ ./output/"$x"/ASSEMBLY/QUAST_out/ ./output/"$x"/READ_QC/final_pure_reads_*	 ./output/"$x"/READ_QC/host_reads_* ./output/"$x"/READ_QC/*report* ./output/"$x"/BIN_REFINEMENT/ ./output/"$x"/INITIAL_BINNING/ ./output/"$x"/BIN_REASSEMBLY/work_files/

done
```
