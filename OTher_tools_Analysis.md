# Other tools used in the study

## Description
This script includes a mix of tools used for various analyses in the study. It integrates different scripts, each tailored for specific aspects analysis.

## Analysis Steps and Script

### Antimicrobal resistnance genes identification using Abricate 1.0.1
```bash
abricate mags_directory --csv --db card --threads 15  >> amr_output.csv
```

### 16S rRNA gene prediction in MAGs using Barrnap 0.9
```bash
for x in $(ls MAGs*);
do
```

```bash
barrnap --threads 5 each_mag.fa > "$x".gff
```

```bash
done
```

### Taxonomy for MAGs was assigned using gtdbtk 1.5.0
```bash
gtdbtk classify_wf --genome_dir mags_directory --out_dir output_directory --cpus 25 --extension fa
```

### Representative MAGs for InStrain database - drep 3.2.0
```bash
dRep dereplicate output_directory  -p 10 -g MAGs* -comp 90 -con 5 -sa 0.98```
