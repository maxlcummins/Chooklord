# Avian Plasmid Read-Mapping Plot

**NOTE: This script has been modified from a pre-existing script available at <https://github.com/cerebis/antarctic_ha>**



It contains scripts supporting an investigation titled - "Whole genome sequence analysis of Australian avian pathogenic Escherichia coli that carry the class 1 integrase gene"



Paper available [here](https://mgen.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000250?fbclid=IwAR2WYZFZVR3B_h4207ndDaOoqoh48oxRt1iXoMUK2SFD7CKlgp7MSa9GIuE)




The python script 'chooklord.py' was used to produce the described heatmap:


## BWA Index
The reference used was pCERC4 (Accession No. KU578032.1).
```
bwa index pCERC4_KU578032.fasta
```
## Array definition
An array was defined as below to facilitate the subsequent step
```
array=($(ls read_directory/*.fastq.gz))
```
## BWA mem
A for loop was used to feed read pairs into an executable 'bwa.qsub'.
```
for ((i=0; i<${#array[@]}; i+=2));
do qsub -v REF=pCERC4_KU578032.fasta,R1=${array[i]},R2=${array[i+1]},OUT=${array[i]%R1_001.fastq.gz} bwa.qsub;
done
```
'bwa.qsub' contained, along with the appropriate job submission syntax, the following
```
bwa mem -t16 -MY $REF $R1 $R2 | samtools view -ubS -F 0x904 - | samtools sort -@8 -T $REF - -o ${R1}.bam           
```
## SAMtools depth 
```
for fn in `ls ./*.bam`;
do samtools depth ${fn} > ${fn%.bam}_coverage.txt;
done
```
## Heatmap generation
Here the binsize and ticksize were 250 and 100, respectively.

```
python chooklord.py -b <binsize> -t <ticksize> pCERC4_KU578032.fasta Heatmap_plot.pdf
```
# Python and Package versions
* Python 2.7
* bwa - 0.7.17
* matplotlib - 2.1.0
* numpy - 1.13.3
* pandas - 0.20.3
* samtools - 0.1.18
* scipy - 0.19.1

# Acknowledgements
Many thanks to Matt DeMaere for his efforts in putting together the original script and helping to adapt and repurpose it for our investigation.

# License
Copyright (c) 2017 Matthew DeMaere

See [LICENSE](../../blob/master/LICENSE) for more information
