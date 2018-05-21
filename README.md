# Avian Plasmid Read-Mapping Plot

NOTE: This script has been modified from a pre-existing script available at <https://github.com/cerebis/antarctic_ha>.


It contains scripts supporting an investigation into Avian Pathogenic Escherichia coli.


Paper available at: ~~URL-to-paper-here~~



The python script 'analysis.py' was used to produce the described heatmap:


## BWA Index
The reference used was pCERC4 (Accession No. KU578032.1).
```
bwa index pCERC4_KU578032.fasta
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
python ~/readdepth_plot/readdepth_plot.py -b <binsize> -t <ticksize> pCERC4_KU578032.fasta Heatmap_plot.pdf
```
# Python and Package versions
* Python 2.7
* bwa - 0.7.17
* matplotlib - 2.1.0
* numpy - 1.13.3
* pandas - 0.20.3
* samtools - 0.1.18
* scipy - 0.19.1
