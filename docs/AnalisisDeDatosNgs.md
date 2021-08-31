Análisis de datos NGS
=====================
...

<p align="center">
    <img width="40%" src="https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/img/RNAseqWorkflow.png">
</p>

## CONTENIDO

- [CONTROL CALIDAD](#control-calidad)
- [FILTRADO DE READS](#filtrado-de-reads)
- [ALINEAMIENTO](#alineamiento)

## CONTROL CALIDAD

```bash
fastqc -t 2 cavtsc_forward_paired.fq.gz cavtsc_reverse_paired.fq.gz -o /mnt/disco2/fascue/cporcellus/results/fastqc/
```
[fastqc link](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

<p align="center">
    <img width="60%" src="https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc.png">
</p>

## FILTRADO DE READS

<p align="center">
    <img width="70%" src="https://usermanual.wiki/Document/TrimmomaticManualV032.1972804677/asset-6.png">
</p>

```
java -jar ../descargas/Trimmomatic-0.39/trimmomatic-0.39.jar PE
     -phred33 
     -threads 2 
     file_1.fastq file_2.fastq 
     file_forward_paired.fq.gz file_forward_unpaired.fq.gz 
     file_reverse_paired.fq.gz file_revers_unpaired.fq.gz 
     ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 
     LEADING:3 
     TRAILING:3 
     SLIDINGWINDOW:4:25 
     MINLEN:25
```

```sh
java -jar ../descargas/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 -threads 2 file_1.fastq file_2.fastq file_forward_paired.fq.gz file_forward_unpaired.fq.gz file_reverse_paired.fq.gz file_revers_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:25 MINLEN:25
```

```   
java -jar ../descargas/Trimmomatic-0.39/trimmomatic-0.39.jar SE
     -phred33 
     -threads 2 
     file.fastq 
     file_trimm.fq 
     ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 
     LEADING:3 
     TRAILING:3 
     SLIDINGWINDOW:4:25 
     MINLEN:25

```

```
ILLUMINACLIP:<fastaWithAdaptersEtc>:<seed mismatches>:<palindrome clip threshold>:<simple clip threshold>
LEADING:<quality> 
TRAILING:<quality> 
SLIDINGWINDOW:<windowSize>:<requiredQuality> 
MINLEN:<length>

```

## ALINEAMIENTO

<p align="center" width="100%">
    <img width="50%" src="https://i.ytimg.com/vi/6BJbEWyO_N0/maxresdefault.jpg">
</p>

```r 
#cargar paquetes para R
library(Rbowtie2)
library(Rsamtools)
library(ape)
library(viridisLite)
library(viridis)
```

Descargar secuencia de referencia [Genoma de A. thaliana](https://www.ncbi.nlm.nih.gov/genome/?term=Arabidopsis%20thaliana) , para tener un panorama del organismo en el que estamos trabajando se debe analizar el genoma anotado, para ello podemos trabajar en R con archivos `gff`

```r
# read gff files with ape

gff_file <- read.gff("sequence.gff3", na.strings = c(".", "?"), GFF3 = TRUE)

# transform in matrix to filter annotations

tab <- as.matrix(table(gff_file$type))
rnames <- as.matrix(rownames(tab))
```

```r
# make a plot of annotation feactures

etiquetas <- paste0(rnames[c(4,7,10,16),],"=",round(100 * tab[c(4,7,10,16),]/sum(tab[c(4,7,10,16),]), 2), "%")                                                       
par(mfrow=c(1,2), adj = TRUE)
pie(tab[c(4,7,10,16),], labels = etiquetas, col = viridis(4))
#pie(tab[c(4,7,10,16),], col = viridis(4), labels = paste0(tab[c(4,7,10,16),], "%"))
barplot(tab[c(4,7,10,16),], col = viridis(4), width = 60)

```
### Preparacion del index 

```r
bowtie2_build("AthalianaChr4.fasta", bt2Index = "index/" , overwrite = TRUE)
```
### Alineamiento de secuencias

```r
bowtie2_build("AthalianaChr4.fasta", bt2Index = "index/" , overwrite = TRUE)

bowtie2(bt2Index = "index/", 
        samOutput = "SRR390310.sam", 
        seq1 = "SRR390310_1.fastq", 
        seq2 = "SRR390310_2.fastq", 
        "--threads=3")
```

### Convertir SAM a BAM

```r
asBam("SRR390310.sam")
```

### Visualizar alineamiento

[igv link](https://software.broadinstitute.org/software/igv/download)
