Análisis Diferencial de Genes
================

## Introducción



## CONTENIDO

- [CONTEO DE GENES](#conteo-de-genes)
- [CONTROL DE CALIDAD](#control-de-calidad)
- [FILTRADO DE SECUENCIAS](#filtrado-de-secuencias)
- [FILTRADO DE rRNA](#filtrado-de-rrna)
- [ALINEAMIENTO](#alineamiento)
- [CONTEO DE SECUENCIAS](#conteo-de-secuencias)
- [IMPORTAR A R](#importar-a-r)

## CONTEO DE GENES
 
 ``` r
library(DESeq2)
library(ggplot2)
library(clusterProfiler)
library(biomaRt)
library(ReactomePA)
library(DOSE)
library(KEGG.db)
library(org.Hs.eg.db)
library(pheatmap)
library(genefilter)
library(RColorBrewer)
library(GO.db)
library(topGO)
library(dplyr)
library(gage)
library(ggsci)
```


``` r
bamfilesSE <- list.files(path = "Star_bam/SE", 
                       pattern = "*.bam$", 
                       full.names = TRUE )

bamfilesPE <- list.files(path = "Star_bam/PE", 
                       pattern = "*.bam$", 
                       full.names = TRUE )
                       
bamfilesSE
```


```r
conteo <- featureCounts( files = bamfilesSE, 
                      annot.ext = "GCF_000001735.4_TAIR10.1_genomic.gff", 
                      isGTFAnnotation = TRUE, 
                      GTF.featureType = "gene", 
                      GTF.attrType = "ID", 
                      isPairedEnd = FALSE, 
                      requireBothEndsMapped = FALSE, 
                      minMQS = 20, 
                      strandSpecific = 2 )
```
           
                             ==========     _____ _    _ ____  _____  ______          _____  
                             =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
                               =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
                                 ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
                                   ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
                             ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
                            Rsubread 2.6.4
                     
                     //========================== featureCounts setting ===========================\\
                     ||                                                                            ||
                     ||             Input files : 10 BAM files                                     ||
                     ||                                                                            ||
                     ||                           SRR390302_sortedByCoord.out.bam                  ||
                     ||                           SRR390303_sortedByCoord.out.bam                  ||
                     ||                           SRR390304_sortedByCoord.out.bam                  ||
                     ||                           SRR390306_sortedByCoord.out.bam                  ||
                     ||                           SRR390307_sortedByCoord.out.bam                  ||
                     ||                           SRR390308_sortedByCoord.out.bam                  ||
                     ||                           SRR390309_sortedByCoord.out.bam                  ||
                     ||                           SRR390311_sortedByCoord.out.bam                  ||
                     ||                           SRR390312_sortedByCoord.out.bam                  ||
                     ||                           SRR390313_sortedByCoord.out.bam                  ||
                     ||                                                                            ||
                     ||              Paired-end : no                                               ||
                     ||        Count read pairs : no                                               ||
                     ||              Annotation : GCF_000001735.4_TAIR10.1_genomic.gff (GTF)       ||
                     ||      Dir for temp files : .                                                ||
                     ||                 Threads : 1                                                ||
                     ||                   Level : meta-feature level                               ||
                     ||      Multimapping reads : counted                                          ||
                     || Multi-overlapping reads : not counted                                      ||
                     ||   Min overlapping bases : 1                                                ||
                     ||                                                                            ||
                     \\============================================================================//
   
```r
data$counts
write.table(data$counts , file = "rawSE_counts.txt", sep = "\t")
```

   
```r
data <- read.table("rawSE_counts.txt", header = TRUE , row.names = 1)
head(data)
```

```r
colnames(data) <- gsub("_sortedByCoord.out.bam", "", colnames(data), fixed = T)
colnames(data) <- gsub("..", "", colnames(data), fixed = T)
row.names(data) <- gsub("gene-", "", rownames(data), fixed = T)
head(data)
```

    
```r
metadata <- read.delim("muestras.txt", row.names = 1)
metadata$sampleid <- row.names(metadata)
metadata <- metadata[match(colnames(data), metadata$sampleid),]
metadata
```
```r
data$counts
write.table(data$counts , file = "rawSE_counts.txt", sep = "\t")
```

```r
ddsMat <- DESeqDataSetFromMatrix(countData = data, 
                                 colData = metadata, 
                                 design = ~Group)
ddsMat <- DESeq(ddsMat)
```

```r
# Get results from testing with FDR adjust pvalues
results <- results(ddsMat, pAdjustMethod = "fdr", alpha = 0.05)

# Generate summary of testing. 
summary(results)
```
```r
mcols(results, use.names = T)
results
```
```r

```
```r

```
```r

```
```r

```
```r

```
```r

```
```r

```


## CONTROL DE CALIDAD

```bash
── results/multiQC/
    └── multiqc_report.html     <-  Estadisticas de calidad de todos los archivos fastqc (.html) 
    └── multiqc_data/           <-  Datos que multiqc encontró de varios archivos de registro (.zip)
```
    
## FILTRADO DE SECUENCIAS

```bash 

#trim_galore 
    trim_galore --quality 20 --fastqc --length 25 --output_dir results/trimmed/ reads/SRR849504/NC_000021.9_1_trimmed.fastq
    
```

## FILTRADO DE rRNA

```bash

```

## ALINEAMIENTO

### Es puede tomar con los datos completos alrededor de 30 minutos
### Para el curso utilizaremos unicamente el cromosoma 21 [NC_000021.9](https://www.ncbi.nlm.nih.gov/nuccore/NC_000021.9)    

Antes de empezar con el alineamiento se va desargar el archivo de anotación gff3 del cromosoma 21 y se utilizará el programa **gffread** para convertir al formato gtf, para mas detalles revisar el repositorio del programa [aqui](https://github.com/gpertea/gffread)

```bash

```

#### Indexar el genoma con gtf-file
#### Descargamos primero el archivo gff3 del NCBI

```bash

```
#### Convertimos el archivo gff3 a gtf
```bash

```
### Indexar el genoma con STAR

| Parámetro | Descripción |
| ---- | ---- |
| `--runMode` | indica el tipo de opcion que utilizará STAR, en este caso queremos generar un índice del genoma por lo que utilizamos la flag `genomeGenerate`|
| `--genomeDir` | indica donde se guardaran los resultados del indice y la ubicación de los archivos del genoma |
| `--genomeFastaFiles` | indica donde estan almacenadas las secuencias del genoma en formato `FASTA` |
| `--sjdbGTFfile` | sj: splice junction db: database GTFfile: archivo GTF indica la ubicación del archivo GTF para mejorar e improvisar el mapeo dado el modelo de los genes |
| `--sjdbOverhang` | Especifica el largo a considerar de la secuencia genómica alrededor del splice junction, este valor esta ligado al largo de los reads y deberia ser `max(ReadLength) - 1`  |
| `--runThreadN` | total de hebras que se ejecutaran en paralelo, este número no debe sobrepasar la cantidad de cores que tiene un computador y pruebas de escalamiento deberian ser ejecutadas para calcular el óptimo |

    STAR \
    --runMode genomeGenerate \
    --genomeDir genome/star_index \
    --genomeFastaFiles genome/NC_000021.9.fna \
    --sjdbGTFfile annotation/NC_000021.9.gtf \
    --runThreadN 2

``` bash 
STAR --runMode genomeGenerate --genomeDir genome/star_index/ --genomeSAindexNbases 7 --genomeFastaFiles genome/NC_000932.1.fasta --sjdbGTFfile genome/NC_000932.1.gtf --runThreadN 2
```

### Alinear el genoma

| Parámetro | Descripción |
| ---- | ---- |
| `--readFilesIn` | archivo de reads a mapear |
| `--genomeDir` | indica donde esta alojado el genoma |
| `--runThreadN`| cantidad de threads |
| `--outSAMType`| tipo de archivo `SAM` o `BAM` en este caso indicamos `BAM`|
| `SortedByCoordinate`| Ordenar el archivo BAM por las coordenadas del genoma |
| `--outFileNamePrefix`| Prefijo para los archivos de salida |

#### Command

    # Help
    STAR -h

    # Run STAR (~3min)
    STAR \
    --genomeDir genome/star_index \
    --readFilesIn results/trimmed/sample_filtered.fq  \
    --runThreadN 2 \
    --outSAMtype BAM SortedByCoordinate \
    --quantMode GeneCounts

```bash
STAR --genomeDir ../../../genome/star_index/ --readFilesIn ../../trimmed/NC_000932.1_6_trimmed.fq --runThreadN 2 --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 12000000000 --quantMode GeneCounts
```
 
``` bash

```

## CONTEO DE SECUENCIAS

    featureCounts -h
    
```bash

```


``` bash

```


## IMPORTAR A R



