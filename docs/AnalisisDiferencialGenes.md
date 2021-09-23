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
bamfilesSE <- list.files(path = "bamfile_SE", 
                       pattern = "*.bam$", 
                       full.names = TRUE )

bamfilesPE <- list.files(path = "bamfile_PE", 
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
conteo$counts
write.table(conteo$counts , file = "rawSE_counts.txt", sep = "\t")
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
