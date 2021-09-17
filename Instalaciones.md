Instalaciones del curso
===========
Antes de empezar con el curso instale los siguientes paquetes:

#### FastQC: quality control for high throughput sequence data [Link fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/){: .btn }

#### Trimmomatic: A flexible read trimming tool for Illumina NGS data [Link trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic){: .btn }  

#### SRA Toolkit:The SRA Toolkit and SDK from NCBI is a collection of tools and libraries for using data in the INSDC Sequence Read Archives [Link sratoolkit](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software){: .btn-purple }

#### R packages: Instalaciones de paquetes para el análisis de RNASeq

Para instalar los paquetes de R usaremos el gestor de paquetes de Bioconductor.

```r
## Instalar gestor de paquetes de Bioconductor
install.packages("BiocManager")

## Activar BiocManager
library(BiocManager)
```

```r
## Instalar los siguientes paquetes:

BiocManager::install("Rbowtie2")
BiocManager::install("Rsubread")
BiocManager::install("Rsamtools")
BiocManager::install("dplyr")
BiocManager::install("org.At.tair.db")
BiocManager::install("pheatmap")
BiocManager::install("GO.db")
BiocManager::install("DESeq2")
BiocManager::install("ggplot2")
BiocManager::install("pathview")
BiocManager::install("KEGGgraph")

```  
