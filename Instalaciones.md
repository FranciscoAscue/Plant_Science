## 2.Instalaciones del curso

Antes de empezar con el curso instale los siguientes programas y paquetes, estos links son para el SO windows.

#### FastQC: quality control for high throughput sequence data [_Descargar fastqc_](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip){: .btn-blue}

#### Trimmomatic: A flexible read trimming tool for Illumina NGS data [_Link trimmomatic_](http://www.usadellab.org/cms/?page=trimmomatic){: .btn-green }  

#### SRA Toolkit:The SRA Toolkit and SDK from NCBI [_Link sratoolkit_](https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.11.1/sratoolkit.2.11.1-win64.zip){: .btn-purple }

### Instalacion de R y Rstudio para Windows.


[Descargar R ](https://cran.r-project.org/src/base/R-4/R-4.1.1.tar.gz){: .btn } [Descargar Rstudio ](https://download1.rstudio.org/desktop/windows/RStudio-1.4.1717.exe){: .btn } 

## Instalaciones de paquetes para el análisis de RNASeq

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
BiocManager::install("ape")
BiocManager::install("ggplot2")
BiocManager::install("pheatmap")
BiocManager::install("pathview")
BiocManager::install("clusterProfiler")
BiocManager::install("biomaRt")
BiocManager::install("genefilter")

```  
