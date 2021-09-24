7.Enrriquecimiento de genes
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
library(org.At.tair.db)
library(GO.db)
library(topGO)
library(clusterProfiler)
library(pathview)
```


``` r
columns(org.At.tair.db)
keytypes(org.At.tair.db)
```


```r
# Add gene full name
results$description <- mapIds(x = org.At.tair.db,
                              keys = row.names(results),
                              column = "GENENAME",
                              keytype = "TAIR",
                              multiVals = "first")

# Add gene symbol
results$symbol <- row.names(results)

# Add ENTREZ ID
results$entrez <- mapIds(x = org.At.tair.db,
                         keys = row.names(results),
                         column = "ENTREZID",
                         keytype = "TAIR",
                         multiVals = "first")

results$evidence <- mapIds(x = org.At.tair.db,
                         keys = row.names(results),
                         column = "EVIDENCE",
                         keytype = "TAIR",
                         multiVals = "first")

results$GO <- mapIds(x = org.At.tair.db,
                         keys = row.names(results),
                         column = "GO",
                         keytype = "TAIR",
                         multiVals = "first")
```

```r
results_sig <- subset(results, padj < 0.05)
results_sig
```

   
```r
# Remove any genes that do not have any entrez identifiers
results_sig_entrez <- subset(results_sig, is.na(GO) == FALSE)

# Create a matrix of gene log2 fold changes
# View the format of the gene matrix
##- Names = ENTREZ ID
##- Values = Log2 Fold changes
results_sig_entrez
```

```r
columns(GO.db)
keytypes(GO.db)
enrich <- AnnotationDbi::select(GO.db, keys = results_sig$GO, columns = c('TERM', 'ONTOLOGY'), keytypes = 'GO')
enrich$symbol <- results_sig$symbol
enrich
```

    
```r
enrich2 <- AnnotationDbi::select(GO.db, keys = results_sig$GO, columns = 'DEFINITION', keytypes = 'GO')
enrich2$symbol <- results_sig$symbol
enrich2
```
```r
gene_matrix <- results_sig_entrez$log2FoldChange
names(gene_matrix) <- results_sig_entrez$entrez
```

```r
pathview(gene.data = gene_matrix, 
         pathway.id = "04130", 
         species = "ath")

```


[tair link](https://www.arabidopsis.org/tools/go_term_enrichment.jsp)

[PHANTER link](http://pantherdb.org/)
