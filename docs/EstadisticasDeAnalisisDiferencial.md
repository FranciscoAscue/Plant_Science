6.Estadísticas del Análisis Diferencial de Genes
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
library(BiocManager)
library(Rsubread)
library(Rbowtie2)
library(Rsamtools)
library(dplyr)
library(org.At.tair.db)
library(pheatmap)
library(GO.db)
library(pathview)
library(KEGGgraph)
library(DESeq2)
library(ggplot2)
library(KEGGREST)
```


``` r
# Convert all samples to rlog
ddsMat_rlog <- rlog(ddsMat, blind = FALSE)

# Plot PCA by column variable
plotPCA(ddsMat_rlog, intgroup = "Group", ntop = 7) +
  theme_bw() + # remove default ggplot2 theme
  geom_point(size = 5) + # Increase point size
  scale_y_continuous(limits = c(-5, 5)) + # change limits to fix figure dimensions
  #scale_color_brewer(palette = "Dark2")+  
  ggtitle(label = "Principal Component Analysis (PCA)", 
          subtitle = "Top 7 most variable genes") 
```


```r
results_sig <- subset(results, padj < 0.05)
results_sig
```
           
                  
```r
library(RColorBrewer)
library(viridisLite)
library(viridis)
```

   
```r
ddsMat_rlog <- rlog(ddsMat, blind = TRUE)
mat <- assay(ddsMat_rlog[row.names(results_sig)])[1:18, ]

##


annotation_col = data.frame(
  Group = factor(colData(ddsMat_rlog)$Group), 
  row.names = colData(ddsMat_rlog)$sampleid
)

##



ann_colors = list(
  Group = c("0day" = "#E65100", "1day" = "#66BB6A", "3day" = "#03A9F4"), 
  sampleid= c(a = "red")  
)

##

pheatmap(mat = mat, 
         color = viridis(180), 
         scale = "row", # Scale genes to Z-score (how many standard deviations)
         annotation_col = annotation_col, # Add multiple annotations to the samples
         annotation_colors = ann_colors,# Change the default colors of the annotations
         fontsize = 7, # Make fonts smaller
         cellwidth = 30, # Make the cells wider
         cutree_cols = 3,
         cutree_rows = 4,
         show_colnames = F)


pheatmap(mat = mat, 
         color = viridis(180), 
         scale = "row", # Scale genes to Z-score (how many standard deviations)
         annotation_col = annotation_col, # Add multiple annotations to the samples
         annotation_colors = ann_colors,# Change the default colors of the annotations
         fontsize = 7, # Make fonts smaller
         cellwidth = 30, # Make the cells wider
         cutree_cols = 3,
         cutree_rows = 4,
         cluster_cols = F,
         show_colnames = F)

```

```r
# Gather Log-fold change and FDR-corrected pvalues from DESeq2 results
## - Change pvalues to -log10 (1.3 = 0.05)
data_vp <- data.frame(gene = row.names(results),
                   pval = -log10(results$padj), 
                   lfc = results$log2FoldChange)

# Remove any rows that have NA as an entry
data_vp <- na.omit(data_vp)

# Color the points which are up or down
## If fold-change > 0 and pvalue > 1.3 (Increased significant)
## If fold-change < 0 and pvalue > 1.3 (Decreased significant)
data_vp <- mutate(data_vp, color = case_when(data_vp$lfc > 0 & data_vp$pval > 1.3 ~ "Increased",
                                       data_vp$lfc < 0 & data_vp$pval > 1.3 ~ "Decreased",
                                       data_vp$pval < 1.3 ~ "nonsignificant"))

# Make a basic ggplot2 object with x-y values
vol <- ggplot(data_vp, aes(x = lfc, y = pval, color = color))

# Add ggplot2 layers
vol +   
  ggtitle(label = "Volcano Plot", subtitle = "Colored by fold-change direction") +
  geom_point(size = 2.5, alpha = 0.8, na.rm = T) +
  scale_color_manual(name = "Directionality",
                     values = c(Increased = "#B8DE29FF", Decreased = "#482677FF", nonsignificant = "darkgray")) +
  theme_bw(base_size = 14) + # change overall theme
  theme(legend.position = "right") + # change the legend
  xlab(expression(log[2]("LoGlu" / "HiGlu"))) + # Change X-Axis label
  ylab(expression(-log[10]("adjusted p-value"))) + # Change Y-Axis label
  geom_hline(yintercept = 1.3, colour = "darkgrey") + # Add p-adj value cutoff line
  scale_y_continuous(trans = "log1p") # Scale yaxis due to large p-v
```

    
```r
plotMA(results, ylim = c(-5, 5))
plotDispEsts(ddsMat)
```
```r
# Convert all samples to rlog
ddsMat_rlog <- rlog(ddsMat, blind = FALSE)

# Get gene with highest expression
top_gene <- rownames(results)[which.min(results$log2FoldChange)]

# Plot single gene
plotCounts(dds = ddsMat, 
           gene = top_gene, 
           intgroup = "Group", 
           normalized = T, 
           transform = T)
```
