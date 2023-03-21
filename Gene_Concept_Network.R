#---- Load Packages ------------------------------------------------------------
library(ggplot2)
library(ggrepel)
library(clusterProfiler)
library(readr)
library(dplyr)
library(enrichplot)
library(org.Mm.eg.db)

#---- Set Directory ------------------------------------------------------------
out_dir <- ''

#---- GO: All Clusters ---------------------------------------------------------
df <- read.csv(paste0(out_dir,'aggr_AAP_vs_YAP','.csv'), row.names = 'X') %>%
  filter(pct.1 > 0.1, pct.2 > 0.1)
df$gene <- row.names(df)
# up-regulated
df_up <- df %>%
  filter(p_val_adj < 0.05, avg_log2FC > 0.25)
ego <- enrichGO(gene = df_up$gene,
                keyType = "SYMBOL",
                OrgDb = org.Mm.eg.db,
                ont = 'BP',
                minGSSize = 10,
                maxGSSize = 500,
                pAdjustMethod = "BH",
                pvalueCutoff = 1)

df_up_select <- df_up %>%
  filter(avg_log2FC > 2)
geneList <- df_up_select$avg_log2FC
names(geneList) <- df_up_select$gene
# only keep genes with avg_log2FC > 2
n = 1
for(i in ego@result$geneID)
{
  genes <- ''
  for (j in as.list(strsplit(i,'/'))[[1]])
  {
    if (j %in% names(geneList))
    {
      genes <- paste0(genes,'/',j)
    }
  }
  ego@result$geneID[n] <- substring(genes,2)
  ego@result$Description[n] <- n
  n = n + 1
}

#---- Gene-Concept Network -----------------------------------------------------
# top 12 pathways
p <- cnetplot(ego, color.params = list(foldChange = geneList), showCategory = 12,
              cex.params = list(category_label = 2, gene_label = 1.5)) +
  scale_colour_gradient(low = "yellow", high = "red", name = "log2FC", breaks <- c(3,4)) +
  theme_bw() +
  theme(legend.text = element_text(size = 35),
        legend.title = element_text(size = 35),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank())
ggsave(paste0(out_dir,'ORA/figure/','gene_concept_network_top12.png'),width=18,height=14)

# top 5 pathways
p <- cnetplot(ego, color.params = list(foldChange = geneList), showCategory = 5,
              cex.params = list(category_label = 2, gene_label = 1.5)) +
  scale_colour_gradient(low = "yellow", high = "red", name = "log2FC", breaks <- c(3,4)) +
  theme_bw() +
  theme(legend.text = element_text(size = 35),
        legend.title = element_text(size = 35),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank())
ggsave(paste0(out_dir,'ORA/figure/','gene_concept_network_top5.png'),width=18,height=14)


