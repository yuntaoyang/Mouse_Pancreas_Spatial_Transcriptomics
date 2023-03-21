#---- Load Packages ------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(stringr)

#---- Set Directory ------------------------------------------------------------
in_dir <- ''
out_dir <- ''

#---- Define Functions ---------------------------------------------------------
# GO-BP
GO <- function(file)
{
  df <- read.csv(paste0(in_dir,file,'.csv'), row.names = 'X') %>%
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
  write.csv(ego,paste0(out_dir,'file/GO/',file,'_GO_Up.csv'),row.names = FALSE)
  # down-regulated
  df_down <- df %>%
    filter(p_val_adj < 0.05, avg_log2FC < -0.25)
  ego <- enrichGO(gene = df_down$gene,
                  keyType = "SYMBOL",
                  OrgDb = org.Mm.eg.db,
                  ont = 'BP',
                  minGSSize = 10,
                  maxGSSize = 500,
                  pAdjustMethod = "BH",
                  pvalueCutoff = 1)
  write.csv(ego,paste0(out_dir,'file/GO/',file,'_GO_Down.csv'),row.names = FALSE)
}

# KEGG
KEGG <- function(file)
{
  df <- read.csv(paste0(in_dir,file,'.csv'), row.names = 'X') %>%
    filter(pct.1 > 0.1, pct.2 > 0.1)
  df$gene <- row.names(df)
  ids <- bitr(df$gene,'SYMBOL','ENTREZID','org.Mm.eg.db')
  data <- merge(df,ids,by.x='gene',by.y='SYMBOL')
  # up-regulated
  data_up <- data %>%
    filter(p_val_adj < 0.05, avg_log2FC > 0.25)
  ekegg <- enrichKEGG(gene = data_up$ENTREZID,
                      organism = "mmu",
                      keyType = "ncbi-geneid",
                      minGSSize = 10,
                      maxGSSize = 500,
                      pAdjustMethod = "BH",
                      pvalueCutoff = 1)
  ekegg <- as.data.frame(ekegg)
  ekegg$geneID <- as.character(ekegg$geneID)
  gene_name = c()
  for(i in ekegg$geneID)
  {
    gene_list = c()
    for (j in as.list(strsplit(i,'/'))[[1]])
    {
      gene_list <- c(gene_list,paste(data[data$ENTREZID==j,]$gene,collapse='/'))
    }
    gene_name <- c(gene_name,paste(gene_list,collapse='/'))
  }
  ekegg <- subset(ekegg, select = -c(geneID))
  ekegg['gene_name'] <- gene_name
  write.csv(ekegg,paste0(out_dir,'file/KEGG/',file,'_KEGG_Up.csv'),row.names = FALSE)
  # down-regulated
  data_down <- data %>%
    filter(p_val_adj < 0.05, avg_log2FC < -0.25)
  ekegg <- enrichKEGG(gene = data_down$ENTREZID,
                      organism = "mmu",
                      keyType = "ncbi-geneid",
                      minGSSize = 10,
                      maxGSSize = 500,
                      pAdjustMethod = "BH",
                      pvalueCutoff = 1)
  ekegg <- as.data.frame(ekegg)
  ekegg$geneID <- as.character(ekegg$geneID)
  gene_name = c()
  for(i in ekegg$geneID)
  {
    gene_list = c()
    for (j in as.list(strsplit(i,'/'))[[1]])
    {
      gene_list <- c(gene_list,paste(data[data$ENTREZID==j,]$gene,collapse='/'))
    }
    gene_name <- c(gene_name,paste(gene_list,collapse='/'))
  }
  ekegg <- subset(ekegg, select = -c(geneID))
  ekegg['gene_name'] <- gene_name
  write.csv(ekegg,paste0(out_dir,'file/KEGG/',file,'_KEGG_Down.csv'),row.names = FALSE)
}

#---- Run ORA ------------------------------------------------------------------
data <- c('aggr_AAP_vs_YAP','aggr_0_AAP_vs_YAP',
          'aggr_1_AAP_vs_YAP','aggr_2_AAP_vs_YAP',
          'aggr_3_AAP_vs_YAP','aggr_4_AAP_vs_YAP',
          'aggr_5_AAP_vs_YAP','aggr_6_AAP_vs_YAP')
# GO
lapply(data, GO)
# KEGG
lapply(data, KEGG)

#---- Barplot ------------------------------------------------------------------
