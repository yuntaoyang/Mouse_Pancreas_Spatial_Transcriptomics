#---- Load Packages ------------------------------------------------------------
library(ggplot2)
library(readr)
library(dplyr)

#---- Set Directory ------------------------------------------------------------
out_dir <- ''

#---- DEG: 30 ------------------------------------------------------------------
degs <- read.csv(paste0(out_dir,'Selected_DEG_30.csv'))
load_data <- function(file_name,cluster)
{
  df <- read.csv(paste0(out_dir,file_name)) %>%
    select(X,avg_log2FC,p_val_adj) %>%
    filter(X %in% degs$genes) %>%
    rename('genes' = 'X', 'FDR' = 'p_val_adj', 'Log2FC' = 'avg_log2FC') %>%
    arrange(factor(genes, levels = degs$genes))
  df$type <- rep(cluster,30)
  return(df)
}
aggr_all <- load_data('aggr_AAP_vs_YAP.csv','All Clusters')
aggr_0 <- load_data('aggr_0_AAP_vs_YAP.csv','Acinar_C0')
aggr_1 <- load_data('aggr_1_AAP_vs_YAP.csv','Acinar_C1')
aggr_2 <- load_data('aggr_2_AAP_vs_YAP.csv','Acinar_C2')
aggr_3 <- load_data('aggr_3_AAP_vs_YAP.csv','Stroma_C3')
aggr_4 <- load_data('aggr_4_AAP_vs_YAP.csv','Acinar_C4')
aggr_5 <- load_data('aggr_5_AAP_vs_YAP.csv','Stroma_C5')
aggr_6 <- load_data('aggr_6_AAP_vs_YAP.csv','Islet_C6')

data <- rbind(aggr_all,aggr_0,aggr_1,aggr_2,aggr_3,aggr_4,aggr_5,aggr_6) %>%
  filter(FDR < 0.05, Log2FC > 0.25 | Log2FC < -0.25)
data$type <- factor(data$type, levels = c('All Clusters','Acinar_C0','Acinar_C1','Acinar_C2',
                                          'Acinar_C4','Stroma_C3','Stroma_C5','Islet_C6'))
data$genes <- factor(data$genes, levels = rev(degs$genes))
write.csv(data, paste0(out_dir,'dotplot_30.csv'), row.names = FALSE)

p <- ggplot(data, aes(x = type, y = genes, color = Log2FC, size = -log10(FDR))) + 
  geom_point() +
  theme_bw() +
  labs(x = expression("Aging"^"AAP" ~ "vs." ~ "Young"^"AAP"), y = 'Differentially Expressed Genes') +
  scale_size_continuous(breaks=c(5,10,50,100)) +
  scale_color_gradient2(low = "blue4", mid = 'white', high = "red2", midpoint = 0) +
  theme(legend.text = element_text(size = 25),
        legend.title = element_text(size = 25),
        legend.key.size = unit(1.5, 'cm'),
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x = element_text(angle = 60, vjust = 0.5, size = 20),
        axis.text.y = element_text(size = 20))
ggsave(paste0(out_dir,'dotplot_30.png'),width=10,height=12)

#---- DEG: 31 ------------------------------------------------------------------
degs <- read.csv(paste0(out_dir,'Selected_DEG_31.csv'))
load_data <- function(file_name,cluster)
{
  df <- read.csv(paste0(out_dir,file_name)) %>%
    select(X,avg_log2FC,p_val_adj) %>%
    filter(X %in% degs$genes) %>%
    rename('genes' = 'X', 'FDR' = 'p_val_adj', 'Log2FC' = 'avg_log2FC') %>%
    arrange(factor(genes, levels = degs$genes))
  df$type <- rep(cluster,31)
  return(df)
}
aggr_all <- load_data('aggr_AAP_vs_YAP.csv','All Clusters')
aggr_0 <- load_data('aggr_0_AAP_vs_YAP.csv','Acinar_C0')
aggr_1 <- load_data('aggr_1_AAP_vs_YAP.csv','Acinar_C1')
aggr_2 <- load_data('aggr_2_AAP_vs_YAP.csv','Acinar_C2')
aggr_3 <- load_data('aggr_3_AAP_vs_YAP.csv','Stroma_C3')
aggr_4 <- load_data('aggr_4_AAP_vs_YAP.csv','Acinar_C4')
aggr_5 <- load_data('aggr_5_AAP_vs_YAP.csv','Stroma_C5')
aggr_6 <- load_data('aggr_6_AAP_vs_YAP.csv','Islet_C6')

data <- rbind(aggr_all,aggr_0,aggr_1,aggr_2,aggr_3,aggr_4,aggr_5,aggr_6) %>%
  filter(FDR < 0.05, Log2FC > 0.25 | Log2FC < -0.25)
data$type <- factor(data$type, levels = c('All Clusters','Acinar_C0','Acinar_C1','Acinar_C2',
                                          'Acinar_C4','Stroma_C3','Stroma_C5','Islet_C6'))
data$genes <- factor(data$genes, levels = rev(degs$genes))
write.csv(data, paste0(out_dir,'dotplot_31.csv'), row.names = FALSE)

p <- ggplot(data, aes(x = type, y = genes, color = Log2FC, size = -log10(FDR))) + 
  geom_point() +
  theme_bw() +
  labs(x = expression("Aging"^"AAP" ~ "vs." ~ "Young"^"AAP"), y = 'Differentially Expressed Genes') +
  scale_size_continuous(breaks=c(5,10,50,100)) +
  scale_color_gradient2(low = "blue4", mid = 'white', high = "red2", midpoint = 0) +
  theme(legend.text = element_text(size = 25),
        legend.title = element_text(size = 25),
        legend.key.size = unit(1.5, 'cm'),
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x = element_text(angle = 60, vjust = 0.5, size = 20),
        axis.text.y = element_text(size = 20))
ggsave(paste0(out_dir,'dotplot_31.png'),width=10,height=12)

#---- GO: top 20 up-regulated pathways -----------------------------------------
pathways <- read.csv(paste0(out_dir,'Top_Pathways.csv'))
colnames(pathways)[1] <- "Description"
pathways$GO <- paste(seq_along(pathways$Description), pathways$Description, sep = ". ") 
load_data <- function(file_name,cluster)
{
  df <- read.csv(paste0(out_dir,'ORA/file/GO/',file_name)) %>%
    select(Description,qvalue,Count) %>%
    filter(Description %in% pathways$Description) %>%
    left_join(pathways, by = 'Description') %>%
    rename('FDR' = 'qvalue') %>%
    arrange(factor(GO, levels = pathways$GO))
  df$type <- rep(cluster,dim(df)[1])
  return(df)
}
aggr_all <- load_data('aggr_AAP_vs_YAP_GO_Up.csv','All Clusters')
aggr_0 <- load_data('aggr_0_AAP_vs_YAP_GO_Up.csv','Acinar_C0')
aggr_1 <- load_data('aggr_1_AAP_vs_YAP_GO_Up.csv','Acinar_C1')
aggr_2 <- load_data('aggr_2_AAP_vs_YAP_GO_Up.csv','Acinar_C2')
aggr_3 <- load_data('aggr_3_AAP_vs_YAP_GO_Up.csv','Stroma_C3')
aggr_4 <- load_data('aggr_4_AAP_vs_YAP_GO_Up.csv','Acinar_C4')
aggr_5 <- load_data('aggr_5_AAP_vs_YAP_GO_Up.csv','Stroma_C5')
aggr_6 <- load_data('aggr_6_AAP_vs_YAP_GO_Up.csv','Islet_C6')

data <- rbind(aggr_all,aggr_0,aggr_1,aggr_2,aggr_3,aggr_4,aggr_5,aggr_6) %>%
  filter(FDR < 0.05)
data$type <- factor(data$type, levels = c('All Clusters','Acinar_C0','Acinar_C1','Acinar_C2',
                                          'Acinar_C4','Stroma_C3','Stroma_C5','Islet_C6'))
data$GO <- factor(data$GO, levels = rev(pathways$GO))
colnames(data)[3] <- "size"
write.csv(data, paste0(out_dir,'ORA/figure/','dotplot_top20_pathways.csv'), row.names = FALSE)

p <- ggplot(data, aes(x = type, y = GO, color = -log10(FDR), size = size)) + 
  geom_point() +
  theme_bw() +
  labs(x = expression("Aging"^"AAP" ~ "vs." ~ "Young"^"AAP"), y = 'GO') +
  scale_size_continuous(breaks=c(10,30,50)) +
  scale_color_gradient(low = "blue", high = "red", breaks = c(5,10,15,20,25)) +
  theme(legend.text = element_text(size = 25),
        legend.title = element_text(size = 25),
        legend.key.size = unit(1.5, 'cm'),
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x = element_text(angle = 60, vjust = 0.5, size = 20),
        axis.text.y = element_text(size = 20))
ggsave(paste0(out_dir,'ORA/figure/','dotplot_top20_pathways.png'),width=16,height=12)

