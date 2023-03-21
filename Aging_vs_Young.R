#---- Load Packages ------------------------------------------------------------
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(glmGamPoi)

#---- Set Directory ------------------------------------------------------------
root_dir <- ""
out_dir <- ''
h5_mat_name <- './filtered_feature_bc_matrix.h5'

#---- Load data into Seurat ----------------------------------------------------
load_data <- function(sample_id,slice)
{
  setwd(paste0(root_dir,sample_id))
  object <- Seurat::Load10X_Spatial(
    data.dir = paste0(root_dir,sample_id), 
    filename = h5_mat_name,
    assay = 'Spatial',
    slice = slice, 
    filter.matrix = TRUE, 
    to.upper = FALSE)
  object <- SCTransform(object, assay = 'Spatial', verbose = FALSE, method="glmGamPoi", min_cells=1)
  object$orig.ident <- sample_id
  return(object)
}
YAP_3859 <- load_data('3859-YAP','YAP3859')
AAP_3956 <- load_data('3956-AAP','AAP3956')
YAP_3858 <- load_data('3858-YAP','YAP3858')
AAP_3955 <- load_data('3955-AAP','AAP3955')

#---- Merge samples ------------------------------------------------------------
st.list = list(YAP_3859 = YAP_3859, AAP_3956 = AAP_3956, 
               YAP_3858 = YAP_3858, AAP_3955 = AAP_3955)
options(future.globals.maxSize = 2000 * 1024^2)
st.features = SelectIntegrationFeatures(st.list, nfeatures = 3000, verbose = FALSE)
st.list <- PrepSCTIntegration(object.list = st.list, anchor.features = st.features,
                              verbose = FALSE)
int.anchors <- FindIntegrationAnchors(object.list = st.list, normalization.method = "SCT",
                                      verbose = FALSE, anchor.features = st.features)
sample_merge <- IntegrateData(anchorset = int.anchors, normalization.method = "SCT",
                              verbose = FALSE)

#---- Dimensional reduction and clustering -------------------------------------
cluster <- function(object, resolution)
{
  object <- RunPCA(object, verbose = FALSE)
  object <- FindNeighbors(object, dims = 1:30)
  object <- FindClusters(object, verbose = FALSE, resolution = resolution)
  object <- RunUMAP(object, reduction = "pca", dims = 1:30)
  object <- RunTSNE(object, reduction = "pca", dims = 1:30)
  return(object)
}
sample_merge <- cluster(sample_merge,0.3)

#---- After rename clusters based on condition ---------------------------------
sample_merge$sample_type <- Idents(sample_merge)
sample_merge$sample_type.ident <- paste(Idents(sample_merge), sample_merge$orig.ident, sep = "_")
Idents(sample_merge) <- "sample_type.ident"

#---- Count spots for each cluster ---------------------------------------------
cluster_count <- as.data.frame(table(sample_merge$sample_type.ident))
colnames(cluster_count) <- c('cluster','count')
cluster_count$group_1 <- c(rep('Acinar_C0',4),rep('Acinar_C1',4),rep('Acinar_C2',4),
                           rep('Stroma_C3',4),rep('Acinar_C4',4),rep('Stroma_C5',4),
                           rep('Islet_C6',4),'Acinar_C7')
cluster_count$group_2 <- c(rep(c('Y','Y','A','A'),7),'A')
# clusters in four samples
cluster_all <- aggregate(cluster_count$count, by=list(cluster=cluster_count$group_1), FUN=sum)
colnames(cluster_all) <- c('cluster','count')
write.csv(cluster_all,paste0(out_dir,'cluster_count_all.csv'),row.names = FALSE)
# clusters in young samples
cluster_young <- filter(cluster_count, group_2=='Y')
cluster_young <- aggregate(cluster_young$count, by=list(cluster=cluster_young$group_1), FUN=sum)
colnames(cluster_young) <- c('cluster','count')
write.csv(cluster_young,paste0(out_dir,'cluster_count_young.csv'),row.names = FALSE)
# clusters in aging samples
cluster_aging <- filter(cluster_count, group_2=='A')
cluster_aging <- aggregate(cluster_aging$count, by=list(cluster=cluster_aging$group_1), FUN=sum)
colnames(cluster_aging) <- c('cluster','count')
write.csv(cluster_aging,paste0(out_dir,'cluster_count_aging.csv'),row.names = FALSE)

#---- Identification of Spatially Variable Features ----------------------------
# differential expression between any clusters
de_markers <- function(object,cluster_1,cluster_2,file_name)
{
  markers <- FindMarkers(object, ident.1 = cluster_1, ident.2 = cluster_2,
                         logfc.threshold =0, min.pct = 0)
  write.csv(markers,paste0(out_dir,'aggr_',file_name,'.csv'))
}

# All clusters in Aging vs. All clusters in Young
YAP_clustes <- c('0_3858-YAP','1_3858-YAP','2_3858-YAP','3_3858-YAP','4_3858-YAP','5_3858-YAP','6_3858-YAP',
                 '0_3859-YAP','1_3859-YAP','2_3859-YAP','3_3859-YAP','4_3859-YAP','5_3859-YAP','6_3859-YAP')
AAP_clustes <- c('0_3955-AAP','1_3955-AAP','2_3955-AAP','3_3955-AAP','4_3955-AAP','5_3955-AAP','6_3955-AAP','7_3955-AAP',
                 '0_3956-AAP','1_3956-AAP','2_3956-AAP','3_3956-AAP','4_3956-AAP','5_3956-AAP','6_3956-AAP')
de_markers(sample_merge,AAP_clustes,YAP_clustes,'AAP_vs_YAP')
# Cluster n in Aging (AAP_1 & AAP_2) vs. Cluster n in Young (YAP_1 & YAP_2)
de_markers(sample_merge,c('0_3955-AAP','0_3956-AAP'),c('0_3858-YAP','0_3859-YAP'),'0_AAP_vs_YAP')
de_markers(sample_merge,c('1_3955-AAP','1_3956-AAP'),c('1_3858-YAP','1_3859-YAP'),'1_AAP_vs_YAP')
de_markers(sample_merge,c('2_3955-AAP','2_3956-AAP'),c('2_3858-YAP','2_3859-YAP'),'2_AAP_vs_YAP')
de_markers(sample_merge,c('3_3955-AAP','3_3956-AAP'),c('3_3858-YAP','3_3859-YAP'),'3_AAP_vs_YAP')
de_markers(sample_merge,c('4_3955-AAP','4_3956-AAP'),c('4_3858-YAP','4_3859-YAP'),'4_AAP_vs_YAP')
de_markers(sample_merge,c('5_3955-AAP','5_3956-AAP'),c('5_3858-YAP','5_3859-YAP'),'5_AAP_vs_YAP')
de_markers(sample_merge,c('6_3955-AAP','6_3956-AAP'),c('6_3858-YAP','6_3859-YAP'),'6_AAP_vs_YAP')
