### cell communication analysis between caf and tumor cells
.libPaths("/hsfscqjf1/ST_CQ/P23Z28400N0255/zhangqian/software/miniconda3/envs/seuratv4/lib/R/library")
library(nichenetr)
library(tidyverse)
library(circlize)
library(Seurat)
library(Matrix)
library(Cairo)
library(MuDataSeurat)
options(bitmapType = "cairo")
library(tidyverse)
library(Seurat)
library(dplyr)
library(MuDataSeurat)
library(reticulate)
library(Matrix)
library(reshape2)
library(ComplexHeatmap)
library(cowplot)

seurat_obj <- ReadH5AD('/hsfscqjf1/ST_CQ/P23Z28400N0255/fanfengpu/3.nichnet/dis100um/Type2_Imm_100um_big30cell.h5ad')


setwd('/hsfscqjf1/ST_CQ/P23Z28400N0255/fanfengpu/3.nichnet/dis100um')

library(SeuratDisk)
# Convert("/hsfscqjf1/ST_CQ/P23Z28400N0255/yangzongzheng/BLCA/public_data/umap/new/figures/Endothelial/subtype/all_MIBC_subset.h5ad", "/hsfscqjf1/ST_CQ/P23Z28400N0255/yangzongzheng/BLCA/public_data/umap/new/figures/Endothelial/subtype/all_MIBC_subset.h5seurat")
# obj <- LoadH5Seurat("/hsfscqjf1/ST_CQ/P23Z28400N0255/yangzongzheng/BLCA/public_data/umap/new/figures/Endothelial/subtype/all_MIBC_subset.h5seurat", assays = "RNA", meta.data = FALSE, misc = FALSE)
# metadata <- read.csv("/hsfscqjf1/ST_CQ/P23Z28400N0255/yangzongzheng/BLCA/public_data/umap/new/figures/Endothelial/subtype/all_MIBC_subset.csv", row.names = 1)
# obj <- AddMetaData(object = obj, metadata = metadata)
seurat_obj <- CreateSeuratObject(counts = seurat_obj@assays$RNA@counts, meta.data = seurat_obj@meta.data)

sender_celltypes = c('cDC2_CLEC10A',  'cDC1_CLEC9A',  'Classical_PlasmaB',  'Naive_B_IGHD', 'PlasmaB_SDC1', 'MemoryB_MS4A1', 'cDC3_LAMP3','Macro_LILRB1', 'Macro_APOE','Macro_LYVE1',  'Naive T cell', 'Treg', 'CTL')
receiver = c('Tum_type2', 'Tum_type3', 'Tum_type4')

ligand_target_matrix <- readRDS("/hsfscqjf1/ST_CQ/P23Z28400N0255/yangzongzheng/BLCA/public_data/umap/new/figures/Endothelial/subtype/nichenet/ligand_target_matrix_nsga2r_final.rds")
lr_network = readRDS("/hsfscqjf1/ST_CQ/P23Z28400N0255/yangzongzheng/BLCA/public_data/umap/new/figures/Endothelial/subtype/nichenet/lr_network_human_21122021.rds")
weighted_networks <- readRDS("/hsfscqjf1/ST_CQ/P23Z28400N0255/yangzongzheng/BLCA/public_data/umap/new/figures/Endothelial/subtype/nichenet/weighted_networks_nsga2r_final.rds")
weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))


# seurat_obj<-obj
Idents(seurat_obj) <- seurat_obj$merged_cluster

list_expressed_genes_receiver = receiver %>% unique() %>% lapply(get_expressed_genes, seurat_obj, 0.3) # lapply to get the expressed genes of every receiver cell type separately here
expressed_genes_receiver = list_expressed_genes_receiver %>% unlist() %>% unique()
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]


list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, seurat_obj, 0.3) # lapply to get the expressed genes of every sender cell type separately here
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()

DE_table_receiver = FindMarkers(object = seurat_obj, ident.1 = receiver, ident.2 = sender_celltypes, min.pct = 0.10) %>% rownames_to_column("gene")
geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.5) %>% pull(gene)
geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

ligands = lr_network %>% pull(from) %>% unique()
receptors = lr_network %>% pull(to) %>% unique()

expressed_ligands = intersect(ligands,expressed_genes_sender)
expressed_receptors = intersect(receptors,expressed_genes_receiver)

potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()

ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)

ligand_activities = ligand_activities %>% arrange(-aupr_corrected) %>% mutate(rank = rank(desc(aupr_corrected)))
best_upstream_ligands = ligand_activities %>% top_n(100, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand) %>% unique()

#select ligands with high expression in sender cells
# x
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 50) %>% bind_rows() %>% drop_na()
active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.6)

order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names()
colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names()

### select target genes
interactions <- active_ligand_target_links[order_targets, ]
interaction_strength <- rowSums(interactions, na.rm = TRUE)
sorted_genes <- order(interaction_strength, decreasing = TRUE)
top_genes <- order_targets[sorted_genes[1:20]]

order_targets1 <- c('KRT13','KRT19','ACKR3','BMPR2','ITGA5','LMCD1','SERPINE1','SPRY1','VEGFA','CCL2','CCN2','CDH1','CXCL8','ID1','CXCR4','ICAM1','IL1B','IL32','ANXA5','PECAM1','ITGA6','APLP2','HSPE1','CRCT1','KDR','SOCS3','CD44')
# order_targets1 <- unique(c(order_targets1,top_genes))
order_targets <- order_targets1[order_targets1 %in% order_targets]
order_targets <- order_targets[!order_targets %in% c('COL4A1','VEGFA','VIM','APP','MYO1B','COL15A1','SOX5','CSGALNACT1','SLCO2A1')]


order_ligands1 <- c('TGFB1','ADAM10','ADAM17','ENG','HMGB1','NAMPT','IGFBP3','JAG1','CMTM8')
order_ligands <- order_ligands1[order_ligands1 %in% order_ligands]

vis_ligand_target <- active_ligand_target_links[order_targets,order_ligands] %>% t()
vis_ligand_target <- vis_ligand_target[,order_targets]
p_ligand_target_network <- vis_ligand_target %>%
  make_heatmap_ggplot("Prioritized ligands", "Predicted target genes", color = "purple", legend_position = "top", x_axis_position = "top", legend_title = "Regulatory potential") +
  theme(
    axis.title.x = element_text(size = 12, face = "italic"),
    axis.title.y = element_text(size = 12, face = "italic"),
    axis.text.x = element_text(size = 10, face = "italic"), 
    axis.text.y = element_text(size = 10),  
    legend.title = element_text(size = 12),  
    legend.text = element_text(size = 10)  
  ) +
  scale_fill_gradient2(low = "whitesmoke", high = "purple", breaks = c(0, 0.15, 0.3))
#Diagram of Interaction Strength
pdf("p_ligand_target_network2_interested_ECtoTumor.pdf",height=3,width=4.5)
print(p_ligand_target_network)
dev.off()


###ligand expression in sender cells
ligands <- c('TGFB1','ADAM10','ADAM17','ENG','HMGB1','NAMPT','IGFBP3','JAG1','CMTM8')
# ligands <- rev(ligands)
seurat_obj<-obj
seurat_obj <- subset(seurat_obj, CellTypeL1 %in% c('SOX5+ EC'))
seurat_obj <- NormalizeData(seurat_obj)
pdf("order_targets_interested_ligands_ECtoTumor.pdf", height=2.8, width=2.2)
p <- DotPlot(seurat_obj, features = ligands, cols = "RdYlBu", group.by = "CellTypeL1") +
  coord_flip() + 
  scale_y_discrete(position = "right") +
  RotatedAxis() +
  theme(legend.box = "vertical",
        legend.key.size = unit(0.2, "cm"), 
        legend.title = element_text(size = 8), 
        legend.text = element_text(size = 6), 
        axis.title.y = element_blank(),
        axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 8.5,angle = 90,hjust = 0), 
        axis.text.y = element_text(size = 8.5))
print(p)
dev.off()

###target gene expression in receiver cells
order_targets <- c('KRT13','KRT19','ACKR3','BMPR2','ITGA5','LMCD1','SERPINE1','SPRY1','VEGFA','CCL2','CCN2','CDH1','CXCL8','ID1','CXCR4','ICAM1','IL1B','IL32','ANXA5','PECAM1','ITGA6','APLP2','HSPE1',
'CRCT1','KDR','SOCS3','CD44')

seurat_obj<-obj
seurat_obj <- subset(seurat_obj, CellTypeL1 %in% c('Basal_tumor_cell'))
seurat_obj <- NormalizeData(seurat_obj)
pdf("order_targets_interested_targets_ECtoTumor.pdf", height=2, width=5.5)
print(DotPlot(seurat_obj, features = order_targets, cols = "RdYlBu", group.by = "CellTypeL1") +
      RotatedAxis() +
      theme(legend.box = "vertical",
            legend.key.size = unit(0.06, "cm"),
            legend.title = element_text(size = 6),
            legend.text = element_text(size = 5),
            axis.title.y = element_blank(),
            axis.text.x = element_text(size = 8),
            axis.text.y = element_text(size = 8),
            axis.title.x = element_text(size=10)) + 
      xlab("Target genes in tumor") 
)
dev.off()


#### ligand-receptor-target network
library(DiagrammeR)
library(DiagrammeRsvg)

weighted_networks = readRDS("weighted_networks_nsga2r_final.rds")
ligand_tf_matrix = readRDS("ligand_target_matrix_nsga2r_final.rds")
lr_network = readRDS("lr_network_human_allInfo_30112033.rds")
sig_network = readRDS("signaling_network_human_21122021.rds")
gr_network = readRDS("gr_network_human_21122021.rds")


ligands_all = c("POSTN")
targets_all = c("FOS")
# extract network
active_signaling_network = get_ligand_signaling_path(
  ligand_tf_matrix = ligand_tf_matrix, 
  ligands_all = ligands_all, 
  targets_all = targets_all,
  weighted_networks = weighted_networks)
lapply(active_signaling_network, head)
## normalize weights
active_signaling_network_min_max = active_signaling_network
active_signaling_network_min_max$sig = active_signaling_network_min_max$sig %>% 
  mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)
active_signaling_network_min_max$gr = active_signaling_network_min_max$gr %>% 
  mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)

graph_min_max = diagrammer_format_signaling_graph(
  signaling_graph_list = active_signaling_network_min_max, 
  ligands_all = ligands_all, 
  targets_all = targets_all, 
  sig_color = "indianred", 
  gr_color = "steelblue")
DiagrammeR::render_graph(graph_min_max, 
                         output="graph", # graph, visNetwork
                         layout = "kk" # nicely, circle, tree, kk, and fr
)
