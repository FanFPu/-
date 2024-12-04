library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(Seurat)
options(stringsAsFactors = FALSE, future.globals.maxSize=20000*1024^2)

od <- '/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/Paper/Spatial/9.cellchat/diff_dis_um/merge_analysis_100um/'
setwd(od)

x <- load('/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/Paper/Spatial/9.cellchat/diff_dis_um/dis100um/cellchat_interaction.RData')
cellchat_high_TIL <- get(x)
y <- load('/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/Paper/Spatial/9.cellchat/diff_dis_um/dis200um/cellchat_interaction.RData')
cellchat_low_TIL <- get(y)

cellchat <- mergeCellChat(list(cellchat_high_TIL, cellchat_low_TIL), add.names = c("100um", "200um"), cell.prefix = TRUE)

pdf(paste0(od, "compare_number.pdf"),width=3,height=3)
compareInteractions(cellchat, show.legend = F, group = c(1,2))
compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
dev.off()

pdf(paste0(od, "compare_number_bycell.pdf"),width=6,height=4)
netVisual_heatmap(cellchat)
netVisual_heatmap(cellchat, measure = "weight")
dev.off()

pdf(paste0(od, "rankNET.pdf"),width=4,height=9)
rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
dev.off()

library(ComplexHeatmap)
cellchat_high_TIL <- netAnalysis_computeCentrality(cellchat_high_TIL, slot.name = "netP")
cellchat_low_TIL <- netAnalysis_computeCentrality(cellchat_low_TIL, slot.name = "netP")
object.list <- list(um100 = cellchat_high_TIL, um200 = cellchat_low_TIL)
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 21)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 21)
pdf(paste0(od, "compare_outgoing_signal_bycell.pdf"),width=8,height=13)
draw(ht1 + ht2, ht_gap = unit(1.5, "cm"))
dev.off()
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 21, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 21, color.heatmap = "GnBu")
pdf(paste0(od, "compare_incoming_signal_bycell.pdf"),width=8,height=13)
draw(ht1 + ht2, ht_gap = unit(1.5, "cm"))
dev.off()
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 21, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 21, color.heatmap = "OrRd")
pdf(paste0(od, "compare_all_signal_bycell.pdf"),width=8,height=13)
draw(ht1 + ht2, ht_gap = unit(1.5, "cm"))
dev.off()

##### Compare L-R pairs
#pdf(paste0(od, "L-R_InvEpi_Immune_compare.pdf"),width=6,height=12)
#netVisual_bubble(cellchat, sources.use = c('Epithelial_Invasive_cancer','myCAF'), targets.use = c('B_cell','Macrophage','Plasmocyte','T_cell'),  comparison = c(1, 2), angle.x = 45)
#dev.off()

pdf(paste0(od, "IncreaseorDecrease_Plasma_tumor_signal.pdf"),width=15,height=24)
gg1 <- netVisual_bubble(cellchat, sources.use = c('Classical_PlasmaB','PlasmaB_SDC1'), targets.use = c('Tum_type2', 'Tum_type3', 'Tum_type4'),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in low_TIL samples", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = c('Classical_PlasmaB','PlasmaB_SDC1'), targets.use = c('Tum_type2', 'Tum_type3', 'Tum_type4'),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in low_TIL samples", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2
dev.off()


#这个地方就是在看200um的东西，subsetCommunication函数只能测试到上调的gene
# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "100um"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
#> Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat, net = net, datasets = "100um",ligand.logFC = 0.3, receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat, net = net, datasets = "200um",ligand.logFC = -0.1, receptor.logFC = -0.1)
gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)
pairLR.use.up = net.up[, "interaction_name", drop = F]
print(pairLR.use.up)
pairLR.use.down = net.down[, "interaction_name", drop = F]

# pairLR.use.up1 <- pairLR.use.up %>% 
#   filter(!grepl("^COL", interaction_name))
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = c('Classical_PlasmaB','PlasmaB_SDC1'), targets.use = c('Tum_type2', 'Tum_type3', 'Tum_type4'),comparison = c(1, 2),  angle.x = 45, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))+
      theme(text = element_text(size = 12)) # 设置字体大小为12

#> Comparing communications on a merged object
#gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = c('IBC cells','Fibroblasts'), targets.use = c('Fibroblasts'), 
                        #comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
pdf(paste0(od, "Up-regulated_Plasma_sender_signal.pdf"),width=7,height=12)
gg1
dev.off()
# 

gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = c('Classical_PlasmaB','PlasmaB_SDC1'), targets.use = c('Tum_type2', 'Tum_type3', 'Tum_type4'),  comparison = c(1, 2),  angle.x = 45, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[1]))
pdf(paste0(od, "Down-regulated_Plasma_sender_signal.pdf"),width=7,height=30)
gg2
dev.off()



gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = c('Macro_LILRB1', 'Macro_APOE','Macro_LYVE1'), targets.use = c('Tum_type2', 'Tum_type3', 'Tum_type4'),comparison = c(1, 2),  angle.x = 45, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[1]))+
      theme(text = element_text(size = 12)) # 设置字体大小为12

#> Comparing communications on a merged object
#gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = c('IBC cells','Fibroblasts'), targets.use = c('Fibroblasts'), 
                        #comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
pdf(paste0(od, "Up-regulated_Mac_sender_signal.pdf"),width=7,height=12)
gg1
dev.off()
# 

gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = c('Macro_LILRB1', 'Macro_APOE','Macro_LYVE1'), targets.use = c('Tum_type2', 'Tum_type3', 'Tum_type4'),  comparison = c(1, 2),  angle.x = 45, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
pdf(paste0(od, "Down-regulated_Mac_sender_signal.pdf"),width=7,height=30)
gg2
dev.off()


# gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = c('Classical_PlasmaB','PlasmaB_SDC1'), targets.use = c('Tum_type2', 'Tum_type3', 'Tum_type4'), 
#                         comparison = c(1, 2),  angle.x = 45, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[1]))
# pdf(paste0(od, "Up-regulated_IBC_sender_signal.pdf"),width=5,height=28)
# gg2
# dev.off()

# gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = c('Classical_PlasmaB','PlasmaB_SDC1'), targets.use = c('Tum_type2', 'Tum_type3', 'Tum_type4'), 
#                         comparison = c(1, 2),  angle.x = 45, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[1]))
# pdf(paste0(od, "Down-regulated_IBC_sender_signal.pdf"),width=5,height=34)
# gg2
# dev.off()

##### plot gene expression
cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("high_TIL", "low_TIL")) # set factor level
pdf(paste0(od, "GeneExpression.pdf"),width=13,height=100)
plotGeneExpression(cellchat, signaling = "Other", split.by = "datasets", colors.ggplot = T)
dev.off()