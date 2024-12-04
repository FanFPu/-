library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(Seurat)
options(stringsAsFactors = FALSE, future.globals.maxSize=20000*1024^2)

od <- '/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/Paper/Spatial/9.cellchat/diff_dis_um/dis200um/'
setwd(od)
load("/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/Paper/Spatial/9.cellchat/diff_dis_um/dis200um/cellchat_interaction.RData")

## visualize
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1, 2), xpd = TRUE)
pdf(paste0(od, "/circlePlot.pdf"))
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge = F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge = F, title.name = "Interaction weights/strength")
dev.off()

# Due to the complicated cell-cell communication network, we can examine the signaling sent from each cell group. Here we also control the parameter edge.weight.max so that we can compare edge weights between differet network.
mat <- cellchat@net$weight
par(mfrow = c(3, 4), xpd = TRUE)
pdf(paste0(od, "/circlePlot_subset.pdf"))
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()

pdf(paste0(od, "/Fibroblast-Tum_type2.pdf"),height = 10, width = 3.5)
netVisual_bubble(cellchat, sources.use = c('Fibroblast'), #发出的 
                 targets.use = c('Tum_type2'), signaling = c("ADGRE5",'THBS','TGFb','SPP1','SEMA5','SEMA3','PERIOSTIN','ncWNT','MPZ','CDH1','COLLAGEN','CXCL','IGF','FN1','MK'), remove.isolate = FALSE) #接收的
dev.off()
# bubble_plot <- bubble_plot + theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) #把横轴左边变为横着的



# pdf(paste0(od, "/fib&DCIS&Peri-MEC.pdf"),height = 8, width = 5)
# netVisual_bubble(cellchat, sources.use = c('Fibroblasts','DCIS cells','Pericytes'), #发出的 
#                  targets.use = c('MEC'), signaling = c("SEMA5","BMP",'THBS','FN1','COLLAGEN',"EGF","SEMA3",'CDH1','ADGRE5','ANGPTL','CSPG4'), remove.isolate = FALSE) #接收的
# dev.off()
# pdf(paste0(od, "/MEC-DCIS.pdf"),height = 5, width = 3)
# netVisual_bubble(cellchat, sources.use = c('MEC'), #发出的 
#                  targets.use = c('DCIS cells'), signaling = c("NECTIN","ncWNT",'AGRN','NOTCH'), remove.isolate = FALSE) #接收的
# dev.off()
# pdf(paste0(od, "/DCIS-MEC.pdf"),height = 5, width = 3)
# netVisual_bubble(cellchat, sources.use = c('DCIS cells'), #发出的 
#                  targets.use = c('MEC'), signaling = c("EGF","SEMA3",'CDH1'), remove.isolate = FALSE) #接收的
# dev.off()

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")#根据每个信号通路，哪些细胞类型在该信号通路互作的
pathways <- cellchat@netP$pathways
dir.create(paste0(od, "/pathway_netVisual_aggregate"), recursive = TRUE)
for (pathways.show in pathways) {
  # Hierarchy plot
  # Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells
  vertex.receiver <- seq(1, 4) # a numeric vector.
  pdf(paste0(od, "/pathway_netVisual_aggregate/", pathways.show, ".pdf"))
  netVisual_aggregate(cellchat, signaling = pathways.show, vertex.receiver = vertex.receiver, layout = "circle", vertex.label.cex = 0.8)
  netVisual_chord_gene(cellchat, signaling = pathways.show, slot.name = "netP", legend.pos.x = 15, legend.pos.y = 0, title.name=paste0(pathways.show," signaling network"))
  netVisual_chord_cell(cellchat, signaling = pathways.show, slot.name = "netP", legend.pos.x = 15, title.name=paste0(pathways.show," signaling network"))
  netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 12, height = 2.5, font.size = 10)
  dev.off()
}


#pdf(paste0(od, "/pathway_netVisual_aggregate/signaling/heatmap.pdf"))
#netVisual_heatmap(cellchat, signaling = 'COLLAGEN', color.heatmap = "Reds")
#netVisual_heatmap(cellchat, signaling = 'FN1', color.heatmap = "Reds")
#netVisual_heatmap(cellchat, signaling = 'IGF', color.heatmap = "Reds")
#dev.off()

#pdf(paste0(od, "/pathway_netVisual_aggregate/signaling/netContribution.pdf"))
#netAnalysis_contribution(cellchat, signaling = 'APRIL')
#dev.off()
# dir.create(paste0(od, "/signaling"), recursive = TRUE)
# pdf(paste0(od, "/signaling/expression.pdf"))
# plotGeneExpression(cellchat, signaling = "NECTIN", split.by = "labels", colors.ggplot = T)
# plotGeneExpression(cellchat, signaling = "ncWNT", split.by = "labels", colors.ggplot = T)
# dev.off()

#par(mfrow = c(1, 2), xpd = TRUE)
#pdf(paste0(od, "/SignalingRole_heatmap.pdf"))
#netAnalysis_SignalingRole_heatmap(cellchat, pattern = "outgoing")
#netAnalysis_SignalingRole_heatmap(cellchat, pattern = "incoming")
#dev.off()

# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# aim_sig <-c("COLLAGEN","LAMININ","FN1","MK","THBS","SPP1","NOTCH","JAM","VEGF","APP","ANGPTL",“,"ICAM","ncWNT","TENASCIN",”,"CXCL","NCAM","GRN","MHC-I","SEMA4","GALECTIN","COMPLEMENT","PERIOSTIN","ITGB2","CDH","PECAM1","EPHB","BMP","HSPG","EPHA","GAS","ACTIVIN","VCAM","CADM","VISFATIN","DESMOSOME","FGF","TGFb","NRG","CHEMERIN","CCL","GDF","MPZ","CDH1","CD46","PARs","CSF","RESISTIN","AGRN","NECTIN","SEMA5","EDN","SEMA6","UGRP1","NEGR","CD45","NPNT","CNTN","PTN","ANGPT","ADGRE5","PROS","CD39","IL16","CSPG4","EGF","CEACAM","CALCR","SEMA7","PSAP","TWEAK","PTPRM","KIT","PDGF","NRXN","VWF","SELL")
ht1 <- netAnalysis_signalingRole_heatmap(cellchat,pattern = "outgoing",width = 10,height = 17.5)
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, signaling =rownames(ht1@matrix)[1:59],pattern = "outgoing",width = 6,height = 17.5)

ht2 <- netAnalysis_signalingRole_heatmap(cellchat, signaling =rownames(ht1@matrix)[1:59], pattern = "incoming",width = 6,height = 17.5)
pdf("netAnalysis_signalingRole_heatmap.pdf",height=17.5,width=12)
ht1 + ht2
dev.off()