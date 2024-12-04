library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(Seurat)
# library(SeuratData)
library(magrittr)
library(MuDataSeurat)
options(stringsAsFactors = FALSE, future.globals.maxSize=20000*1024^2)


library(doParallel)
future::plan("multisession", workers = 6)  

#成纤维、T、内皮、周细胞、巨噬、DC、B
od <- '/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/Paper/Spatial/9.cellchat/diff_dis_um/dis100um/'
setwd(od)

sc <- MuDataSeurat::ReadH5AD('/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/Paper/Spatial/9.cellchat/diff_dis_um/dis100um/Type2_Imm_100um_big30cell.h5ad')
sc$annotated_cluster <- sc$merged_cluster
str(sc)
sc <- NormalizeData(sc, normalization.method = "LogNormalize", scale.factor = 10000)
#sc <- sc[, ! sc@meta.data$annotated_cluster %in% c('Atypical hyperplasia epithelial cells')]
# sc <- sc[, ! sc@meta.data$annotated_cluster %in% c('low_quality')]

#load('/jdfsbjcas1/ST_BJ/P21H28400N0232/wangdi/xieguixiang/BRCAsc_Bcell_plasmocyte/cell_communication/sc/B_other.RData')
#在这里注意细胞类型
sc$annotated_cluster < as.factor(sc$annotated_cluster)
table(sc$annotated_cluster)
data.input  <- sc@assays$RNA@data
identity = data.frame(group=sc$annotated_cluster, row.names = names(sc$annotated_cluster)) # create a dataframe consisting of the cell labels
head(identity)
cellchat <- createCellChat(data.input)
cellchat <- addMeta(cellchat, meta = identity, meta.name = "labels")
head(cellchat@meta)
str(cellchat)

cellchat <- setIdent(cellchat, ident.use = "labels")
groupSize <- as.numeric(table(cellchat@idents))
table(cellchat@idents)
save(cellchat, file="cellchat.RData")

#CellChatDB <- CellChatDB.human
CellChatDB <- readRDS('/jdfsbjcas1/ST_BJ/P21H28400N0232/gongchanghao/database/cellchat/human/CellChatDB.human_nichenet.rds')
CellChatDB.use <- CellChatDB
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling for cell-cell communication analysis
cellchat@DB <- CellChatDB.use # set the used database in the object
unique(CellChatDB$interaction$annotation)
save(cellchat, file = "cellchat_allDB.RData")

## pre-processing
cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
# future::plan("multiprocess", workers = 10) # do parallel  
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
save(cellchat, file = "cellchat_preprocess.RData")

## interaction inference
cellchat <- computeCommunProb(cellchat, population.size = TRUE)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat@netP$pathways
head(cellchat@LR$LRsig)
save(cellchat, file = "cellchat_interaction.RData")
