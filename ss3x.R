###导入如下包
library(limma)
library(Seurat)
library(dplyr)
library(magrittr)
library(ggrepel)
library(reticulate)
library(Matrix)
library(cowplot)
library(R.utils)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(SingleR)
library(celldex)
library(immunarch)
library(corrplot)
library(tidyverse)
library(maftools)
library(future)
library(data.table)

plan("multicore")
options(future.globals.maxSize = 1000 * 1024^5)
###DOWNLOAD AND CHECK DATA######################################################

source('/realspace/pe/liyong/sc_standard_pipline/s1.DCAF.r')

# Call the function with your parameters
ss3x.align(
  #source_path = "oss://skyseq-product/C1579367615888093184/BAZH-20240206-L-01-2024-02-191849/",
  destination_path = "/realspace/project/proj_SS3X_ABio_ZHM_2023_01/Proj_Abio_human_SPM/SH-test/",
  fastq_path = '/realspace/project/proj_SS3X_ABio_ZHM_2023_01/Proj_Abio_human_SPM/SH-test/fastq/',
  #include_pattern = "*SH-test*",
  species = 'human',
  data_type = 'scRNA',
  download = F,
  align = TRUE
)

###数据读取及去重
rt=read.table("/realspace/pe/liyong/ZHM/HLX/results/report//count_matrix(0).csv",sep=",",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

###创建seurat对象
mca <- CreateSeuratObject(counts = data,project = "test", 
                           min.cells = 3, min.features = 50, names.delim = "-",)

###数据过滤和质控
mca[["percent.mt"]] <- PercentageFeatureSet(object = mca, pattern = "^MT-")

pdf(file="./plot/01.featureViolin.pdf",width=10,height=6)
VlnPlot(object = mca, features = c("nFeature_RNA", "nCount_RNA", "percent_mt"), 
        ncol = 3)
dev.off()

mca <- subset(mca, nFeature_RNA > 1000 & percent.mt < 20) ###可根据图01进行调整
# srt <- subset(srt,
#               nFeature_RNA > 1000 &
#                 nFeature_RNA < 10000 &
#                 nCount_RNA > 1000 & 
#                 percent.mt < 20);dim(srt)

pdf(file="./plot/02.featureCor.pdf",width=10,height=6)              
plot1 <- FeatureScatter(object = mca, feature1 = "nCount_RNA", 
                        feature2 = "percent.mt",pt.size=1.5)
plot2 <- FeatureScatter(object = mca, feature1 = "nCount_RNA", 
                        feature2 = "nFeature_RNA",pt.size=1.5)
CombinePlots(plots = list(plot1, plot2))
dev.off()

###数据标准化及PCA降维
mca <- NormalizeData(object = mca)
#mca[['batch']]<-CreateAssayObject(data=mca[['RNA']]@data)
#DefaultAssay(mca)<-"batch"
mca <- FindVariableFeatures(object = mca, nfeatures = 1500) 
###nfeatures可调整，默认值为2000

top10 <- head(x = VariableFeatures(object = mca), 10)

pdf(file="./plot/03.featureVar.pdf",width=10,height=6)              
plot1 <- VariableFeaturePlot(object = mca)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0,
                     ynudge = 0)
CombinePlots(plots = list(plot1, plot2))
dev.off()

mca = CellCycleScoring(mca,s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)

VariableFeatures(mca) <- grep('^RP[SL]',VariableFeatures(mca),
                              value = T,invert = T)
mca <- ScaleData(mca,vars.to.regress = c('nCount_RNA','percent.mt','S.Score','G2M.Score'))
##mca <- ScaleData(mca, features = VariableFeatures(mca), 
##                 vars.to.regress=c('nCount_RNA', 'percent.mt'))

mca <- RunPCA(mca, features = VariableFeatures(mca), npcs=15)
mca <- RunICA(mca, features = VariableFeatures(mca))

pdf(file="./plot/04-1.pcaGene.pdf",width=10,height=8)
VizDimLoadings(object = mca, dims = 1:4, reduction = "pca",nfeatures = 20)
dev.off()

pdf(file="./plot/04-2.icaGene.pdf",width=10,height=8)
VizDimLoadings(object = mca, dims = 1:4, reduction = "ica",nfeatures = 20)
dev.off()

pdf(file="./plot/05-1.PCA.pdf",width=6.5,height=6)
DimPlot(object = mca, reduction = "pca")
dev.off()

pdf(file="./plot/05-2.ICA.pdf",width=6.5,height=6)
DimPlot(object = mca, reduction = "ica")
dev.off()

pdf(file="./plot/06-1.pcaHeatmap.pdf",width=10,height=8)
DimHeatmap(object = mca, dims = 1:4, cells = 500, balanced = TRUE,
           nfeatures = 30,ncol=2)
dev.off()

pdf(file="./plot/06-2.icaHeatmap.pdf",width=10,height=8)
DimHeatmap(object = mca, dims = 1:4, cells = 500, balanced = TRUE,
           nfeatures = 30, ncol=2, reduction = "ica")
dev.off()

pdf(file="./plot/07.elbowpot.pdf",width=10,height=8)
ElbowPlot(mca, ndims=15)
dev.off()

###t-SNE及umap可视化
pcSelect=15 ##可调整15,20,30均可
mca <- FindNeighbors(object = mca, dims = 1:pcSelect)

r = seq(0.1,to=2,by=0.1)
mca <- FindClusters(object = mca, resolution = r)
##resolution根据聚类结果调整
mca <- RunTSNE(object = mca, dims = 1:pcSelect)##若出现Error in .check_tsne_params(nrow(X), 
                                               ##dims = dims, perplexity = perplexity,  : 
                                               ##perplexity is too large for the number of samples, 
                                               ##将perplexity调成5~50之间数值
mca <- RunUMAP(object = mca, dims = 1:pcSelect, umap.method='umap-learn')

logFCfilter=0.5
adjPvalFilter=0.05

markers <- FindAllMarkers(object = mca, only.pos = TRUE,logfc.threshold = logFCfilter)
sig.markers=markers[(abs(as.numeric(as.vector(markers$avg_log2FC)))>logFCfilter & 
                       as.numeric(as.vector(markers$p_val_adj))<adjPvalFilter),]
################split marker save in excel#########################
readRDS("seurat.object.rds") %>%
  FindAllMarkers(only.pos = T) %>%
  split(.$cluster) %>%
  openxlsx::write.xlsx("deg.xlsx")
##################################################################
write.table(sig.markers,file="./table/01.markers.xls",sep="\t",row.names=F,quote=F)

top10 <- sig.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

saveRDS(mca, './mca.rds')

###umap
pdf(file="./plot/08.UMAP.pdf",width=6.5,height=6)
DimPlot(mca, pt.size=1, reduction="umap", label=T, label.size=4) + NoLegend()
dev.off()

###tsne
pdf(file="./plot/09.TSNE.pdf",width=6.5,height=6)
DimPlot(mca, pt.size=1, reduction="tsne", label=T, label.size=4) + NoLegend()
dev.off()

###marker热图
pdf(file="./plot/10.markerHeatmap.pdf",width=12,height=9)
DoHeatmap(object = mca, features = top10$gene) + NoLegend()
dev.off()

###marker散点图
features = c("CD3E", "CD3D", "CD8A", "CD4", "CX3CR1", "GZMH",
             "IL7R", "FCGR3A", "FGFBP2") ##可调整替换为其他marker
pdf(file="./plot/11.markerScatter.pdf",width=10,height=6)
FeaturePlot(object = mca, features = features, cols = c("gray", "purple"))
dev.off()

###GO分析
gene.go <- bitr(rownames(sig.markers), fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = org.Hs.eg.db)
gene.go=gene.go[is.na(gene.go[,"ENTREZID"])==F,]
go.gene = gene.go$ENTREZID

kk <- enrichGO(gene = go.gene,
               OrgDb = org.Hs.eg.db, 
               pvalueCutoff =0.05, 
               qvalueCutoff = 0.05,
               ont="all",
               readable =T)
write.table(kk,file="./table/02.GO.txt",sep="\t",quote=F,row.names = F)

go.data = kk@result
go.data$GO.ID = row.names(go.data)
rownames(go.data) = 1:nrow(go.data)
GO_term_order=factor(as.integer(rownames(go.data)),labels=go.data$Description)
COLS <- c("#66C3A5", "#8DA1CB", "#FD8D62")

###GO统计图
pdf(file="./plot/11.GO_count_bar.pdf",width = 30,height = 15)
ggplot(data=go.data, aes(x=GO_term_order,y=Count, fill=ONTOLOGY)) +
  geom_bar(stat="identity", width=0.8)  + 
  scale_fill_manual(values = COLS) + theme_bw() +
  xlab("GO term") + ylab("Num of Genes") + 
  labs(title = "The Most Enriched GO Terms")+ 
  theme(axis.text.x=element_text(face = "bold",
                                 angle = 80,vjust = 1, hjust = 1 ))
dev.off()

###GO柱形图
pdf(file="./plot/12.GO_barplot.pdf",width = 12,height = 10)
barplot(kk, drop = TRUE, showCategory =10,split="ONTOLOGY", label_format = 100)
dev.off()

###GO气泡图
pdf(file="./plot/13.GO_bubble.pdf",width = 10,height = 8)
dotplot(kk,showCategory =10,split="ONTOLOGY", label_format = 100)
dev.off()

###GO DAG图
kk.bp <- enrichGO(gene = go.gene,
               OrgDb = org.Hs.eg.db, 
               pvalueCutoff =0.05, 
               qvalueCutoff = 0.05,
               ont="BP",
               readable =T)
pdf(file="./plot/14-1.GO_DAG_bp.pdf",width = 8,height = 6)
plotGOgraph(kk.bp)
dev.off()

kk.mf <- enrichGO(gene = go.gene,
                  OrgDb = org.Hs.eg.db, 
                  pvalueCutoff =0.05, 
                  qvalueCutoff = 0.05,
                  ont="MF",
                  readable =T)
pdf(file="./plot/14-2.GO_DAG_mf.pdf",width = 8,height = 6)
plotGOgraph(kk.mf)
dev.off()

kk.cc <- enrichGO(gene = go.gene,
                  OrgDb = org.Hs.eg.db, 
                  pvalueCutoff =0.05, 
                  qvalueCutoff = 0.05,
                  ont="CC",
                  readable =T)
pdf(file="./plot/14-3.GO_DAG_cc.pdf",width = 8,height = 6)
plotGOgraph(kk.cc)
dev.off()

###singleR自动注释
load("/realspace/pe/liyong/XXH_clinical_sequence/SS3X-test/HumanPrimaryCellAtlas_hpca.se_human.RData")
load("/realspace/pe/liyong/XXH_clinical_sequence/SS3X-test/DatabaseImmuneCellExpressionData.RData")

mcaSingleR <- GetAssayData(mca, slot = "data")
mca.hesc <- SingleR(test = mcaSingleR, ref = hpca.se, 
                    labels = hpca.se$label.main) ##大类注释用hpca.se$label.main
                                                 ##小类注释用hpca.se$label.fine
                                                 ##免疫细胞考虑使用dice
mca@meta.data$sR_clusters <- mca.hesc$labels

pdf(file="./plot/15.singleR_anno_umap.pdf",width = 8,height = 6)
DimPlot(mca, group.by = c("seurat_clusters", "sR_clusters"),reduction = "umap")
dev.off()

pdf(file="./plot/16.singleR_anno_tsne.pdf",width = 8,height = 6)
DimPlot(mca, group.by = c("seurat_clusters", "sR_clusters"),reduction = "tsne")
dev.off()

col3 = hcl.colors(12, "Blues 3", rev = TRUE)
per.table = table(mca@meta.data$meta_cluster, mca@meta.data$sR_cluster)
per.table = as.data.frame.array(per.table)
corrplot(per.table, method="circle", is.corr = FALSE, col.lim = c(0,1), tl.col = 'black', col = col3)
###最终结果存为./plot/17.singleR_compare.pdf

###monocle伪时间分析
monocle.matrix=as.matrix(mca@assays$RNA@data)
monocle.matrix=cbind(id=row.names(monocle.matrix),monocle.matrix)
write.table(monocle.matrix,file="./table/03.monocleMatrix.txt",quote=F,sep="\t",row.names=F)
monocle.sample=as.matrix(mca@meta.data)
monocle.sample=cbind(id=row.names(monocle.sample),monocle.sample)
write.table(monocle.sample,file="./table/04.monocleSample.txt",quote=F,sep="\t",row.names=F)
monocle.geneAnn=data.frame(gene_short_name = row.names(monocle.matrix), 
                           row.names = row.names(monocle.matrix))
monocle.geneAnn=cbind(id=row.names(monocle.geneAnn),monocle.geneAnn)
write.table(monocle.geneAnn,file="./table/05.monocleGene.txt",quote=F,sep="\t",row.names=F)
write.table(singler$other,file="./table/06.monocleClusterAnn.txt",quote=F,sep="\t",col.names=F)
write.table(sig.markers,file="./table/07.monocleMarkers.txt",sep="\t",row.names=F,quote=F)

monocle.matrix=read.table("./table/03.monocleMatrix.txt",sep="\t",
                          header=T,row.names=1,check.names=F)
monocle.sample=read.table("./table/04.monocleSample.txt",sep="\t",
                          header=T,row.names=1,check.names=F)
monocle.geneAnn=read.table("./table/05.monocleGene.txt",sep="\t",
                           header=T,row.names=1,check.names=F)
marker=read.table("./table/07.monocleMarkers.txt",sep="\t",
                  header=T,check.names=F)

data.m <- as(as.matrix(monocle.matrix), 'sparseMatrix')
pd<-new("AnnotatedDataFrame", data = monocle.sample)
fd<-new("AnnotatedDataFrame", data = monocle.geneAnn)
cds <- newCellDataSet(data.m, phenoData = pd, featureData = fd)

names(pData(cds))[names(pData(cds))=="seurat_clusters"]="Cluster"
pData(cds)[,"Cluster"]=paste0("cluster",pData(cds)[,"Cluster"])

clusterRt=read.table("./table/06.monocleClusterAnn.txt",
                     header=F,sep="\t",check.names=F)
clusterAnn=as.character(clusterRt[,2])
names(clusterAnn)=paste0("cluster",clusterRt[,1])
pData(cds)$cell_type2 <- plyr::revalue(as.character(pData(cds)$Cluster),clusterAnn)

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

cds <- setOrderingFilter(cds, marker$gene)
plot_ordering_genes(cds)
cds <- reduceDimension(cds, max_components = 2,reduction_method = 'DDRTree')
cds <- orderCells(cds)
pdf(file="./plot/18.cluster.trajectory.pdf",width=6.5,height=6)
plot_cell_trajectory(cds, color_by = "Cluster") ##颜色分类可换
dev.off()

###cellChat
cellchat = createCellChat(object = mca, group.by = 'sR_clusters')
cellchat = setIdent(cellchat, ident.use = 'sR_clusters')
groupSize = as.numeric(table(cellchat@idents))

CellChatDB = CellChatDB.human
colnames(CellChatDB$interaction)
showDatabaseCategory(CellChatDB)
unique(CellChatDB$interaction$annotation)
CellChatDB.use = subsetDB(CellChatDB, search = 'Secreted Signaling')
cellchat@DB = CellChatDB.use
cellchat = subsetData(cellchat)
future::plan("multisession", workers = 4) ##workers数量可根据细胞数量上调
cellchat = identifyOverExpressedGenes(cellchat)
cellchat = identifyOverExpressedInteractions(cellchat)
cellchat = projectData(cellchat, PPI.human)

cellchat = computeCommunProb(cellchat, raw.use = FALSE, population.size = TRUE)
cellchat = filterCommunication(cellchat, min.cells = 10)
df.net = subsetCommunication(cellchat)
write.csv(df.net, file = './table/08.cellchat_net.csv')

cellchat = computeCommunProbPathway(cellchat)
df.netp = subsetCommunication(cellchat, slot.name = 'netP')
write.csv(df.netp, file = './table/09.cellchat_net_pathway.csv')

cellchat = aggregateNet(cellchat)
par(mfrow = c(1,2), xpd=TRUE)

pdf(file="./plot/19.cellchat_number_of_interactions.pdf",width = 8,height = 6)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
                 weight.scale = T, label.edge = F, 
                 title.name = 'Number of interactions')
dev.off()

pdf(file="./plot/20.cellchat_weight_of_interactions.pdf",width = 8,height = 6)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, 
                 weight.scale = T, label.edge = F, 
                 title.name = 'Interactions weights')
dev.off()

mat = cellchat@net$count
par(mfrow = c(1,2), xpd=TRUE)
for(i in 1:nrow(mat)) {
  mat2 = matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] = mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, 
                   arrow.width = 0.2, arrow.size = 0.1, 
                   edge.weight.max = max(mat), 
                   title.name = rownames(mat)[i])
}
##根据数量调整子图，需预览调整，最终保存为./plot/21.cellchat_number_of_cell.pdf

par("mar")
par(mar=c(1,1,1,1))

mat = cellchat@net$weight
par(mfrow = c(1,2), xpd=TRUE)###根据mat数值调整
for(i in 1:nrow(mat)) {
  mat2 = matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] = mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, arrow.width = 0.2,
                   arrow.size = 0.1, edge.weight.max = max(mat), 
                   title.name = rownames(mat)[i])
}
##根据数量调整子图，需预览调整，最终保存为./plot/22.cellchat_weight_of_cell.pdf

cellchat@netP$pathways
pathways.show = c('IL16')
levels(cellchat@idents)
vertex.receiver = c(1,2)
netVisual_aggregate(cellchat, signaling = pathways.show, vertex.receiver = vertex.receiver)
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = 'Reds')
netVisual_bubble(cellchat, sources.use = c(1,2), targets.use = c(1,2), 
                 remove.isolate = FALSE)
###更多画图请参照教程，最终保存为./plot/23-x.cellchat_pathway.pdf，x为对应通路序号

###TCR
immdata = repLoad('/realspace/pe/liyong/XXH_clinical_sequence/P135/meta/test.txt')

###基本数据统计和可视化
repExplore(immdata$data, "lens") %>% vis() 
repExplore(immdata$data, .method = "volume") %>% vis() 
repExplore(immdata$data, .method = "count") %>% vis()

###克隆
repClonality(immdata$data, "homeo") %>% vis() 
repClonality(immdata$data, .method = "clonal.prop") %>% vis()
repClonality(immdata$data, .method = "rare") %>% vis()
repClonality(immdata$data, .method = "top") %>% vis()

###其他
geneUsage(immdata$data[[1]]) %>% vis()
repOverlap(immdata$data, .method = 'overlap') %>% vis()###探索和比较T细胞/B细胞的序列，多样本

repDiversity(immdata$data) %>% vis 
###method('chao1', 'hill', 'div', 'gini.simp','inv.simp', 'gini', 'raref', 'd50', 'dxx')
###根据需求选择
###循环计算
#method = c('chao1', 'hill', 'div', 'gini.simp','inv.simp', 'gini', 'raref', 'd50', 'dxx')
#diversity = list()
#for(i in 1:length(method)){
#  value = repDiversity(immdata$data, .method = method[i])
#  diversity[[method[i]]] = value
#}


###TCR联合转录组
imm_meta = read.table('/realspace/pe/liyong/XXH_clinical_sequence/pipline/immun/TRUST_report.tsv', header = T)
immun = separate(imm_meta, col = cid, into = c('cell.id', 'assemble'), sep = '_a')
immun$chain = substring(immun$V, 1, 3)

immun_A = subset(immun, chain == 'TRA')
immun_B = subset(immun, chain == 'TRB')
immun_D = subset(immun, chain == 'TRD')
immun_G = subset(immun, chain == 'TRG')
immun_A_D = rbind(immun_A, immun_D)
immun_B_G = rbind(immun_B, immun_G)
row.names(immun_A_D) = immun_A_D$cell.id
row.names(immun_B_G) = immun_B_G$cell.id
immun.t = merge(immun_A_D, immun_B_G, all=TRUE, by='row.names', 
                suffixes = c("_A","_B"))

clones = c()
for(i in 1:nrow(immun.t)){
  value = min(immun.t$count.x[i], immun.t$count.y[i])
  clones = append(clones,value)
}
immun.t$Clones = clones
immun.t$cell = 'T'
immun.t = unite(immun.t, col = 'VDJtype', c('cell',"chain.A", 'chain.B'), sep = '_')
imm = unite(imm, col = 'CDR3nt', c('CDR3nt.x',"CDR3nt.y"), sep = ';')
imm = unite(imm, col = 'CDR3aa', c('CDR3aa.x',"CDR3aa.y"), sep = ';')
imm = unite(imm, col = 'V', c('V.x',"V.y"), sep = ';')
imm = unite(imm, col = 'D', c('D.x',"D.y"), sep = ';')
imm = unite(imm, col = 'J', c('J.x',"J.y"), sep = ';')
imm = unite(imm, col = 'C', c('C.x',"C.y"), sep = ';')
imm = unite(imm, col = 'assemble', c('assemble.x',"assemble.y"), sep = ';')
imm = unite(imm, col = 'chain', c('chain.x',"chain.y"), sep = ';')

imm.f = imm[,c('Row.names', 'Clones', 'CDR3nt', 'CDR3aa', 'V', 'D', 'J', 'C', 'assemble',
               'chain')]

clono_seu <- AddMetaData(object=mca, metadata=immun_t)

cloneTypes = c(None = 0, Single = 1, Small = 5, Medium = 20, 
               Large = 100, Hyperexpanded = 500)
clono_seu@meta.data$cloneTypes[clono_seu@meta.data$Clones == 0 | 
                                 clono_seu@meta.data$Clones == 'NA'] <- "None"
clono_seu@meta.data$cloneTypes[0 < clono_seu@meta.data$Clones & 
                                 clono_seu@meta.data$Clones < 1 |
                                 clono_seu@meta.data$Clone == 1] <- "Single"
clono_seu@meta.data$cloneTypes[1 < clono_seu@meta.data$Clones & 
                                 clono_seu@meta.data$Clones < 5 | 
                                 clono_seu@meta.data$Clone == 5] <- "Small"
clono_seu@meta.data$cloneTypes[5 < clono_seu@meta.data$Clones & 
                                 clono_seu@meta.data$Clones < 20 | 
                                 clono_seu@meta.data$Clone == 20] <- "Mediem"
clono_seu@meta.data$cloneTypes[20 < clono_seu@meta.data$Clones & 
                                 clono_seu@meta.data$Clones < 100 | 
                                 clono_seu@meta.data$Clone == 100] <- "Large"
clono_seu@meta.data$cloneTypes[100 < clono_seu@meta.data$Clones] <- "Hyperexpanded"

###STARTRAC示例
vdj <- readRDS('/workspace/zhengliangtao/work/proj_LYM_NJ/zenodo/data/tcr/cHL.tcr.ext.flt.tb.rds')
cell.data<-data.frame(Cell_Name=vdj$cellID,
                      clone.id=vdj$cloneID,
                      patient=vdj$patient,
                      majorCluster=vdj$majorCluster,
                      loc="L")
dim(cell.data)
st<-Startrac.run(cell.data, proj = "HL")
st@cluster.sig.data


###SNP
##first step
#cd SNP_annot_path
#bash merge.snp.sh

##second step
var.file = '/realspace/pe/liyong/ZHM/SNP/result/snp/vcf_anno/sample_merge.txt'
var.annovar.maf = annovarToMaf(annovar = var.file, 
                               Center = 'NA', 
                               refBuild = 'hg38', 
                               tsbCol = 'Sample_Barcode', 
                               table = 'refGene',
                               sep = "\t")

var_maf = read.maf(maf = var.annovar.maf)

plotmafSummary(maf = var_maf, rmOutlier = TRUE, addStat = 'median')
oncoplot(maf = var_maf, top = 30, fontSize = 12 ,showTumorSampleBarcodes = F )
laml.titv = titv(maf = var_maf, plot = FALSE, useSyn = TRUE)
plotTiTv(res = laml.titv)
somaticInteractions(maf = var_maf, top = 25, pvalue = c(0.05, 0.1))
##更多图片需参考maftools画图教程