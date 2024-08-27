
##########################################################################################
#### The entire process of analyzing single-cell RNAseq data of TNBC in GSE161529
##########################################################################################


####merge data####
library(Matrix)
library(Seurat)
library(dplyr)
library(patchwork)
library(foreach)
GSE161529_TNBC_exprMat           <- foreach(i = 81:88,.combine = cbind)%dopar%{
  dir=paste0("Pro_TNBC/data/scRNASeq/52/GSE161529_RAW/GSM49092",i)
  GSM49092i_counts               <- Read10X(
    data.dir=dir,
    gene.column = 1,
    cell.column = 2,
    unique.features = TRUE,
    strip.suffix = FALSE)
  GSM49092i_exprMat              <- as.matrix(GSM49092i_counts)
}
Sample=c(rep("TNBC_0126",3666),rep("TNBC_0135",15870),rep("TNBC_0106",1065),rep("TNBC_0114",2015),rep("TNBC_4031",5581),rep("TNBC_0131",6456),rep("TNBC_0554",9593),rep("TNBC_0177",21130))
GSE161529_TNBC                           <- data.frame(cell.id=colnames(GSE161529_TNBC_exprMat),Sample=Sample)
GSE161529_TNBC$cell.name                 <- apply(GSE161529_TNBC,1,function(x){paste(x[2],x[1],sep = "_")})
colnames(GSE161529_TNBC_exprMat)         <- GSE161529_TNBC$cell.name
GSE161529_TNBC_scRNA                     <- CreateSeuratObject(counts = GSE161529_TNBC_exprMat)
GSE161529_TNBC_scRNA$orig.ident          <- GSE161529_TNBC$Sample
#calculate the proportion of ribosomal genes in cells
GSE161529_TNBC_scRNA[["percent.mt"]]     <- PercentageFeatureSet(GSE161529_TNBC_scRNA,pattern = "^MT-")
head(GSE161529_TNBC_scRNA@meta.data,5)
VlnPlot(GSE161529_TNBC_scRNA,features = c("nCount_RNA","nFeature_RNA","percent.mt"),ncol=3)
GSE161529_TNBC_scRNA                     <- subset(GSE161529_TNBC_scRNA,nFeature_RNA >0 & nFeature_RNA<8000)
#normalized
GSE161529_TNBC_scRNA                     <- NormalizeData(GSE161529_TNBC_scRNA)
#features choosing
GSE161529_TNBC_scRNA                     <- FindVariableFeatures(GSE161529_TNBC_scRNA,selection.method = "vst",nfeatures = 2000)#choose high variability between cells.
#identify the 10 most highly variable genes
top10                                    <- head(VariableFeatures(GSE161529_TNBC_scRNA),10)
#genes related to cell cycle
cellcycle_gene                           <- c(cc.genes$s.genes,cc.genes$g2m.genes)
head(cellcycle_gene)
CaseMatch(cellcycle_gene,VariableFeatures(GSE161529_TNBC_scRNA))#check which cell cycle-related genes are among the hypervariable genes we selected.
# results:character(0)
#data scaling 
all_gene                                   <- rownames(GSE161529_TNBC_scRNA)
GSE161529_TNBC_scRNA                       <- ScaleData(GSE161529_TNBC_scRNA,features = all_gene)#scaled data存放在 pbmc[["RNA"]]@scale.data
#linear dimensionality reduction
GSE161529_TNBC_scRNA                       <- RunPCA(GSE161529_TNBC_scRNA,features = VariableFeatures(object=GSE161529_TNBC_scRNA))
#dimension selection
GSE161529_TNBC_scRNA                       <- JackStraw(GSE161529_TNBC_scRNA,num.replicate = 100)
GSE161529_TNBC_scRNA                       <- ScoreJackStraw(GSE161529_TNBC_scRNA,dims = 1:20)
JackStrawPlot(GSE161529_TNBC_scRNA,dims=1:20)
ElbowPlot(GSE161529_TNBC_scRNA,ndims = 50)#choose dims of 30
#cell clustering
GSE161529_TNBC_scRNA                       <- FindNeighbors(GSE161529_TNBC_scRNA,dims = 1:40)
GSE161529_TNBC_scRNA                       <- FindClusters(GSE161529_TNBC_scRNA,resolution = 0.1)
table(Idents(GSE161529_TNBC_scRNA))
#unlinear dimension reduction(UMAP or tSNE)
GSE161529_TNBC_scRNA                        <- RunUMAP(GSE161529_TNBC_scRNA,dims = 1:40)
GSE161529_TNBC_scRNA                        <- RunTSNE(GSE161529_TNBC_scRNA,dims = 1:40,check_duplicates = FALSE)
DimPlot(GSE161529_TNBC_scRNA,reduction = "umap",label = T)
DimPlot(GSE161529_TNBC_scRNA,reduction = "tsne",label = T)
GSE161529_TNBC_uamp                         <- as.data.frame(GSE161529_TNBC_scRNA@reductions$umap@cell.embeddings)
GSE161529_TNBC_uamp$sample.id               <- rownames(GSE161529_TNBC_uamp)
GSE161529_TNBC_uamp                         <- merge(GSE161529_TNBC_uamp,GSE161529_TNBC_celltype,by="sample.id")
ggplot(data = GSE161529_TNBC_uamp,aes(x=UMAP_1,y=UMAP_2,col=celltype))+
  geom_point()+ggtitle("UAMP plot of GSE161529_TNBC by using all genes")
GSE161529_TNBC_uamp                         <- merge(GSE161529_TNBC_uamp,GSE161529_TNBC_metadata,by="sample.id")



####Cell annotation####
p12                                  <- DimPlot(GSE161529_TNBC_scRNA,reduction = "umap",label = T)
p14                                  <- FeaturePlot(GSE161529_TNBC_scRNA,features = "ENSG00000119888",reduction = "umap",pt.size = 1)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())+ggtitle("EPCAM:epithelial")
p15                                  <- FeaturePlot(GSE161529_TNBC_scRNA,features = "ENSG00000164692",reduction = "umap",pt.size = 1)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())+ggtitle("COL1A2:fibroblasts")
p16                                  <- FeaturePlot(GSE161529_TNBC_scRNA,features = "ENSG00000134853",reduction = "umap",pt.size = 1)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())+ggtitle("PDGFRA:fibroblasts")
p17                                  <- FeaturePlot(GSE161529_TNBC_scRNA,features = "ENSG00000177575",reduction = "umap",pt.size = 1)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())+ggtitle("CD163:macrophage")#macrophage marker
p18                                  <- FeaturePlot(GSE161529_TNBC_scRNA,features = "ENSG00000010610",reduction = "umap",pt.size = 1)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())+ggtitle("CD4:T-cell")#T-cell
p22                                  <-  FeaturePlot(GSE161529_TNBC_scRNA,features = "ENSG00000153563",reduction = "umap",pt.size = 1)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())+ggtitle("CD8A:T-cell")#T-cell
p23                                  <- FeaturePlot(GSE161529_TNBC_scRNA,features = "ENSG00000116824",reduction = "umap",pt.size = 1)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())+ggtitle("CD2:T-cell")#T-cell
p19                                  <- FeaturePlot(GSE161529_TNBC_scRNA,features = "ENSG00000261371",reduction = "umap",pt.size = 1)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())+ggtitle("PECAM1:endothlial")#endothlial marker
p20                                  <- FeaturePlot(GSE161529_TNBC_scRNA,features = "ENSG00000177455",reduction = "umap",pt.size = 1)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())+ggtitle("CD19:B-cell")#B-cell marker
p21                                  <- FeaturePlot(GSE161529_TNBC_scRNA,features = "ENSG00000139193",reduction = "umap",pt.size = 1)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())+ggtitle("CD27:B-cell")#B-cell marker
p24                                  <- FeaturePlot(GSE161529_TNBC_scRNA,features = "ENSG00000211685",reduction = "umap",pt.size = 1)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())+ggtitle("IGLC7:plasma")#plasma cells marker
p25                                  <- FeaturePlot(GSE161529_TNBC_scRNA,features = "ENSG00000211895",reduction = "umap",pt.size = 1)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())+ggtitle("IGHA1:plasma")#plasma cells marker
p26                                  <- FeaturePlot(GSE161529_TNBC_scRNA,features = "ENSG00000110799",reduction = "umap",pt.size = 1)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())+ggtitle("VWF:endothelial")#plasma cells marker
p27                                  <- FeaturePlot(GSE161529_TNBC_scRNA,features = "ENSG00000174059",reduction = "umap",pt.size = 1)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())+ggtitle("CD34:TAMs")#TAMs  marker
p28                                  <- FeaturePlot(GSE161529_TNBC_scRNA,features = "ENSG00000129226",reduction = "umap",pt.size = 1)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())+ggtitle("CD68:TAMs")#TAMs  marker
p29                                  <- FeaturePlot(GSE161529_TNBC_scRNA,features = "ENSG00000170458",reduction = "umap",pt.size = 1)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())+ggtitle("CD14:TAMs")#TAMs  marker
p30                                  <- FeaturePlot(GSE161529_TNBC_scRNA,features = "ENSG00000076706",reduction = "umap",pt.size = 1)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())+ggtitle("MCAM:pericytes")#
p31                                  <- FeaturePlot(GSE161529_TNBC_scRNA,features = "ENSG00000148773",reduction = "umap",pt.size = 1)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())+ggtitle("MKI67:cycling_cell")#
p32                                  <- FeaturePlot(GSE161529_TNBC_scRNA,features = "ENSG00000175084",reduction = "umap",pt.size = 1)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())+ggtitle("DES:pericytes")#
p33                                  <- FeaturePlot(GSE161529_TNBC_scRNA,features = "ENSG00000149591",reduction = "umap",pt.size = 1)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())+ggtitle("TAGLN:pericytes")#
p34                                  <- FeaturePlot(GSE161529_TNBC_scRNA,features = "ENSG00000074181",reduction = "umap",pt.size = 1)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())+ggtitle("Notch3:pericytes")#
figure_1                             <- p12+p14+p15+p16+p17+p18+p22+p23+p19+p26+p20+p21+p24+p25+p27+p28+p29+p31+p30+p32+p33+p34 
ggsave(figure_1,filename = "Pro_TNBC/output/data/scRNASeq/38sample/GSE161529_TNBC/celltype_marker_gene.pdf",height = 30,width = 30  )

####cell anatation by using SingleR
library(clusterProfiler)
library(SingleR)
load("~/Pro_TNBC/output/data/scRNASeq/cellclassify.refdata.RData")
testdata                  <- GetAssayData(GSE161529_TNBC_scRNA, slot="data")
testdata                  <- as.data.frame(testdata)
testdata$ENSEMBL          <- rownames(testdata)
gene                      <- bitr(testdata$ENSEMBL,fromType = "ENSEMBL",toType = "SYMBOL",OrgDb = "org.Hs.eg.db")
testdata                  <- merge(gene,testdata,by="ENSEMBL")
rownames(testdata)        <- make.names(testdata[,2],TRUE)#Allow duplicate row names.
library(dplyr)
testdata                  <- testdata[,-c(1,2)] %>%as.matrix()
testdata                  <- CreateSeuratObject(counts = testdata)
testdata                  <- GetAssayData(testdata, slot="data")
clusters                  <- GSE161529_TNBC_scRNA@meta.data$seurat_clusters
cellpred                  <- SingleR(test = testdata, ref = refdata, labels = refdata$label.main, 
                                     method = "cluster", clusters = clusters, 
                                     assay.type.test = "logcounts", assay.type.ref = "logcounts")
celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)
celltype[12,2]            <- "T_cells"
celltype[13,2]            <- "Pericyte"
celltype[3,2]             <- "TAM"
GSE161529_TNBC_scRNA@meta.data$celltype = "NA"
#annotation
for(i in 1:nrow(celltype)){
  GSE161529_TNBC_scRNA@meta.data[which(GSE161529_TNBC_scRNA@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
#visualization
DimPlot(GSE161529_TNBC_scRNA, group.by="celltype", label=T, label.size=5, reduction='umap')
save(GSE161529_TNBC_scRNA,file = "Pro_TNBC/output/data/scRNASeq/38sample/GSE161529_TNBC/GSE161529_TNBC_scRNA.RData")


####infercnv####
library(Seurat)
library(tidyverse)
library(ggplot2)
library(infercnv)
library(ComplexHeatmap)
library(ggpubr)
#inferCNV needs three files:1.count expression matrix,2.group information,3.Gene chromosome information
#make gene chromosome position information and extract expression matrix.
dat             <- GetAssayData(GSE161529_TNBC_scRNA,assay = "RNA",slot = "count") 
dat             <- as.data.frame(dat)
library(AnnoProbe)  
geneInfor=annoGene(rownames(dat),"ENSEMBL",'human')  
colnames(geneInfor)
geneInfor=geneInfor[with(geneInfor, order(chr, start)),c(3,4:6)]      
geneInfor=geneInfor[!duplicated(geneInfor[,1]),]
length(unique(geneInfor[,1]))
head(geneInfor)
dat                 <- dat[match(geneInfor[,1],rownames(dat)),]
geneInfor_rown      <- geneInfor$ENSEMBL
geneInfor           <- geneInfor[,-1]
rownames(geneInfor) <- geneInfor_rown
meta                <- subset(GSE161529_TNBC_scRNA@meta.data,select = c("celltype"))

#inferCNV
#Two-step construction object
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=dat,
                                    annotations_file=meta,
                                    delim="\t",
                                    gene_order_file=geneInfor,
                                    ref_group_names = c("T_cells"))   #Select the basic cell or sample to see the input type of meta, or do not select the algorithm to calculate by yourself according to the average value.
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir="Pro_TNBC/output/data/scRNASeq/38sample/GSE161529_TNBC/infercnv", 
                             cluster_by_groups=TRUE,  # If TRUE is selected, group by sample will be changed to FALSE, and another parameter k will be selected_ obs_ The number of groups given by groups (the default is 1) for grouping.
                             denoise=T,     
                             HMM=F,no_plot = T,num_threads=20)   # #Whether to predict CNV based on HMM? True is a long time
#Finally, many files will be output in out_dir. You can directly use the heat map inside.

####choose tumor cell
library(infercnv)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
run.final.infercnv_obj = readRDS("Pro_TNBC/output/data/scRNASeq/38sample/GSE161529_TNBC/infercnv/run.final.infercnv_obj")
expr          <- run.final.infercnv_obj@expr.data
normal_loc    <- run.final.infercnv_obj@reference_grouped_cell_indices
normal_loc    <- normal_loc$T_cells
test_loc      <- run.final.infercnv_obj@observation_grouped_cell_indices
test_loc      <- c(test_loc$B_cell,test_loc$CMP,test_loc$Pericyte,test_loc$Endothelial_cells,test_loc$Erythroblast,test_loc$Epithelial_cells,test_loc$Fibroblasts,test_loc$TAM)
anno.df=data.frame(
  CB=c(colnames(expr)[normal_loc],colnames(expr)[test_loc]),
  class=c(rep("normal",length(normal_loc)),rep("test",length(test_loc)))
)
head(anno.df)
gn            <- rownames(expr)
geneFile      <- geneInfor
sub_geneFile  <- geneFile[intersect(gn,rownames(geneFile)),]
expr=expr[intersect(gn,rownames(geneFile)),]
head(sub_geneFile,4)
expr[1:4,1:4]
#clustering
set.seed(20210418)
kmeans.result      <- kmeans(t(expr), 7)
kmeans_df          <- data.frame(kmeans_class=kmeans.result$cluster)
kmeans_df$CB=rownames(kmeans_df)
kmeans_df=kmeans_df%>%inner_join(anno.df,by="CB") #combination
kmeans_df_s=arrange(kmeans_df,kmeans_class) #ordering
rownames(kmeans_df_s)=kmeans_df_s$CB
kmeans_df_s$CB=NULL
kmeans_df_s$kmeans_class=as.factor(kmeans_df_s$kmeans_class) #将kmeans_class转换为因子，作为热图的一个注释，最终在热图上就会按照1:7排序
head(kmeans_df_s)

#Define annotations and color schemes for heat maps
top_anno          <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "NA",col="NA"), labels = 1:22,labels_gp = gpar(cex = 1.5)))
color_v=RColorBrewer::brewer.pal(8, "Dark2")[1:7] #类别数
names(color_v)=as.character(1:7)
left_anno         <- rowAnnotation(df = kmeans_df_s,col=list(class=c("test"="red","normal" = "blue"),kmeans_class=color_v))
#draw
pdf("Pro_TNBC/output/data/scRNASeq/38sample/GSE161529_TNBC/infercnv/try1.pdf",width = 15,height = 10)
ht = Heatmap(t(expr)[rownames(kmeans_df_s),],
             col = colorRamp2(c(0.4,1,1.6), c("#377EB8","#F0F0F0","#E41A1C")), 
             cluster_rows = F,cluster_columns = F,show_column_names = F,show_row_names = F,
             column_split = factor(sub_geneFile$chr, paste("chr",1:22,sep = "")), 
             
             heatmap_legend_param = list(title = "Modified expression",direction = "vertical",title_position = "leftcenter-rot",at=c(0.4,1,1.6),legend_height = unit(3, "cm")),
             
             top_annotation = top_anno,left_annotation = left_anno, #添加注释
             row_title = NULL,column_title = NULL)
draw(ht, heatmap_legend_side = "right")
dev.off()
#
write.csv(kmeans_df_s,file="Pro_TNBC/output/data/scRNASeq/38sample/GSE161529_TNBC/infercnv/kmeans_df_s.csv")
GSE161529_TNBC_scRNA@meta.data$kmeans_cluster <- "NA"
for (i in 1:nrow(kmeans_df_s)){
  GSE161529_TNBC_scRNA@meta.data[which(rownames(GSE161529_TNBC_scRNA@meta.data)==rownames(kmeans_df_s)[i]),'kmeans_cluster'] <- kmeans_df_s$kmeans_class[i]                            
}

GSE161529_TNBC_scRNA@meta.data$cell <- "NA"
GSE161529_TNBC_scRNA@meta.data$cell[
                                      GSE161529_TNBC_scRNA@meta.data$kmeans_cluster== "1" |
                                      GSE161529_TNBC_scRNA@meta.data$kmeans_cluster=="6"|
                                      GSE161529_TNBC_scRNA@meta.data$kmeans_cluster=="7"|
                                      GSE161529_TNBC_scRNA@meta.data$kmeans_cluster=="2"|
                                      GSE161529_TNBC_scRNA@meta.data$kmeans_cluster=="4"
] <- "tumorcell"
GSE161529_TNBC_scRNA@meta.data$cell[GSE161529_TNBC_scRNA@meta.data$kmeans_cluster== "5"|
                                      GSE161529_TNBC_scRNA@meta.data$kmeans_cluster=="3"                                   
] <- "normalcell"
View(GSE161529_TNBC_scRNA@meta.data)
tumorcell     <- GSE161529_TNBC_scRNA@meta.data[GSE161529_TNBC_scRNA@meta.data$cell=="tumorcell",]
table(tumorcell$celltype)
save(GSE161529_TNBC_scRNA,file = "Pro_TNBC/output/data/scRNASeq/38sample/GSE161529_TNBC/GSE161529_TNBC_scRNA.RData")
