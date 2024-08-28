####################################################################################
####Identification of Claudin-low cancer cells in Basal-like cell line and patient
####################################################################################


##########################################
####analyzing HDQP1 in GSE173634
##########################################

library(Matrix)
library(Seurat)
library(dplyr)
library(patchwork)
library(GSVA)
load("~/Pro_TNBC/output/data/CCLE/GSE173634_scRNA.RData")
load("Pro_TNBC/paper/data/results/section_1/UBS93.gene.df.RData")
load("Pro_TNBC/paper/data/method/UBS93.data.RData")
HDQP1_scRNA                     <- subset(GSE173634_scRNA,orig.ident=="HDQP1")
#calculate the proportion of ribosomal genes in cells
HDQP1_scRNA[["percent.mt"]]     <- PercentageFeatureSet(HDQP1_scRNA,pattern = "^MT-")
head(HDQP1_scRNA@meta.data,5)
VlnPlot(HDQP1_scRNA,features = c("nCount_RNA","nFeature_RNA","percent.mt"),ncol=3)
HDQP1_scRNA                     <- subset(HDQP1_scRNA,nFeature_RNA >1000 & nFeature_RNA<10000 )
#normalized
HDQP1_scRNA                     <- NormalizeData(HDQP1_scRNA)
HDQP1_scRNA[["RNA"]]@data[c("ENSG00000167286", "ENSG00000100721", "ENSG00000156738"), 1:10]#data after normalizing
#features choosing
HDQP1_scRNA                     <- FindVariableFeatures(HDQP1_scRNA,selection.method = "vst",nfeatures = 2000)#choose high variability between cells.
#identify the 10 most highly variable genes
top10                           <- head(VariableFeatures(HDQP1_scRNA),10)
#genes related to cell cycle
cellcycle_gene                  <- c(cc.genes$s.genes,cc.genes$g2m.genes)
head(cellcycle_gene)
CaseMatch(cellcycle_gene,VariableFeatures(HDQP1_scRNA))#check which cell cycle-related genes are among the hypervariable genes we selected.
# results:character(0)
#data scaling 
all_gene                        <- rownames(HDQP1_scRNA)
HDQP1_scRNA                     <- ScaleData(HDQP1_scRNA,features = all_gene)
#linear dimensionality reduction
HDQP1_scRNA                     <- RunPCA(HDQP1_scRNA,features = VariableFeatures(object=HDQP1_scRNA))
#dimension selection
HDQP1_scRNA                     <- JackStraw(HDQP1_scRNA,num.replicate = 100)
HDQP1_scRNA                     <- ScoreJackStraw(HDQP1_scRNA,dims = 1:20)
JackStrawPlot(HDQP1_scRNA,dims=1:20)
ElbowPlot(HDQP1_scRNA,ndims = 50)#choose dims of 30
#cell clustering
HDQP1_scRNA                      <- FindNeighbors(HDQP1_scRNA,dims = 1:50)
HDQP1_scRNA                      <- FindClusters(HDQP1_scRNA,resolution = 0.8)
table(Idents(HDQP1_scRNA))

####*TSNE####
HDQP1_exp                        <- as.matrix(HDQP1_scRNA@assays$RNA@data)
HDQP1_exp_UBS93                  <- HDQP1_exp[rownames(HDQP1_exp) %in% UBS93.gene.df$ENSEMBL,]
library(Rtsne)
cor_matrix                       <- cor(HDQP1_exp_UBS93,method = "spearman")
dissimilarity_matrix             <- 1 - cor_matrix
normalized_dissimilarity_matrix  <- dissimilarity_matrix / max(dissimilarity_matrix)
HDQP1_scRNA                      <- RunTSNE(
  HDQP1_scRNA,
  reduction = "pca",
  cells = NULL,
  dims = 1:30,
  features = NULL,
  seed.use = 1,
  tsne.method = "Rtsne",
  dim.embed = 2,
  distance.matrix = normalized_dissimilarity_matrix,
  reduction.name = "tsne",
  reduction.key = "tSNE_"
)

####*run UBS93 predicting for HDQP1 data####
HDQP1_exp                           <- as.matrix(HDQP1_scRNA@assays$RNA@data) #normalied data
HDQP1_cell_subtype                  <- breast.cancer.predictor(expr.of.sample = HDQP1_exp,
                                                               expr.of.centroid = UBS93.data$UBS93.centroid,
                                                               marker.gene = UBS93.data$UBS93.gene.df$ENSEMBL,
                                                               HER2.amp.signature.genes = UBS93.data$HER2.amp.signature.genes,
                                                               ssgsea.cutoff = UBS93.data$ssgsea.cutoff)

HDQP1_subtype                       <- HDQP1_cell_subtype$subtype
HDQP1_subtype$sample.id             <- rownames(HDQP1_subtype)
GSE173634_HDQP1_metadata            <- HDQP1_scRNA@meta.data
GSE173634_HDQP1_metadata$sample.id  <- rownames(GSE173634_HDQP1_metadata)
GSE173634_HDQP1_metadata            <- merge(GSE173634_HDQP1_metadata,HDQP1_subtype,by="sample.id")
rown                                <- GSE173634_HDQP1_metadata$sample.id
GSE173634_HDQP1_metadata            <- GSE173634_HDQP1_metadata[,-1]
rownames(GSE173634_HDQP1_metadata)  <- rown
HDQP1_scRNA@meta.data               <- GSE173634_HDQP1_metadata
save(HDQP1_scRNA,file = "Pro_TNBC/paper/data/section_4/GSE173634_HDQP1_scRNA.RData")

####*Figure 4a: the t-SNE plot for HDQP1 data####
df1                              <- as.data.frame(HDQP1_scRNA@reductions$tsne@cell.embeddings)
df1$cell.id                      <- rownames(df1)
subtype                          <- as.data.frame(HDQP1_scRNA@meta.data)
subtype$cell.id                  <- rownames(subtype)
df1                              <- merge(df1,subtype,by="cell.id")
fig_4a                           <- ggplot(data = df1,aes(x=tSNE_1,y=tSNE_2,col=subtype))+
  geom_point(size=5)+
  ggplot.style+theme(plot.title = element_text(hjust=0.5,size = 44))+scale_color_manual(values = c(Basal="Red",Claudin_low="Blue"))
ggsave(fig_4a,filename = "Pro_TNBC/paper/plot/section_4/tsne.plot.of.HDQP1.by.using.all.genes.pdf",width=23,height=15)


####*Figure 4b: a heatmap for the expression values of 5 marker genes in HDQP1 data ####
HDQP1_exp_df                <- as.data.frame(HDQP1_exp)
HDQP1_exp_df$ENSEMBL        <- rownames(HDQP1_exp_df) 
symbol                      <- bitr(HDQP1_exp_df$ENSEMBL,fromType = "ENSEMBL",toType = "SYMBOL",OrgDb="org.Hs.eg.db")
HDQP1_exp_df                <- merge(symbol,HDQP1_exp_df,by="ENSEMBL")
marker.gene                 <- c("CLDN3","CLDN4","CLDN7","VIM","EPCAM")
HDQP1_exp_df                <- HDQP1_exp_df[HDQP1_exp_df$SYMBOL %in% marker.gene,]
rown                        <- HDQP1_exp_df$SYMBOL
HDQP1_exp_df                <- HDQP1_exp_df[,-c(1,2)]
rownames(HDQP1_exp_df)      <- rown
HDQP1_exp_matrix            <- as.matrix(HDQP1_exp_df)
identical(HDQP1_subtype$sample.id,colnames(HDQP1_exp_matrix))
annotation.col              <- HDQP1_subtype[,1,drop=F]
rownames(annotation.col)    <- HDQP1_subtype$sample.id 
annot.colors                <- list(subtype=c("Basal"="Red","Claudin_low"= "Blue"))
library(pheatmap)
fig4b          <- pheatmap(HDQP1_exp_matrix,cluster_rows = F,
                           cluster_cols = T,
                           color = colorRampPalette(c("navy","white","firebrick3"))(100),
                           show_colnames = F,border_color = NA,
                           scale = "row",
                           show_rownames = T,annotation_col = annotation.col,
                           fontsize=30,cellwidth = 8, cellheight = 55,
                           annotation_colors = annot.colors
)
#fontsize_row = 8,fontsize_col=12 ,
ggsave(fig4b,filename = "Pro_TNBC/paper/plot/section_4/heat.map.of.marker.gene.in.HDQP1.GSE173634.pdf",width = 20,height = 15)

####*Figure S10: box plots for the expression values of marker genes in every subtype#### 
HDQP1_exp_df             <- as.matrix(HDQP1_exp_df) %>% t() %>% as.data.frame()
HDQP1_exp_df$sample.id   <- rownames(HDQP1_exp_df)
HDQP1_exp_df             <- merge(HDQP1_exp_df,HDQP1_subtype,by="sample.id")
supfig_11a               <- ggplot(HDQP1_exp_df,aes(x=subtype,y=VIM))+
  geom_boxplot(outlier.shape=NA,lwd=1.2)+
  geom_point(size=4,position="jitter") + ggplot.style + 
  xlab('')+
  theme(panel.grid.major=element_line(colour=NA),
        panel.grid.minor = element_blank())+
  scale_x_discrete("", labels = c("Basal" = paste("Basal",sep = "\n","(n=58)"),"Claudin_low" = paste("Cladin-low",sep = "\n","(n=58)")))
wilcox.test(VIM ~subtype,data = HDQP1_exp_df)#p-value = 1.83e-14
ggsave(supfig_11a,filename = "Pro_TNBC/paper/plot/section_4/HDQP1_GSE173634/boxplot.VIM.pdf",width = 20,height = 15)
supfig_11b               <- ggplot(HDQP1_exp_df,aes(x=subtype,y=EPCAM))+
  geom_boxplot(outlier.shape=NA,lwd=1.2)+
  geom_point(size=4,position="jitter") + ggplot.style + 
  xlab('')+
  theme(panel.grid.major=element_line(colour=NA),
        panel.grid.minor = element_blank())+
  scale_x_discrete("", labels = c("Basal" = paste("Basal",sep = "\n","(n=58)"),"Claudin_low" = paste("Cladin-low",sep = "\n","(n=58)")))
wilcox.test(EPCAM ~subtype,data = HDQP1_exp_df)#p-value = 4.123e-13
ggsave(supfig_11b,filename = "Pro_TNBC/paper/plot/section_4/HDQP1_GSE173634/boxplot.EPCAM.pdf",width = 20,height = 15)
supfig_11c               <- ggplot(HDQP1_exp_df,aes(x=subtype,y=CLDN3))+
  geom_point(size=4,position="jitter") + ggplot.style + 
  xlab('')+
  theme(panel.grid.major=element_line(colour=NA),
        panel.grid.minor = element_blank())+
  scale_x_discrete("", labels = c("Basal" = paste("Basal",sep = "\n","(n=58)"),"Claudin_low" = paste("Cladin-low",sep = "\n","(n=58)")))
wilcox.test(CLDN3 ~subtype,data = HDQP1_exp_df)#p-value = 9.982e-05
ggsave(supfig_11c,filename = "Pro_TNBC/paper/plot/section_4/HDQP1_GSE173634/boxplot.CLDN3.pdf",width = 20,height = 15)

supfig_11d               <- ggplot(HDQP1_exp_df,aes(x=subtype,y=CLDN4))+
  geom_point(size=4,position="jitter") + ggplot.style + 
  xlab('')+
  theme(panel.grid.major=element_line(colour=NA),
        panel.grid.minor = element_blank())+
  scale_x_discrete("", labels = c("Basal" = paste("Basal",sep = "\n","(n=58)"),"Claudin_low" = paste("Cladin-low",sep = "\n","(n=58)")))
wilcox.test(CLDN4 ~subtype,data = HDQP1_exp_df)#p-value = 5.771e-09
ggsave(supfig_11d,filename = "Pro_TNBC/paper/plot/section_4/HDQP1_GSE173634/boxplot.CLDN4.pdf",width = 20,height = 15)

supfig_11e               <- ggplot(HDQP1_exp_df,aes(x=subtype,y=CLDN7))+
  geom_point(size=4,position="jitter") + ggplot.style + 
  xlab('')+
  theme(panel.grid.major=element_line(colour=NA),
        panel.grid.minor = element_blank())+
  scale_x_discrete("", labels = c("Basal" = paste("Basal",sep = "\n","(n=58)"),"Claudin_low" = paste("Cladin-low",sep = "\n","(n=58)")))
wilcox.test(CLDN7 ~subtype,data = HDQP1_exp_df)#p-value = 6.104e-16
ggsave(supfig_11e,filename = "Pro_TNBC/paper/plot/section_4/HDQP1_GSE173634/boxplot.CLDN7.pdf",width = 20,height = 15)


###########################################
####analyzing HDQP1 in GSE173634
###########################################

load("Pro_TNBC/output/data/scRNASeq/GSE202771/GSE202771_scRNA.RData")
HDQP1_scRNA                     <- subset(GSE202771_scRNA,orig.ident=="HDQP1_B2_scRNA")
#calculate the proportion of ribosomal genes in cells
HDQP1_scRNA[["percent.mt"]]     <- PercentageFeatureSet(HDQP1_scRNA,pattern = "^MT-")
head(HDQP1_scRNA@meta.data,5)
VlnPlot(HDQP1_scRNA,features = c("nCount_RNA","nFeature_RNA","percent.mt"),ncol=3)
HDQP1_scRNA                     <- subset(HDQP1_scRNA,nFeature_RNA >1000 & nFeature_RNA<5000&percent.mt < 20 )
#normalized
HDQP1_scRNA                     <- NormalizeData(HDQP1_scRNA)
HDQP1_scRNA[["RNA"]]@data[c("ENSG00000167286", "ENSG00000100721", "ENSG00000156738"), 1:10]#data after normalizing
#features choosing
HDQP1_scRNA                     <- FindVariableFeatures(HDQP1_scRNA,selection.method = "vst",nfeatures = 2000)#choose high variability between cells.
#identify the 10 most highly variable genes
top10                           <- head(VariableFeatures(HDQP1_scRNA),10)
#genes related to cell cycle
cellcycle_gene                  <- c(cc.genes$s.genes,cc.genes$g2m.genes)
head(cellcycle_gene)
CaseMatch(cellcycle_gene,VariableFeatures(HDQP1_scRNA))#check which cell cycle-related genes are among the hypervariable genes we selected.
# results:character(0)
#data scaling 
all_gene                        <- rownames(HDQP1_scRNA)
HDQP1_scRNA                     <- ScaleData(HDQP1_scRNA,features = all_gene)#scaled data存放在 pbmc[["RNA"]]@scale.data
#linear dimensionality reduction
HDQP1_scRNA                     <- RunPCA(HDQP1_scRNA,features = VariableFeatures(object=HDQP1_scRNA))
#dimension selection
HDQP1_scRNA                     <- JackStraw(HDQP1_scRNA,num.replicate = 100)
HDQP1_scRNA                     <- ScoreJackStraw(HDQP1_scRNA,dims = 1:20)
JackStrawPlot(HDQP1_scRNA,dims=1:20)
ElbowPlot(HDQP1_scRNA,ndims = 50)#choose dims of 30
#cell clustering
HDQP1_scRNA                     <- FindNeighbors(HDQP1_scRNA,dims = 1:40)
HDQP1_scRNA                     <- FindClusters(HDQP1_scRNA,resolution = 0.1)
table(Idents(HDQP1_scRNA))


####*TSNE####
HDQP1_exp                        <- as.matrix(HDQP1_scRNA@assays$RNA@data)
HDQP1_exp_UBS93                  <- HDQP1_exp[rownames(HDQP1_exp) %in% UBS93.gene.df$SYMBOL,]
library(Rtsne)
cor_matrix                       <- cor(HDQP1_exp_UBS93,method = "spearman")
dissimilarity_matrix             <- 1 - cor_matrix
normalized_dissimilarity_matrix  <- dissimilarity_matrix / max(dissimilarity_matrix)
HDQP1_scRNA      <- RunTSNE(
  HDQP1_scRNA,
  reduction = "pca",
  cells = NULL,
  dims = 1:30,
  features = NULL,
  seed.use = 1,
  tsne.method = "Rtsne",
  dim.embed = 2,
  distance.matrix = normalized_dissimilarity_matrix,
  reduction.name = "tsne",
  reduction.key = "tSNE_"
)

####*run UBS93 predicting for HDQP1 data####
HDQP1_exp                         <- as.matrix(HDQP1_scRNA@assays$RNA@data) #normalied data
HDQP1_exp_UBS93                   <- as.data.frame(HDQP1_exp)
HDQP1_exp_UBS93$SYMBOL            <- rownames(HDQP1_exp_UBS93)
HDQP1_exp_UBS93                   <- merge(UBS93.data$UBS93.gene.df,HDQP1_exp_UBS93,by="SYMBOL")
rown                              <- HDQP1_exp_UBS93$ENSEMBL
HDQP1_exp_UBS93                   <- HDQP1_exp_UBS93[,-c(1:3)]
rownames(HDQP1_exp_UBS93)         <- rown
HDQP1_exp_UBS93                   <- as.matrix(HDQP1_exp_UBS93)
HDQP1_cell_subtype                <- breast.cancer.predictor(expr.of.sample = HDQP1_exp_UBS93,
                                                             expr.of.centroid = UBS93.data$UBS93.centroid,
                                                             marker.gene = UBS93.data$UBS93.gene.df$ENSEMBL,
                                                             HER2.amp.signature.genes = UBS93.data$HER2.amp.signature.genes,
                                                             ssgsea.cutoff = UBS93.data$ssgsea.cutoff)
HDQP1_subtype                     <- HDQP1_cell_subtype$subtype
HDQP1_subtype$sample.id           <- rownames(HDQP1_subtype)
table(HDQP1_subtype$subtype)
GSE202771_HDQP1_metadata            <- HDQP1_scRNA@meta.data
GSE202771_HDQP1_metadata$sample.id  <- rownames(GSE202771_HDQP1_metadata)
GSE202771_HDQP1_metadata            <- merge(GSE202771_HDQP1_metadata,HDQP1_subtype,by="sample.id")
rown                                <- GSE202771_HDQP1_metadata$sample.id
GSE202771_HDQP1_metadata            <- GSE202771_HDQP1_metadata[,-1]
rownames(GSE202771_HDQP1_metadata)  <- rown
HDQP1_scRNA@meta.data               <- GSE202771_HDQP1_metadata
save(HDQP1_scRNA,file = "Pro_TNBC/paper/data/results/section_4/GSE202771_HDQP1_scRNA.RData")

####*Figure S11a: the t-SNE plot for HDQP1 data####
df1             <- as.data.frame(HDQP1_scRNA@reductions$tsne@cell.embeddings)
df1$cell.id     <- rownames(df1)
subtype         <- as.data.frame(HDQP1_scRNA@meta.data)
subtype$cell.id <- rownames(subtype)
df1             <- merge(df1,subtype,by="cell.id")
p8              <- ggplot(data = df1,aes(x=tSNE_1,y=tSNE_2,col=subtype))+
  geom_point(size=5)+xlab("")+ylab("")+
  ggplot.style+theme(plot.title = element_text(hjust=0.5,size = 44))+scale_color_manual(values = c("Red","Blue","Purple","#FFCC00"))
ggsave(p8,filename = "Pro_TNBC/paper/plot/section_4/HDQP1_GSE202771/tsne.plot.of.GSE202771.HDQP1.by.using.all.genes.pdf",width=23,height=15)

####*Figure S11b: a heatmap for the expression values of 5 marker genes in HDQP1 data ####
HDQP1_exp_df                <- HDQP1_exp
marker.gene                 <- c("CLDN3","CLDN4","CLDN7","VIM","EPCAM")
HDQP1_exp_df                <- HDQP1_exp_df[rownames(HDQP1_exp_df) %in% marker.gene,]
HDQP1_subtype               <- HDQP1_scRNA@meta.data
HDQP1_subtype$sample.id     <- rownames(HDQP1_subtype)
identical(HDQP1_subtype$sample.id,colnames(HDQP1_exp_df))
annotation.col              <- HDQP1_subtype[,1,drop=F]
rownames(annotation.col)    <- HDQP1_subtype$sample.id 
annot.colors                <- list(subtype=c("Basal"="Red","Claudin_low"= "Blue","HER2_amp"="Purple","Luminal"="#FFCC00"))
library(pheatmap)
fig4b          <- pheatmap(HDQP1_exp_df,cluster_rows = F,
                  cluster_cols = T,
                  color = colorRampPalette(c("navy","white","firebrick3"))(100),
                  show_colnames = F,border_color = NA,
                  scale = "row",
                  show_rownames = T,annotation_col = annotation.col,
                  fontsize=30,cellwidth = 0.1, cellheight = 55,
                  annotation_colors = annot.colors
                  
)
#fontsize_row = 8,fontsize_col=12 ,
ggsave(fig4b,filename = "Pro_TNBC/paper/plot/section_4/HDQP1_GSE202771/heat.map.of.marker.gene.in.HDQP1.GSE202771.pdf",width = 20,height = 15)
save(HDQP1_scRNA,file = "Pro_TNBC/paper/data/results/section_4/GSE202771_HDQP1_scRNA.RData")

####*Figure S12: box plots for the expression values of marker genes in every subtype####
HDQP1_exp_df              <- as.matrix(HDQP1_exp_df) %>% t() %>% as.data.frame()
HDQP1_exp_df$sample.id    <- rownames(HDQP1_exp_df)
HDQP1_exp_df              <- merge(HDQP1_exp_df,HDQP1_subtype,by="sample.id")
supfig_11a                <- ggplot(HDQP1_exp_df,aes(x=subtype,y=VIM))+
  geom_boxplot(outlier.shape=NA,lwd=1.2)+
  geom_point(size=4,position="jitter") + ggplot.style + 
  xlab('')+
  theme(panel.grid.major=element_line(colour=NA),
        panel.grid.minor = element_blank())+
  scale_x_discrete("", labels = c("Basal" = paste("Basal",sep = "\n","(n=6729)"),"Claudin_low" = paste("Cladin-low",sep = "\n","(n=239)"),"HER2_amp"=paste("HER2_amp",sep = "\n","(n=3)"),"Luminal"=paste("Luminal",sep = "\n","(n=1)")))
HDQP1_exp_bcl             <- subset(HDQP1_exp_df,subtype=="Basal"|subtype=="Claudin_low")
wilcox.test(VIM ~subtype,data = HDQP1_exp_bcl)#p-value < 2.2e-16
ggsave(supfig_11a,filename = "Pro_TNBC/paper/plot/section_4/HDQP1_GSE202771/boxplot.VIM.pdf",width = 20,height = 15)
supfig_11b                <- ggplot(HDQP1_exp_df,aes(x=subtype,y=EPCAM))+
  geom_boxplot(outlier.shape=NA,lwd=1.2)+
  geom_point(size=4,position="jitter") + ggplot.style + 
  xlab('')+
  theme(panel.grid.major=element_line(colour=NA),
        panel.grid.minor = element_blank())+
  scale_x_discrete("", labels = c("Basal" = paste("Basal",sep = "\n","(n=6729)"),"Claudin_low" = paste("Cladin-low",sep = "\n","(n=239)"),"HER2_amp"=paste("HER2_amp",sep = "\n","(n=3)"),"Luminal"=paste("Luminal",sep = "\n","(n=1)")))
wilcox.test(EPCAM ~subtype,data = HDQP1_exp_bcl)#p-value < 2.2e-16
ggsave(supfig_11b,filename = "Pro_TNBC/paper/plot/section_4/HDQP1_GSE202771/boxplot.EPCAM.pdf",width = 20,height = 15)
supfig_11c                 <- ggplot(HDQP1_exp_df,aes(x=subtype,y=CLDN3))+
  geom_boxplot(outlier.shape=NA,lwd=1.2)+
  geom_point(size=4,position="jitter") + ggplot.style + 
  xlab('')+
  theme(panel.grid.major=element_line(colour=NA),
        panel.grid.minor = element_blank())+
  scale_x_discrete("", labels = c("Basal" = paste("Basal",sep = "\n","(n=6729)"),"Claudin_low" = paste("Cladin-low",sep = "\n","(n=239)"),"HER2_amp"=paste("HER2_amp",sep = "\n","(n=3)"),"Luminal"=paste("Luminal",sep = "\n","(n=1)")))
wilcox.test(CLDN3 ~subtype,data = HDQP1_exp_bcl)#p-value < 2.2e-16
ggsave(supfig_11c,filename = "Pro_TNBC/paper/plot/section_4/HDQP1_GSE202771//boxplot.CLDN3.pdf",width = 20,height = 15)

supfig_11d                 <- ggplot(HDQP1_exp_df,aes(x=subtype,y=CLDN4))+
  geom_boxplot(outlier.shape=NA,lwd=1.2)+
  geom_point(size=4,position="jitter") + ggplot.style + 
  xlab('')+
  theme(panel.grid.major=element_line(colour=NA),
        panel.grid.minor = element_blank())+
  scale_x_discrete("", labels = c("Basal" = paste("Basal",sep = "\n","(n=6729)"),"Claudin_low" = paste("Cladin-low",sep = "\n","(n=239)"),"HER2_amp"=paste("HER2_amp",sep = "\n","(n=3)"),"Luminal"=paste("Luminal",sep = "\n","(n=1)")))
wilcox.test(CLDN4 ~subtype,data = HDQP1_exp_bcl)#p-value < 2.2e-16
ggsave(supfig_11d,filename = "Pro_TNBC/paper/plot/section_4/HDQP1_GSE202771/boxplot.CLDN4.pdf",width = 20,height = 15)

supfig_11e                 <- ggplot(HDQP1_exp_df,aes(x=subtype,y=CLDN7))+
  geom_boxplot(outlier.shape=NA,lwd=1.2)+
  geom_point(size=4,position="jitter") + ggplot.style + 
  xlab('')+
  theme(panel.grid.major=element_line(colour=NA),
        panel.grid.minor = element_blank())+
  scale_x_discrete("", labels = c("Basal" = paste("Basal",sep = "\n","(n=6729)"),"Claudin_low" = paste("Cladin-low",sep = "\n","(n=239)"),"HER2_amp"=paste("HER2_amp",sep = "\n","(n=3)"),"Luminal"=paste("Luminal",sep = "\n","(n=1)")))
wilcox.test(CLDN7 ~subtype,data = HDQP1_exp_bcl)#p-value < 2.2e-16
ggsave(supfig_11e,filename = "Pro_TNBC/paper/plot/section_4/HDQP1_GSE202771/boxplot.CLDN7.pdf",width = 20,height = 15)


#########################################################################################################
####Figure S13: correlation coefficients between Claudin-low centroid and Claudin-low cells in TNBC 
#########################################################################################################
load("Pro_TNBC/paper/data/results/section_4/TNBC/TNBC_CL_cor.RData")
fig_1b                           <- ggplot(TNBC_CL_cor,aes(x=reorder(patient.id,Claudin_low,FUN = median),y=Claudin_low))+
  geom_boxplot() + 
  theme(axis.text = element_text(size=12, face="bold",angle = 90),
        axis.title = element_text(size=45, face="bold"))+
  geom_point()+
  theme(panel.grid.major=element_line(colour=NA),
        panel.grid.minor = element_blank())+xlab("Patient id")+ylab("Correlation coefficients")+
  scale_x_discrete("Patient ID", labels = c(
    "BIOKEY_1" = paste("BIOKEY_1", sep = "\n", "(n=55)"),
    "BIOKEY_10" = paste("BIOKEY_10", sep = "\n", "(n=456)"),
    "BIOKEY_11" = paste("BIOKEY_11", sep = "\n", "(n=29)"),
    "BIOKEY_14" = paste("BIOKEY_14", sep = "\n", "(n=25)"),
    "BIOKEY_15" = paste("BIOKEY_15", sep = "\n", "(n=29)"),
    "BIOKEY_16" = paste("BIOKEY_16", sep = "\n", "(n=49)"),
    "BIOKEY_19" = paste("BIOKEY_19", sep = "\n", "(n=38)"),
    "BIOKEY_2" = paste("BIOKEY_2", sep = "\n", "(n=44)"),
    "BIOKEY_25" = paste("BIOKEY_25", sep = "\n", "(n=2)"),
    "BIOKEY_26" = paste("BIOKEY_26", sep = "\n", "(n=465)"),
    "BIOKEY_31" = paste("BIOKEY_31", sep = "\n", "(n=35)"),
    "BIOKEY_33" = paste("BIOKEY_33", sep = "\n", "(n=13)"),
    "BIOKEY_35" = paste("BIOKEY_35", sep = "\n", "(n=159)"),
    "BIOKEY_36" = paste("BIOKEY_36", sep = "\n", "(n=4)"),
    "BIOKEY_39" = paste("BIOKEY_39", sep = "\n", "(n=5)"),
    "BIOKEY_41" = paste("BIOKEY_41", sep = "\n", "(n=127)"),
    "BIOKEY_8" = paste("BIOKEY_8", sep = "\n", "(n=2)"),
    "BIOKEY_9" = paste("BIOKEY_9", sep = "\n", "(n=305)"),
    "CID4495" = paste("CID4495", sep = "\n", "(n=93)"),
    "CID44971" = paste("CID44971", sep = "\n", "(n=11)"),
    "CID44991" = paste("CID44991", sep = "\n", "(n=7)"),
    "CID4513" = paste("CID4513", sep = "\n", "(n=979)"),
    "CID4523" = paste("CID4523", sep = "\n", "(n=8)"),
    "PT039" = paste("PT039", sep = "\n", "(n=2)"),
    "PT084" = paste("PT084", sep = "\n", "(n=3)"),
    "PT089" = paste("PT089", sep = "\n", "(n=10)"),
    "PT126" = paste("PT126", sep = "\n", "(n=15)"),
    "sc5rJUQ033" = paste("sc5rJUQ033", sep = "\n", "(n=40)"),
    "sc5rJUQ039" = paste("sc5rJUQ039", sep = "\n", "(n=1)"),
    "sc5rJUQ042" = paste("sc5rJUQ042", sep = "\n", "(n=6)"),
    "sc5rJUQ045" = paste("sc5rJUQ045", sep = "\n", "(n=347)"),
    "sc5rJUQ053" = paste("sc5rJUQ053", sep = "\n", "(n=36)"),
    "SU4" = paste("SU4", sep = "\n", "(n=21)"),
    "SU58" = paste("SU58", sep = "\n", "(n=1)"),
    "TNBC_0126" = paste("TNBC_0126", sep = "\n", "(n=12)"),
    "TNBC_0131" = paste("TNBC_0131", sep = "\n", "(n=25)"),
    "TNBC_0135" = paste("TNBC_0135", sep = "\n", "(n=3)"),
    "TNBC_0554" = paste("TNBC_0554", sep = "\n", "(n=215)")
  ))
ggsave(fig_1b,filename = "Pro_TNBC/paper/plot/section_4/Box.plot.of.cor.cl.cells.and.centroid.pdf",height = 15,width = 20)


#####################################
####analyzing TNBC 0554 patient
#####################################

load("Pro_TNBC/output/data/scRNASeq/38sample/GSE161529_TNBC/GSE161529_TNBC_scRNA.RData")
TNBC_BRCA1_0554_scRNA                           <- subset(GSE161529_TNBC_scRNA,orig.ident=="TNBC_0554")
TNBC_BRCA1_0554_tumor_scRNA                     <- subset(TNBC_BRCA1_0554_scRNA,celltype=="Epithelial_cells"&cell=="tumorcell")
TNBC_BRCA1_0554_tumor_EXP                       <- GetAssayData(TNBC_BRCA1_0554_tumor_scRNA,assay = "RNA",slot = "count") 
TNBC_BRCA1_0554_tumor_scRNA                     <- CreateSeuratObject(counts =TNBC_BRCA1_0554_tumor_EXP )
#calculate the proportion of ribosomal genes in cells
TNBC_BRCA1_0554_tumor_scRNA[["percent.mt"]]     <- PercentageFeatureSet(TNBC_BRCA1_0554_tumor_scRNA,pattern = "^MT-")
head(TNBC_BRCA1_0554_tumor_scRNA@meta.data,5)
VlnPlot(TNBC_BRCA1_0554_tumor_scRNA,features = c("nCount_RNA","nFeature_RNA","percent.mt"),ncol=3)
TNBC_BRCA1_0554_tumor_scRNA                     <- subset(TNBC_BRCA1_0554_tumor_scRNA,nFeature_RNA >0 & nFeature_RNA<7500)
#normalized
TNBC_BRCA1_0554_tumor_scRNA                     <- NormalizeData(TNBC_BRCA1_0554_tumor_scRNA)
TNBC_BRCA1_0554_tumor_scRNA[["RNA"]]@data[c("ENSG00000167286", "ENSG00000100721", "ENSG00000156738"), 1:10]#data after normalizing
#features choosing
TNBC_BRCA1_0554_tumor_scRNA                     <- FindVariableFeatures(TNBC_BRCA1_0554_tumor_scRNA,selection.method = "vst",nfeatures = 2000)#choose high variability between cells.
#identify the 10 most highly variable genes
top10                                           <- head(VariableFeatures(TNBC_BRCA1_0554_tumor_scRNA),10)
#genes related to cell cycle
cellcycle_gene                                  <- c(cc.genes$s.genes,cc.genes$g2m.genes)
head(cellcycle_gene)
CaseMatch(cellcycle_gene,VariableFeatures(TNBC_BRCA1_0554_tumor_scRNA))#check which cell cycle-related genes are among the hypervariable genes we selected.
# results:character(0)
#data scaling 
all_gene                                        <- rownames(TNBC_BRCA1_0554_tumor_scRNA)
TNBC_BRCA1_0554_tumor_scRNA                     <- ScaleData(TNBC_BRCA1_0554_tumor_scRNA,features = all_gene)#scaled data存放在 pbmc[["RNA"]]@scale.data
#linear dimensionality reduction
TNBC_BRCA1_0554_tumor_scRNA                     <- RunPCA(TNBC_BRCA1_0554_tumor_scRNA,features = VariableFeatures(object=TNBC_BRCA1_0554_tumor_scRNA))
#dimension selection
TNBC_BRCA1_0554_tumor_scRNA                     <- JackStraw(TNBC_BRCA1_0554_tumor_scRNA,num.replicate = 100)
TNBC_BRCA1_0554_tumor_scRNA                     <- ScoreJackStraw(TNBC_BRCA1_0554_tumor_scRNA,dims = 1:20)
JackStrawPlot(TNBC_BRCA1_0554_tumor_scRNA,dims=1:20)
ElbowPlot(TNBC_BRCA1_0554_tumor_scRNA,ndims = 50)#choose dims of 30
#cell clustering
TNBC_BRCA1_0554_tumor_scRNA                     <- FindNeighbors(TNBC_BRCA1_0554_tumor_scRNA,dims = 1:30)
TNBC_BRCA1_0554_tumor_scRNA                     <- FindClusters(TNBC_BRCA1_0554_tumor_scRNA,resolution = 0.3)
table(Idents(TNBC_BRCA1_0554_tumor_scRNA))


####*TSNE####
TNBC_BRCA1_0554_tumor_exp                       <- as.matrix(TNBC_BRCA1_0554_tumor_scRNA@assays$RNA@data)
TNBC_BRCA1_0554_tumor_exp_UBS83                 <- TNBC_BRCA1_0554_tumor_exp[rownames(TNBC_BRCA1_0554_tumor_exp) %in% UBS93.gene.df$ENSEMBL,]
library(Rtsne)
cor_matrix                                      <- cor(TNBC_BRCA1_0554_tumor_exp_UBS83,method = "spearman")
dissimilarity_matrix                            <- 1 - cor_matrix
normalized_dissimilarity_matrix                 <- dissimilarity_matrix / max(dissimilarity_matrix)
TNBC_BRCA1_0554_tumor_scRNA                     <- RunTSNE(
  TNBC_BRCA1_0554_tumor_scRNA,
  reduction = "pca",
  cells = NULL,
  dims = 1:30,
  features = NULL,
  seed.use = 1,
  tsne.method = "Rtsne",
  dim.embed = 2,
  distance.matrix = normalized_dissimilarity_matrix,
  reduction.name = "tsne",
  reduction.key = "tSNE_"
)

####*run UBS93 predicting for TNBC 0554 data####
TNBC_BRCA1_0554_tumor_exp                         <- as.matrix(TNBC_BRCA1_0554_tumor_scRNA@assays$RNA@data) #normalied data
TNBC_BRCA1_0554_tumor_cell_subtype                <- breast.cancer.predictor(expr.of.sample = TNBC_BRCA1_0554_tumor_exp,
                                                                             expr.of.centroid = UBS93.data$UBS93.centroid,
                                                                             marker.gene = UBS93.data$UBS93.gene.df$ENSEMBL,
                                                                             HER2.amp.signature.genes = UBS93.data$HER2.amp.signature.genes,
                                                                             ssgsea.cutoff = UBS93.data$ssgsea.cutoff)
TNBC_BRCA1_0554_tumor_subtype                     <- TNBC_BRCA1_0554_tumor_cell_subtype$subtype
TNBC_BRCA1_0554_tumor_subtype$sample.id           <- rownames(TNBC_BRCA1_0554_tumor_subtype)
TNBC_BRCA1_0554_tumor_scRNA_metadata              <- TNBC_BRCA1_0554_tumor_scRNA@meta.data
TNBC_BRCA1_0554_tumor_scRNA_metadata$sample.id    <- rownames(TNBC_BRCA1_0554_tumor_scRNA_metadata)
TNBC_BRCA1_0554_tumor_scRNA_metadata              <- merge(TNBC_BRCA1_0554_tumor_scRNA_metadata,TNBC_BRCA1_0554_tumor_subtype,by="sample.id")
rown                                              <- TNBC_BRCA1_0554_tumor_scRNA_metadata$sample.id
TNBC_BRCA1_0554_tumor_scRNA_metadata              <- TNBC_BRCA1_0554_tumor_scRNA_metadata[,-1]
rownames(TNBC_BRCA1_0554_tumor_scRNA_metadata)    <- rown
TNBC_BRCA1_0554_tumor_scRNA@meta.data             <- TNBC_BRCA1_0554_tumor_scRNA_metadata
save(TNBC_BRCA1_0554_tumor_scRNA,file = "Pro_TNBC/output/data/scRNASeq/38sample/TNBC_BRCA1_0554/TNBC_BRCA1_0554_tumor_scRNA.RData")

####*Figure 4c: the t-SNE plot for TNBC 0554 data####
df1             <- as.data.frame(TNBC_BRCA1_0554_tumor_scRNA@reductions$tsne@cell.embeddings)
df1$cell.id     <- rownames(df1)
subtype         <- as.data.frame(TNBC_BRCA1_0554_tumor_scRNA@meta.data)
subtype$cell.id <- rownames(subtype)
df1             <- merge(df1,subtype,by="cell.id")
fig_4c              <- ggplot(data = df1,aes(x=tSNE_1,y=tSNE_2,col=subtype))+
  geom_point(size=5)+xlab("")+ylab("")+
  ggplot.style+theme(plot.title = element_text(hjust=0.5,size = 44))+scale_color_manual(values = c("Red","Blue","Purple","#FFCC00"))
ggsave(fig_4c,filename = "Pro_TNBC/paper/plot/section_4/TNBC_BRCA1_0554/tsne.plot.of.TNBC_BRCA1_0554_tumor.by.using.all.genes.pdf",width=23,height=15)



####*Figure S14: a T-SNE visualization of these 5 marker genes from TNBC patient 0554.####
library(ggplot2)
p1 <- FeaturePlot(TNBC_BRCA1_0554_tumor_scRNA,features = c("ENSG00000165215"),#CLDN3
                  reduction = "tsne",pt.size = 5)+
  theme( 
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
    axis.ticks = element_blank(),axis.text.x = element_blank(),axis.title.y = element_blank(),
    legend.position = "none", 
    plot.title = element_text(size=55,hjust=0.5))+ggtitle("CLDN3")+theme_bw(base_size = 55) + theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))+xlab("")+ylab("")
ggsave(p1,filename = "Pro_TNBC/paper/plot/section_4/CLDN3.in.TNBC_BRCA1_0554_tumor.of.GSE202771.pdf",width=20,height=15)

p2 <- FeaturePlot(TNBC_BRCA1_0554_tumor_scRNA,features = c("ENSG00000189143"),#CLDN3
                  reduction = "tsne",pt.size = 5)+
  theme( 
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
    axis.ticks = element_blank(),axis.text.x = element_blank(),axis.title.y = element_blank(),
    legend.position = "none", 
    plot.title = element_text(size=55,hjust=0.5),
    plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))+ggtitle("CLDN4")+theme_bw(base_size = 55) +
  xlab("")+ylab("")
ggsave(p2,filename = "Pro_TNBC/paper/plot/section_4/CLDN4.in.TNBC_BRCA1_0554_tumor.of.GSE202771.pdf",width=20,height=15)

p3 <- FeaturePlot(TNBC_BRCA1_0554_tumor_scRNA,features = c("ENSG00000181885"),#CLDN3
                  reduction = "tsne",pt.size = 5)+
  theme( 
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
    axis.ticks = element_blank(),axis.text.x = element_blank(),axis.title.y = element_blank(),
    legend.position = "none", 
    plot.title = element_text(size=55,hjust=0.5),
    plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))+ggtitle("CLDN7")+theme_bw(base_size = 55) +
  xlab("")+ylab("")
ggsave(p3,filename = "Pro_TNBC/paper/plot/section_4/CLDN7.in.TNBC_BRCA1_0554_tumor.of.GSE202771.pdf",width=20,height=15)

p4 <- FeaturePlot(TNBC_BRCA1_0554_tumor_scRNA,features = c("ENSG00000026025"),#CLDN3
                  reduction = "tsne",pt.size = 5)+
  theme( 
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
    axis.ticks = element_blank(),axis.text.x = element_blank(),axis.title.y = element_blank(),
    legend.position = "none", 
    plot.title = element_text(size=55,hjust=0.5),
    plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))+ggtitle("VIM")+theme_bw(base_size = 55) +
  xlab("")+ylab("")
ggsave(p4,filename = "Pro_TNBC/paper/plot/section_4/VIM.in.TNBC_BRCA1_0554_tumor.of.GSE202771.pdf",width=20,height=15)

p5 <- FeaturePlot(TNBC_BRCA1_0554_tumor_scRNA,features = c("ENSG00000119888"),#CLDN3
                  reduction = "tsne",pt.size = 5)+
  theme( 
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
    axis.ticks = element_blank(),axis.text.x = element_blank(),axis.title.y = element_blank(),
    legend.position = "none", 
    plot.title = element_text(size=55,hjust=0.5),
    plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))+ggtitle("EPCAM")+theme_bw(base_size = 55) +
  xlab("")+ylab("")
ggsave(p5,filename = "Pro_TNBC/paper/plot/section_4/EPCAM.in.TNBC_BRCA1_0554_tumor.of.GSE202771.pdf",width=20,height=15)
save(TNBC_BRCA1_0554_tumor_scRNA,file = "Pro_TNBC/paper/data/section_4/TNBC_BRCA1_0554_tumor_scRNA(tsne_by_cor).RData")


####*Figure 4d: a heatmap for the expression values of 5 marker genes in TNBC 0554 data ####
TNBC_BRCA1_0554_tumor_exp_df          <- as.data.frame(TNBC_BRCA1_0554_tumor_exp)
TNBC_BRCA1_0554_tumor_exp_df$ENSEMBL  <- rownames(TNBC_BRCA1_0554_tumor_exp_df) 
symbol                                <- bitr(TNBC_BRCA1_0554_tumor_exp_df$ENSEMBL,fromType = "ENSEMBL",toType = "SYMBOL",OrgDb="org.Hs.eg.db")
TNBC_BRCA1_0554_tumor_exp_df          <- merge(symbol,TNBC_BRCA1_0554_tumor_exp_df,by="ENSEMBL")
marker.gene                           <- c("CLDN3","CLDN4","CLDN7","VIM","EPCAM")
TNBC_BRCA1_0554_tumor_exp_df          <- TNBC_BRCA1_0554_tumor_exp_df[TNBC_BRCA1_0554_tumor_exp_df$SYMBOL %in% marker.gene,]
rown                                  <- TNBC_BRCA1_0554_tumor_exp_df$SYMBOL
TNBC_BRCA1_0554_tumor_exp_df          <- TNBC_BRCA1_0554_tumor_exp_df[,-c(1,2)]
rownames(TNBC_BRCA1_0554_tumor_exp_df)      <- rown
TNBC_BRCA1_0554_tumor_exp_matrix            <- as.matrix(TNBC_BRCA1_0554_tumor_exp_df)
TNBC_BRCA1_0554_tumor_subtype               <- TNBC_BRCA1_0554_tumor_scRNA@meta.data
TNBC_BRCA1_0554_tumor_subtype$sample.id     <- rownames(TNBC_BRCA1_0554_tumor_subtype)
identical(TNBC_BRCA1_0554_tumor_subtype$sample.id,colnames(TNBC_BRCA1_0554_tumor_exp_matrix))
annotation.col              <- TNBC_BRCA1_0554_tumor_subtype[,7,drop=F]
rownames(annotation.col)    <- TNBC_BRCA1_0554_tumor_subtype$sample.id 
annot.colors                <- list(subtype=c("Basal"="Red","Claudin_low"= "Blue","Her2_amp"= "Purple","Luminal"= "#FFCC00"))
library(pheatmap)
fig4d            <- pheatmap(TNBC_BRCA1_0554_tumor_exp_matrix,cluster_rows = F,
                   cluster_cols = T,
                   color = colorRampPalette(c("navy","white","firebrick3"))(100),
                   show_colnames = F,border_color = NA,
                   scale = "row",
                   show_rownames = T,annotation_col = annotation.col,
                   fontsize=30,cellwidth = 0.4, cellheight = 55,
                   annotation_colors = annot.colors
)
#fontsize_row = 8,fontsize_col=12 ,
ggsave(fig4d,filename = "Pro_TNBC/paper/plot/section_4/TNBC_BRCA1_0554/heat.map.of.marker.gene.in.TNBC_BRCA1_0554_tumor.pdf",width = 20,height = 15)

####*Figure S15: box plots for the expression values of marker genes in every subtype####
TNBC_BRCA1_0554_tumor_exp_df           <- as.matrix(TNBC_BRCA1_0554_tumor_exp_df) %>% t() %>% as.data.frame()
TNBC_BRCA1_0554_tumor_exp_df$sample.id <- rownames(TNBC_BRCA1_0554_tumor_exp_df)
TNBC_BRCA1_0554_tumor_exp_df           <- merge(TNBC_BRCA1_0554_tumor_exp_df,TNBC_BRCA1_0554_tumor_subtype,by="sample.id")
supfig_11a                             <- ggplot(TNBC_BRCA1_0554_tumor_exp_df,aes(x=subtype,y=VIM))+
  geom_boxplot(outlier.shape=NA,lwd=1.2)+
  geom_point(size=1,position="jitter") + ggplot.style + 
  xlab('')+
  theme(panel.grid.major=element_line(colour=NA),
        panel.grid.minor = element_blank())+
  scale_x_discrete("", labels = c("Basal" = paste("Basal",sep = "\n","(n=1961)"),"Claudin_low" = paste("Cladin-low",sep = "\n","(n=170)"),"Her2_amp"=paste("HER2-amp",sep = "\n","(n=18)"),"Luminal"=paste("Luminal",sep = "\n","(n=4)")))
kruskal.test(VIM ~subtype,data = TNBC_BRCA1_0554_tumor_exp_df)#p-value = 1.83e-14
ggsave(supfig_11a,filename = "Pro_TNBC/paper/plot/section_4/TNBC_BRCA1_0554/TNBC_BRCA1_0554_tumor_boxplot.VIM.pdf",width = 20,height = 15)
supfig_11b                             <- ggplot(TNBC_BRCA1_0554_tumor_exp_df,aes(x=subtype,y=EPCAM))+
  geom_boxplot(outlier.shape=NA,lwd=1.2)+
  geom_point(size=1,position="jitter") + ggplot.style + 
  xlab('')+
  theme(panel.grid.major=element_line(colour=NA),
        panel.grid.minor = element_blank())+
  scale_x_discrete("", labels = c("Basal" = paste("Basal",sep = "\n","(n=1961)"),"Claudin_low" = paste("Cladin-low",sep = "\n","(n=170)"),"Her2_amp"=paste("HER2-amp",sep = "\n","(n=18)"),"Luminal"=paste("Luminal",sep = "\n","(n=4)")))
kruskal.test(EPCAM ~subtype,data = TNBC_BRCA1_0554_tumor_exp_df)#p-value = 4.123e-13
ggsave(supfig_11b,filename = "Pro_TNBC/paper/plot/section_4/TNBC_BRCA1_0554/TNBC_BRCA1_0554_tumor_boxplot.EPCAM.pdf",width = 20,height = 15)
supfig_11c                             <- ggplot(TNBC_BRCA1_0554_tumor_exp_df,aes(x=subtype,y=CLDN3))+
  geom_boxplot(outlier.shape=NA,lwd=1.2)+
  geom_point(size=1,position="jitter") + ggplot.style + 
  xlab('')+
  theme(panel.grid.major=element_line(colour=NA),
        panel.grid.minor = element_blank())+
  scale_x_discrete("", labels = c("Basal" = paste("Basal",sep = "\n","(n=1961)"),"Claudin_low" = paste("Cladin-low",sep = "\n","(n=170)"),"Her2_amp"=paste("HER2-amp",sep = "\n","(n=18)"),"Luminal"=paste("Luminal",sep = "\n","(n=4)")))
kruskal.test(CLDN3 ~subtype,data = TNBC_BRCA1_0554_tumor_exp_df)#p-value = 9.982e-05
ggsave(supfig_11c,filename = "Pro_TNBC/paper/plot/section_4/TNBC_BRCA1_0554/TNBC_BRCA1_0554_tumor_boxplot.CLDN3.pdf",width = 20,height = 15)

supfig_11d                             <- ggplot(TNBC_BRCA1_0554_tumor_exp_df,aes(x=subtype,y=CLDN4))+
  geom_boxplot(outlier.shape=NA,lwd=1.2)+
  geom_point(size=1,position="jitter") + ggplot.style + 
  xlab('')+
  theme(panel.grid.major=element_line(colour=NA),
        panel.grid.minor = element_blank())+
  scale_x_discrete("", labels = c("Basal" = paste("Basal",sep = "\n","(n=1961)"),"Claudin_low" = paste("Cladin-low",sep = "\n","(n=170)"),"Her2_amp"=paste("HER2-amp",sep = "\n","(n=18)"),"Luminal"=paste("Luminal",sep = "\n","(n=4)")))
kruskal.test(CLDN4 ~subtype,data = TNBC_BRCA1_0554_tumor_exp_df)#p-value = 5.771e-09
ggsave(supfig_11d,filename = "Pro_TNBC/paper/plot/section_4/TNBC_BRCA1_0554/TNBC_BRCA1_0554_tumor_boxplot.CLDN4.pdf",width = 20,height = 15)

supfig_11e                             <- ggplot(TNBC_BRCA1_0554_tumor_exp_df,aes(x=subtype,y=CLDN7))+
  geom_boxplot(outlier.shape=NA,lwd=1.2)+
  geom_point(size=1,position="jitter") + ggplot.style + 
  xlab('')+
  theme(panel.grid.major=element_line(colour=NA),
        panel.grid.minor = element_blank())+
  scale_x_discrete("", labels = c("Basal" = paste("Basal",sep = "\n","(n=1961)"),"Claudin_low" = paste("Cladin-low",sep = "\n","(n=170)"),"Her2_amp"=paste("HER2-amp",sep = "\n","(n=18)"),"Luminal"=paste("Luminal",sep = "\n","(n=4)")))
kruskal.test(CLDN7 ~subtype,data = TNBC_BRCA1_0554_tumor_exp_df)#p-value = 6.104e-16
ggsave(supfig_11e,filename = "Pro_TNBC/paper/plot/section_4/TNBC_BRCA1_0554/TNBC_BRCA1_0554_tumor_boxplot.CLDN7.pdf",width = 20,height = 15)


####*infercnv of all cell types(dividing tumor cells into 4 subtypes by using UBS93) ####
TNBC_BRCA1_0554_metadata             <- TNBC_BRCA1_0554_scRNA@meta.data
TNBC_BRCA1_0554_metadata$sample.id   <- rownames(TNBC_BRCA1_0554_metadata)
basal.id                             <- TNBC_BRCA1_0554_tumor_subtype[TNBC_BRCA1_0554_tumor_subtype$subtype=="Basal",]$sample.id
claudin_low_id                       <- TNBC_BRCA1_0554_tumor_subtype[TNBC_BRCA1_0554_tumor_subtype$subtype=="Claudin_low",]$sample.id                                                                 
her2_id                              <- TNBC_BRCA1_0554_tumor_subtype[TNBC_BRCA1_0554_tumor_subtype$subtype=="Her_amp",]$sample.id
Luminal_id                           <- TNBC_BRCA1_0554_tumor_subtype[TNBC_BRCA1_0554_tumor_subtype$subtype=="Luminal",]$sample.id
B_cell_id                            <- TNBC_BRCA1_0554_metadata[TNBC_BRCA1_0554_metadata$celltype=="B_cell",]$sample.id
CMP_id                               <- TNBC_BRCA1_0554_metadata[TNBC_BRCA1_0554_metadata$celltype=="CMP",]$sample.id
DC_id                                <- TNBC_BRCA1_0554_metadata[TNBC_BRCA1_0554_metadata$celltype=="DC",]$sample.id
Monocyte_id                          <- TNBC_BRCA1_0554_metadata[TNBC_BRCA1_0554_metadata$celltype=="Monocyte",]$sample.id
T_cells_id                           <- TNBC_BRCA1_0554_metadata[TNBC_BRCA1_0554_metadata$celltype=="T_cells",]$sample.id
fibroblasts_id                       <- TNBC_BRCA1_0554_metadata[TNBC_BRCA1_0554_metadata$celltype=="Fibroblasts",]$sample.id
normal_epithlial_id                  <- TNBC_BRCA1_0554_metadata[TNBC_BRCA1_0554_metadata$celltype=="Epithelial_cells"&TNBC_BRCA1_0554_metadata$cell=="normalcell",]$sample.id
Pericytes_id                         <- TNBC_BRCA1_0554_metadata[TNBC_BRCA1_0554_metadata$celltype=="Pericytes",]$sample.id
TNBC_BRCA1_0554_metadata             <- within(TNBC_BRCA1_0554_metadata,{
  cell_subtype                                 <- NA
  cell_subtype[sample.id %in% basal.id]        <- "Basal"
  cell_subtype[sample.id %in% claudin_low_id]  <- "Claudin_low"
  cell_subtype[sample.id %in% her2_id]         <- "HER2_amp"
  cell_subtype[sample.id %in% Luminal_id]      <- "Luminal"
  cell_subtype[sample.id %in% B_cell_id]       <- "B_cell"
  cell_subtype[sample.id %in% CMP_id]          <- "CMP"
  cell_subtype[sample.id %in% DC_id]           <- "DC"
  cell_subtype[sample.id %in% Monocyte_id]     <- "Monocyte"
  cell_subtype[sample.id %in% fibroblasts_id]  <- "Fibroblasts"
  cell_subtype[sample.id %in% T_cells_id]      <- "T_cells"
  cell_subtype[sample.id %in% normal_epithlial_id]      <- "normal_epithlial"
  cell_subtype[sample.id %in% Pericytes_id]             <- "Pericytes"
  
})
TNBC_BRCA1_0554_scRNA@meta.data                <- TNBC_BRCA1_0554_metadata
save(TNBC_BRCA1_0554_scRNA,file = "Pro_TNBC/output/data/scRNASeq/38sample/TNBC_BRCA1_0554/TNBC_BRCA1_0554_scRNA.RData")


library(Seurat)
library(tidyverse)
library(ggplot2)
library(infercnv)
library(ComplexHeatmap)
library(ggpubr)
library(HiddenMarkov)

#inferCNV needs three files:1.count expression matrix,2.group information,3.Gene chromosome information
#make gene chromosome position information and extract expression matrix.
dat             <- GetAssayData(TNBC_BRCA1_0554_scRNA,assay = "RNA",slot = "count") 
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
meta                <- subset(TNBC_BRCA1_0554_scRNA@meta.data,select = c("cell_subtype"))
#inferCNV
#Two-step construction object
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=dat,
                                    annotations_file=meta,
                                    delim="\t",
                                    gene_order_file=geneInfor,
                                    ref_group_names = c("T_cells"))   #Select the basic cell or sample to see the input type of meta, or do not select the algorithm to calculate by yourself according to the average value.
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir="Pro_TNBC/output/data/scRNASeq/38sample/TNBC_BRCA1_0554/infercnv_allcells_new/", 
                             cluster_by_groups=T,  # If TRUE is selected, group by sample will be changed to FALSE, and another parameter k will be selected_ obs_ The number of groups given by groups (the default is 1) for grouping.
                             denoise=T,     
                             HMM=F)   # #Whether to predict CNV based on HMM? True is a long time
#Finally, many files will be output in out_dir. You can directly use the heat map inside.
infercnv::plot_cnv(infercnv_obj, #上两步得到的infercnv对象
                   output_filename = "Pro_TNBC/output/data/scRNASeq/38sample/TNBC_BRCA1_0554/infercnv_allcells_new/infercnv_heatmap",output_format = "pdf", #保存为pdf文件
                   cluster_by_groups = T
) 
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir="Pro_TNBC/output/data/scRNASeq/38sample/TNBC_BRCA1_0554/infercnv_basalandCL_1", 
                             cluster_by_groups=F,  # If TRUE is selected, group by sample will be changed to FALSE, and another parameter k will be selected_ obs_ The number of groups given by groups (the default is 1) for grouping.
                             denoise=T,
                             HMM=F
) 


####*Figure 4e: a heatmap for showing the average cnv in every cell type####
load("~/Pro_TNBC/output/data/scRNASeq/38sample/TNBC_BRCA1_0554/TNBC_BRCA1_0554_scRNA.RData")
run.final.infercnv_obj = readRDS("Pro_TNBC/output/data/scRNASeq/38sample/TNBC_BRCA1_0554/infercnv_basalandCL/run.final.infercnv_obj")
expr              <- run.final.infercnv_obj@expr.data
dim(expr)
expr[1:4,1:4]
gn                <- rownames(expr)
library(AnnoProbe)  
geneInfor=annoGene(gn,"ENSEMBL",'human')  
colnames(geneInfor)
geneInfor=geneInfor[with(geneInfor, order(chr, start)),c(3,4:6)]      
geneInfor=geneInfor[!duplicated(geneInfor[,1]),]
length(unique(geneInfor[,1]))
geneFile          <- geneInfor
geneFile$chr      <- factor(geneFile$chr,levels = paste0("chr", seq(1, 22)))
geneFile          <- arrange(geneFile,chr)
identical(rownames(expr),geneFile$ENSEMBL)
expr              <- expr[match(geneFile$ENSEMBL,rownames(expr)),]
geneFile$new_cluster    <- sapply(as.character(geneFile$chr), function(x) substr(x, 4, nchar(x)))
geneFile$new_cluster    <- factor(geneFile$new_cluster,levels = seq(1, 22))
library(ComplexHeatmap)
metadata          <- TNBC_BRCA1_0554_scRNA@meta.data[,11,drop=F]
identical(colnames(expr),rownames(metadata))
metadata          <- metadata[match(colnames(expr),rownames(metadata)),,drop=F]
library(dplyr)
library(tibble)
expr_df            <- as.data.frame(t(expr))
identical(rownames(expr_df),rownames(metadata))
expr_df            <- cbind(expr_df,metadata)
expr_df[,1:6281]   <- apply(expr_df[,1:6281],2,as.numeric)
averages           <- expr_df %>% group_by(cell_subtype) %>% summarise_all(mean)
averages           <- column_to_rownames(averages,"cell_subtype")
averages                 <- t(averages)
identical(geneFile$ENSEMBL,rownames(averages))
new_cluster              <- geneFile$new_cluster
top_labels               <- HeatmapAnnotation(
  cluster = anno_block( labels = levels(new_cluster), 
                       labels_gp = gpar(cex = 0.9, col = "black"))) 
pdf("Pro_TNBC/paper/plot/section_4/TNBC_BRCA1_0554/heatmap.of.cnv.in.TNBC.0554.pdf",width = 15,height = 10)
averages_heatmap  <- Heatmap(t(averages),
                             cluster_rows = T,
                             cluster_columns = F,
                             show_column_names = F,
                             show_row_names = T,
                             column_split = new_cluster,
                             show_heatmap_legend=T,
                             top_annotation = top_labels,
                             column_title = "Genomic Region",
                             column_title_side = c("bottom"),
                             column_title_gp = gpar(fontsize = 25),
                             heatmap_legend_param = list(
                               title = "Modified Expression",
                               title_position = "leftcenter-rot", 
                               at=c(0.8,1.2), 
                               legend_height = unit(3, "cm") 
                             ))
draw(averages_heatmap, heatmap_legend_side = "right") # 图例位置
dev.off()
range(averages)
                                  
