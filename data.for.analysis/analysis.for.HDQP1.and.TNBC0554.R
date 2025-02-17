

##########################################
####analyzing HDQP1 in GSE173634####
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

####################################
####analyzing HDQP1 in GSE202771####
####################################

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
#unlinear dimension reduction(UMAP or tSNE)
HDQP1_scRNA                     <- RunUMAP(HDQP1_scRNA,dims = 1:40)
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


#####################################
####analyzing TNBC 0554 patient######
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
#unlinear dimension reduction(UMAP or tSNE)
TNBC_BRCA1_0554_tumor_scRNA                     <- RunUMAP(TNBC_BRCA1_0554_tumor_scRNA,dims = 1:30)

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
