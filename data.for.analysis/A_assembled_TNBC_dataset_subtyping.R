######################################################################################################
#### calculate correlation coefficients between Claudin-low centroid and Claudin-low cells in TNBC
######################################################################################################

library(GSVA)
load("Pro_TNBC/paper/data/method/UBS93.data.RData")
####GSE176078####
library(readr)
library(Seurat)
load("~/Pro_TNBC/output/data/scRNASeq/26_sample/GSE176078_scRNA.RData")
GSE176078_metadata                <- read_csv("Pro_TNBC/data/scRNASeq/26/Wu_etal_2021_BRCA_scRNASeq/metadata.csv")
GSE176078_tumorcell               <- GSE176078_metadata[GSE176078_metadata$celltype_major=="Cancer Epithelial",]
colnames(GSE176078_tumorcell)[1]  <- "cell.id"
sample.id                         <-  c("CID4465","CID4495","CID44971","CID44991","CID4513","CID4515","CID4523")
GSE176078_tumorcell_tnbc          <-  subset(GSE176078_tumorcell,orig.ident %in% sample.id)
GSE176078_exp                     <-  as.matrix(scRNA@assays$RNA@data)
GSE176078_tnbc_tumorcell_exp      <-  GSE176078_exp[,colnames(GSE176078_exp) %in% GSE176078_tumorcell_tnbc$cell.id]
exp_ubs93                         <- as.data.frame(GSE176078_tnbc_tumorcell_exp)
exp_ubs93$SYMBOL                  <- rownames(exp_ubs93)
exp_ubs93                         <- merge(UBS93.data$UBS93.gene.df,exp_ubs93,by="SYMBOL")
rown                              <- exp_ubs93$ENSEMBL
exp_ubs93                         <- exp_ubs93[,-c(1:3)] %>% as.matrix()
rownames(exp_ubs93)               <- rown
predictor.res                     <- breast.cancer.predictor(expr.of.sample = exp_ubs93,
                                                expr.of.centroid = UBS93.data$UBS93.centroid,
                                                marker.gene = UBS93.data$UBS93.gene.df$ENSEMBL,
                                                HER2.amp.signature.genes = UBS93.data$HER2.amp.signature.genes,
                                                ssgsea.cutoff = UBS93.data$ssgsea.cutoff)
GSE176078_TNBC_tumorcell_subtype                   <- predictor.res$subtype
table(GSE176078_TNBC_tumorcell_subtype$subtype)
GSE176078_TNBC_tumorcell_subtype$cell.id           <- rownames(GSE176078_TNBC_tumorcell_subtype)
GSE176078_TNBC_tumorcell_CL                        <- GSE176078_TNBC_tumorcell_subtype[GSE176078_TNBC_tumorcell_subtype=="Claudin_low",]$cell.id
GSE176078_TNBC_tumorcell_cormat                    <- predictor.res$cor.matrix  %>% as.data.frame()
GSE176078_TNBC_tumorcell_cormat$cell.id            <- rownames(GSE176078_TNBC_tumorcell_cormat)
GSE176078_TNBC_tumorcell_CLcormat                  <- subset(GSE176078_TNBC_tumorcell_cormat,cell.id %in% GSE176078_TNBC_tumorcell_CL)
GSE176078_TNBC_tumorcell_CLcormat                  <- merge(GSE176078_TNBC_tumorcell_CLcormat,GSE176078_tumorcell_tnbc,by="cell.id")
table(GSE176078_TNBC_tumorcell_CLcormat$orig.ident)
GSE176078_TNBC_CLcormat                            <- GSE176078_TNBC_tumorcell_CLcormat[,c(1,3,5)]
rown                                               <- GSE176078_TNBC_CLcormat$cell.id
GSE176078_TNBC_CLcormat                            <- GSE176078_TNBC_CLcormat[,-1]
rownames(GSE176078_TNBC_CLcormat)                  <- rown
colnames(GSE176078_TNBC_CLcormat)[2]               <- "patient.id"


####GSE161529####
load("Pro_TNBC/output/data/scRNASeq/38sample/GSE161529_TNBC/GSE161529_TNBC_scRNA.RData")
GSE161529_TNBC_metadata                            <- GSE161529_TNBC_scRNA@meta.data
GSE161529_TNBC_tumorcell                           <- subset(GSE161529_TNBC_metadata,celltype=="Epthelial_cells"|cell=="tumorcell")
table(GSE161529_TNBC_tumorcell$orig.ident)
GSE161529_TNBC_tumorcell_scRNA                     <- subset(GSE161529_TNBC_scRNA,celltype=="Epthelial_cells"|cell=="tumorcell")
GSE161529_TNBC_tumorcell_exprmat                   <- as.matrix(GSE161529_TNBC_tumorcell_scRNA@assays$RNA@data) # data after normalizing
predictor.res                                      <- breast.cancer.predictor(expr.of.sample = GSE161529_TNBC_tumorcell_exprmat,
                                                     expr.of.centroid = UBS93.data$UBS93.centroid,
                                                     marker.gene = UBS93.data$UBS93.gene.df$ENSEMBL,
                                                     HER2.amp.signature.genes = UBS93.data$HER2.amp.signature.genes,
                                                     ssgsea.cutoff = UBS93.data$ssgsea.cutoff)

GSE161529_TNBC_tumorcell_subtype                   <- predictor.res$subtype
table(GSE161529_TNBC_tumorcell_subtype$subtype)
GSE161529_TNBC_tumorcell_subtype$cell.id           <- rownames(GSE161529_TNBC_tumorcell_subtype)
GSE161529_TNBC_tumorcell_CL                        <- GSE161529_TNBC_tumorcell_subtype[GSE161529_TNBC_tumorcell_subtype=="Claudin_low",]$cell.id
GSE161529_TNBC_tumorcell_cormat                    <- predictor.res$cor.matrix  %>% as.data.frame()
GSE161529_TNBC_tumorcell_cormat$cell.id            <- rownames(GSE161529_TNBC_tumorcell_cormat)
GSE161529_TNBC_tumorcell_CLcormat                  <- subset(GSE161529_TNBC_tumorcell_cormat,cell.id %in% GSE161529_TNBC_tumorcell_CL)
GSE161529_TNBC_tumorcell$cell.id                   <- rownames(GSE161529_TNBC_tumorcell)
GSE161529_TNBC_tumorcell_CLcormat                  <- merge(GSE161529_TNBC_tumorcell_CLcormat,GSE161529_TNBC_tumorcell,by="cell.id")
table(GSE161529_TNBC_tumorcell_CLcormat$orig.ident)

GSE161529_TNBC_CLcormat                            <- GSE161529_TNBC_tumorcell_CLcormat[,c(1,3,5)]
rown                                               <- GSE161529_TNBC_CLcormat$cell.id
GSE161529_TNBC_CLcormat                            <- GSE161529_TNBC_CLcormat[,-1]
rownames(GSE161529_TNBC_CLcormat)                  <- rown
colnames(GSE161529_TNBC_CLcormat)[2]               <- "patient.id"

####GSE148673:single cell RNAseq data of patients####
####*TNBC1####
library(readr)
TNBC1_matrix                    <- read.delim("Pro_TNBC/data/scRNASeq/Data_Gao2021_Breast/GSE148673_RAW/GSM4476486_combined_UMIcount_CellTypes_TNBC1.txt",sep = "\t")
TNBC1_matrix                    <- as.matrix(TNBC1_matrix)  %>% t() %>% as.data.frame()
table(TNBC1_matrix$copykat.pred)
TNBC1_tumor_matrix              <- TNBC1_matrix[TNBC1_matrix$copykat.pred=="T",]
TNBC1_tumor_matrix              <- TNBC1_tumor_matrix[,-c(1,2)]  %>% as.matrix() %>% t()
TNBC1_scRNA                     <- CreateSeuratObject(counts = TNBC1_tumor_matrix)
#calculate the proportion of ribosomal genes in cells
TNBC1_scRNA[["percent.mt"]]     <- PercentageFeatureSet(TNBC1_scRNA,pattern = "^MT-")
head(TNBC1_scRNA@meta.data,5)
VlnPlot(TNBC1_scRNA,features = c("nCount_RNA","nFeature_RNA","percent.mt"),ncol=3)
TNBC1_scRNA                     <- subset(TNBC1_scRNA,nFeature_RNA >0 & nFeature_RNA<7500)
#normalized
TNBC1_scRNA                     <- NormalizeData(TNBC1_scRNA)
#features choosing
TNBC1_scRNA                     <- FindVariableFeatures(TNBC1_scRNA,selection.method = "vst",nfeatures = 2000)#choose high variability between cells.
#identify the 10 most highly variable genes
top10                           <- head(VariableFeatures(TNBC1_scRNA),10)
#data scaling 
all_gene                        <- rownames(TNBC1_scRNA)
TNBC1_scRNA                     <- ScaleData(TNBC1_scRNA,features = all_gene)#scaled data存放在 pbmc[["RNA"]]@scale.data
#linear dimensionality reduction
TNBC1_scRNA                     <- RunPCA(TNBC1_scRNA,features = VariableFeatures(object=TNBC1_scRNA))
#dimension selection
TNBC1_scRNA                     <- JackStraw(TNBC1_scRNA,num.replicate = 100)
TNBC1_scRNA                     <- ScoreJackStraw(TNBC1_scRNA,dims = 1:20)
JackStrawPlot(TNBC1_scRNA,dims=1:20)
ElbowPlot(TNBC1_scRNA,ndims = 50)#choose dims of 30
#cell clustering
TNBC1_scRNA                     <- FindNeighbors(TNBC1_scRNA,dims = 1:30)
TNBC1_scRNA                     <- FindClusters(TNBC1_scRNA,resolution = 0.1)
table(Idents(TNBC1_scRNA))
#unlinear dimension reduction(UMAP or tSNE)
TNBC1_scRNA                     <- RunUMAP(TNBC1_scRNA,dims = 1:30)
TNBC1_scRNA                     <- RunTSNE(TNBC1_scRNA,dims = 1:30,check_duplicates = FALSE)
DimPlot(TNBC1_scRNA,reduction = "umap",label = T)
DimPlot(TNBC1_scRNA,reduction = "tsne",label = T)
save(TNBC1_scRNA,file = "Pro_TNBC/output/data/scRNASeq/Gao.et.al.2021/TNBC1/TNBC1_scRNA.RData")
#run our predicter
TNBC1_exp                         <- as.matrix(TNBC1_scRNA@assays$RNA@data) #normalied data
TNBC1_exp_bct93                   <- as.data.frame(TNBC1_exp)
TNBC1_exp_bct93$SYMBOL            <- rownames(TNBC1_exp_bct93)
TNBC1_exp_bct93                   <- merge(UBS93.data$UBS93.gene.df,TNBC1_exp_bct93,by="SYMBOL")
rown                              <- TNBC1_exp_bct93$ENSEMBL
TNBC1_exp_bct93                   <- TNBC1_exp_bct93[,-c(1:3)] %>% as.matrix()
rownames(TNBC1_exp_bct93)         <- rown
TNBC1_cell_subtype                <- breast.cancer.predictor(expr.of.sample = TNBC1_exp_bct93,
                                                             expr.of.centroid = UBS93.data$UBS93.centroid,
                                                             marker.gene = UBS93.data$UBS93.gene.df$ENSEMBL,
                                                             HER2.amp.signature.genes = UBS93.data$HER2.amp.signature.genes,
                                                             ssgsea.cutoff = UBS93.data$ssgsea.cutoff)
TNBC1_subtype                     <- TNBC1_cell_subtype$subtype
TNBC1_subtype$sample.id           <- rownames(TNBC1_subtype)
TNBC1_uamp                        <- as.data.frame(TNBC1_scRNA@reductions$umap@cell.embeddings)
TNBC1_uamp$sample.id              <- rownames(TNBC1_uamp)
TNBC1_uamp                        <- merge(TNBC1_uamp,TNBC1_subtype,by="sample.id")
ggplot(data = TNBC1_uamp,aes(x=UMAP_1,y=UMAP_2,col=subtype))+
  geom_point()+ggtitle("UAMP plot of TNBC1 by using all genes")
table(TNBC1_subtype$subtype)

####*TNBC2####
library(readr)
TNBC2_matrix                    <- read.delim("Pro_TNBC/data/scRNASeq/Data_Gao2021_Breast/GSE148673_RAW/GSM4476487_combined_UMIcount_CellTypes_TNBC2.txt",sep = "\t")
TNBC2_matrix                    <- as.matrix(TNBC2_matrix)  %>% t() %>% as.data.frame()
table(TNBC2_matrix$copykat.pred)
TNBC2_tumor_matrix              <- TNBC2_matrix[TNBC2_matrix$copykat.pred=="T",]
TNBC2_tumor_matrix              <- TNBC2_tumor_matrix[,-c(1,2)]  %>% as.matrix() %>% t()
TNBC2_scRNA                     <- CreateSeuratObject(counts = TNBC2_tumor_matrix)
#calculate the proportion of ribosomal genes in cells
TNBC2_scRNA[["percent.mt"]]     <- PercentageFeatureSet(TNBC2_scRNA,pattern = "^MT-")
head(TNBC2_scRNA@meta.data,5)
VlnPlot(TNBC2_scRNA,features = c("nCount_RNA","nFeature_RNA","percent.mt"),ncol=3)
TNBC2_scRNA                     <- subset(TNBC2_scRNA,nFeature_RNA >0 & nFeature_RNA<7500)
#normalized
TNBC2_scRNA                     <- NormalizeData(TNBC2_scRNA)
#features choosing
TNBC2_scRNA                     <- FindVariableFeatures(TNBC2_scRNA,selection.method = "vst",nfeatures = 2000)#choose high variability between cells.
#identify the 10 most highly variable genes
top10                           <- head(VariableFeatures(TNBC2_scRNA),10)
#data scaling 
all_gene                        <- rownames(TNBC2_scRNA)
TNBC2_scRNA                     <- ScaleData(TNBC2_scRNA,features = all_gene)#scaled data存放在 pbmc[["RNA"]]@scale.data
#linear dimensionality reduction
TNBC2_scRNA                     <- RunPCA(TNBC2_scRNA,features = VariableFeatures(object=TNBC2_scRNA))
#dimension selection
TNBC2_scRNA                     <- JackStraw(TNBC2_scRNA,num.replicate = 100)
TNBC2_scRNA                     <- ScoreJackStraw(TNBC2_scRNA,dims = 1:20)
JackStrawPlot(TNBC2_scRNA,dims=1:20)
ElbowPlot(TNBC2_scRNA,ndims = 50)#choose dims of 30
#cell clustering
TNBC2_scRNA                     <- FindNeighbors(TNBC2_scRNA,dims = 1:40)
TNBC2_scRNA                     <- FindClusters(TNBC2_scRNA,resolution = 0.1)
table(Idents(TNBC2_scRNA))
#unlinear dimension reduction(UMAP or tSNE)
TNBC2_scRNA                     <- RunUMAP(TNBC2_scRNA,dims = 1:40)
TNBC2_scRNA                     <- RunTSNE(TNBC2_scRNA,dims = 1:40,check_duplicates = FALSE)
DimPlot(TNBC2_scRNA,reduction = "umap",label = T)
DimPlot(TNBC2_scRNA,reduction = "tsne",label = T)
save(TNBC2_scRNA,file = "Pro_TNBC/output/data/scRNASeq/Gao.et.al.2021/TNBC2/TNBC2_scRNA.RData")
#run our predicter
TNBC2_exp                         <- as.matrix(TNBC2_scRNA@assays$RNA@data) #normalied data
TNBC2_exp_bct93                   <- as.data.frame(TNBC2_exp)
TNBC2_exp_bct93$SYMBOL            <- rownames(TNBC2_exp_bct93)
TNBC2_exp_bct93                   <- merge(UBS93.data$UBS93.gene.df,TNBC2_exp_bct93,by="SYMBOL")
rown                              <- TNBC2_exp_bct93$ENSEMBL
TNBC2_exp_bct93                   <- TNBC2_exp_bct93[,-c(1:3)] %>% as.matrix()
rownames(TNBC2_exp_bct93)         <- rown
TNBC2_cell_subtype                <- breast.cancer.predictor(expr.of.sample = TNBC2_exp_bct93,
                                                             expr.of.centroid = UBS93.data$UBS93.centroid,
                                                             marker.gene = UBS93.data$UBS93.gene.df$ENSEMBL,
                                                             HER2.amp.signature.genes = UBS93.data$HER2.amp.signature.genes,
                                                             ssgsea.cutoff = UBS93.data$ssgsea.cutoff)

TNBC2_subtype                     <- TNBC2_cell_subtype$subtype
TNBC2_subtype$sample.id           <- rownames(TNBC2_subtype)
TNBC2_uamp                        <- as.data.frame(TNBC2_scRNA@reductions$umap@cell.embeddings)
TNBC2_uamp$sample.id              <- rownames(TNBC2_uamp)
TNBC2_uamp                        <- merge(TNBC2_uamp,TNBC2_subtype,by="sample.id")
ggplot(data = TNBC2_uamp,aes(x=UMAP_1,y=UMAP_2,col=subtype))+
  geom_point()+ggtitle("UAMP plot of TNBC2 by using all genes")
table(TNBC2_subtype$subtype)

####*TNBC3####
library(readr)
TNBC3_matrix                    <- read.delim("Pro_TNBC/data/scRNASeq/Data_Gao2021_Breast/GSE148673_RAW/GSM4476488_combined_UMIcount_CellTypes_TNBC3.txt",sep = "\t")
TNBC3_matrix                    <- as.matrix(TNBC3_matrix)  %>% t() %>% as.data.frame()
table(TNBC3_matrix$copykat.pred)
TNBC3_tumor_matrix              <- TNBC3_matrix[TNBC3_matrix$copykat.pred=="T",]
TNBC3_tumor_matrix              <- TNBC3_tumor_matrix[,-c(1,2)]  %>% as.matrix() %>% t()
TNBC3_scRNA                     <- CreateSeuratObject(counts = TNBC3_tumor_matrix)
#calculate the proportion of ribosomal genes in cells
TNBC3_scRNA[["percent.mt"]]     <- PercentageFeatureSet(TNBC3_scRNA,pattern = "^MT-")
head(TNBC3_scRNA@meta.data,5)
VlnPlot(TNBC3_scRNA,features = c("nCount_RNA","nFeature_RNA","percent.mt"),ncol=3)
TNBC3_scRNA                     <- subset(TNBC3_scRNA,nFeature_RNA >0 & nFeature_RNA<7500)
#normalized
TNBC3_scRNA                     <- NormalizeData(TNBC3_scRNA)
TNBC3_scRNA[["RNA"]]@data[c("ENSG00000167286", "ENSG00000100721", "ENSG00000156738"), 1:10]#data after normalizing
#features choosing
TNBC3_scRNA                     <- FindVariableFeatures(TNBC3_scRNA,selection.method = "vst",nfeatures = 2000)#choose high variability between cells.
#identify the 10 most highly variable genes
top10                           <- head(VariableFeatures(TNBC3_scRNA),10)
#data scaling 
all_gene                        <- rownames(TNBC3_scRNA)
TNBC3_scRNA                     <- ScaleData(TNBC3_scRNA,features = all_gene)#scaled data存放在 pbmc[["RNA"]]@scale.data
#linear dimensionality reduction
TNBC3_scRNA                     <- RunPCA(TNBC3_scRNA,features = VariableFeatures(object=TNBC3_scRNA))
#dimension selection
TNBC3_scRNA                     <- JackStraw(TNBC3_scRNA,num.replicate = 100)
TNBC3_scRNA                     <- ScoreJackStraw(TNBC3_scRNA,dims = 1:20)
JackStrawPlot(TNBC3_scRNA,dims=1:20)
ElbowPlot(TNBC3_scRNA,ndims = 50)#choose dims of 30
#cell clustering
TNBC3_scRNA                     <- FindNeighbors(TNBC3_scRNA,dims = 1:50)
TNBC3_scRNA                     <- FindClusters(TNBC3_scRNA,resolution = 0.1)
table(Idents(TNBC3_scRNA))
#unlinear dimension reduction(UMAP or tSNE)
TNBC3_scRNA                     <- RunUMAP(TNBC3_scRNA,dims = 1:50)
TNBC3_scRNA                     <- RunTSNE(TNBC3_scRNA,dims = 1:50,check_duplicates = FALSE)
DimPlot(TNBC3_scRNA,reduction = "umap",label = T)
DimPlot(TNBC3_scRNA,reduction = "tsne",label = T)
save(TNBC3_scRNA,file = "Pro_TNBC/output/data/scRNASeq/Gao.et.al.2021/TNBC3/TNBC3_scRNA.RData")

#run our predicter
TNBC3_exp                         <- as.matrix(TNBC3_scRNA@assays$RNA@data) #normalied data
TNBC3_exp_bct93                   <- as.data.frame(TNBC3_exp)
TNBC3_exp_bct93$SYMBOL            <- rownames(TNBC3_exp_bct93)
TNBC3_exp_bct93                   <- merge(UBS93.data$UBS93.gene.df,TNBC3_exp_bct93,by="SYMBOL")
rown                              <- TNBC3_exp_bct93$ENSEMBL
TNBC3_exp_bct93                   <- TNBC3_exp_bct93[,-c(1:3)] %>% as.matrix()
rownames(TNBC3_exp_bct93)         <- rown
TNBC3_cell_subtype                <- breast.cancer.predictor(expr.of.sample = TNBC3_exp_bct93,
                                                             expr.of.centroid = UBS93.data$UBS93.centroid,
                                                             marker.gene = UBS93.data$UBS93.gene.df$ENSEMBL,
                                                             HER2.amp.signature.genes = UBS93.data$HER2.amp.signature.genes,
                                                             ssgsea.cutoff = UBS93.data$ssgsea.cutoff)

TNBC3_subtype                     <- TNBC3_cell_subtype$subtype
TNBC3_subtype$sample.id           <- rownames(TNBC3_subtype)
TNBC3_uamp                        <- as.data.frame(TNBC3_scRNA@reductions$umap@cell.embeddings)
TNBC3_uamp$sample.id              <- rownames(TNBC3_uamp)
TNBC3_uamp                        <- merge(TNBC3_uamp,TNBC3_subtype,by="sample.id")
ggplot(data = TNBC3_uamp,aes(x=UMAP_1,y=UMAP_2,col=subtype))+
  geom_point()+ggtitle("UAMP plot of TNBC3 by using all genes")
table(TNBC3_subtype$subtype)

####*TNBC4####
library(readr)
TNBC4_matrix                    <- read.delim("Pro_TNBC/data/scRNASeq/Data_Gao2021_Breast/GSE148673_RAW/GSM4476489_combined_UMIcount_CellTypes_TNBC4.txt",sep = "\t")
TNBC4_matrix                    <- as.matrix(TNBC4_matrix)  %>% t() %>% as.data.frame()
table(TNBC4_matrix$copykat.pred)
TNBC4_tumor_matrix              <- TNBC4_matrix[TNBC4_matrix$copykat.pred=="T",]
TNBC4_tumor_matrix              <- TNBC4_tumor_matrix[,-c(1,2)]  %>% as.matrix() %>% t()
TNBC4_scRNA                     <- CreateSeuratObject(counts = TNBC4_tumor_matrix)
#calculate the proportion of ribosomal genes in cells
TNBC4_scRNA[["percent.mt"]]     <- PercentageFeatureSet(TNBC4_scRNA,pattern = "^MT-")
head(TNBC4_scRNA@meta.data,5)
VlnPlot(TNBC4_scRNA,features = c("nCount_RNA","nFeature_RNA","percent.mt"),ncol=3)
TNBC4_scRNA                     <- subset(TNBC4_scRNA,nFeature_RNA >0 & nFeature_RNA<7500)
#normalized
TNBC4_scRNA                     <- NormalizeData(TNBC4_scRNA)
#features choosing
TNBC4_scRNA                     <- FindVariableFeatures(TNBC4_scRNA,selection.method = "vst",nfeatures = 2000)#choose high variability between cells.
#identify the 10 most highly variable genes
top10                           <- head(VariableFeatures(TNBC4_scRNA),10)
#data scaling 
all_gene                        <- rownames(TNBC4_scRNA)
TNBC4_scRNA                     <- ScaleData(TNBC4_scRNA,features = all_gene)#scaled data存放在 pbmc[["RNA"]]@scale.data
#linear dimensionality reduction
TNBC4_scRNA                     <- RunPCA(TNBC4_scRNA,features = VariableFeatures(object=TNBC4_scRNA))
#dimension selection
TNBC4_scRNA                     <- JackStraw(TNBC4_scRNA,num.replicate = 100)
TNBC4_scRNA                     <- ScoreJackStraw(TNBC4_scRNA,dims = 1:20)
JackStrawPlot(TNBC4_scRNA,dims=1:20)
ElbowPlot(TNBC4_scRNA,ndims = 50)#choose dims of 30
#cell clustering
TNBC4_scRNA                     <- FindNeighbors(TNBC4_scRNA,dims = 1:40)
TNBC4_scRNA                     <- FindClusters(TNBC4_scRNA,resolution = 0.1)
table(Idents(TNBC4_scRNA))
#unlinear dimension reduction(UMAP or tSNE)
TNBC4_scRNA                     <- RunUMAP(TNBC4_scRNA,dims = 1:40)
TNBC4_scRNA                     <- RunTSNE(TNBC4_scRNA,dims = 1:40,check_duplicates = FALSE)
DimPlot(TNBC4_scRNA,reduction = "umap",label = T)
DimPlot(TNBC4_scRNA,reduction = "tsne",label = T)
save(TNBC4_scRNA,file = "Pro_TNBC/output/data/scRNASeq/Gao.et.al.2021/TNBC4/TNBC4_scRNA.RData")

#run our predicter
TNBC4_exp                         <- as.matrix(TNBC4_scRNA@assays$RNA@data) #normalied data
TNBC4_exp_bct93                   <- as.data.frame(TNBC4_exp)
TNBC4_exp_bct93$SYMBOL            <- rownames(TNBC4_exp_bct93)
TNBC4_exp_bct93                   <- merge(UBS93.data$UBS93.gene.df,TNBC4_exp_bct93,by="SYMBOL")
rown                              <- TNBC4_exp_bct93$ENSEMBL
TNBC4_exp_bct93                   <- TNBC4_exp_bct93[,-c(1:3)] %>% as.matrix()
rownames(TNBC4_exp_bct93)         <- rown
TNBC4_cell_subtype                <- breast.cancer.predictor(expr.of.sample = TNBC4_exp_bct93,
                                                             expr.of.centroid = UBS93.data$UBS93.centroid,
                                                             marker.gene = UBS93.data$UBS93.gene.df$ENSEMBL,
                                                             HER2.amp.signature.genes = UBS93.data$HER2.amp.signature.genes,
                                                             ssgsea.cutoff = UBS93.data$ssgsea.cutoff)

TNBC4_subtype                     <- TNBC4_cell_subtype$subtype
TNBC4_subtype$sample.id           <- rownames(TNBC4_subtype)
TNBC4_uamp                        <- as.data.frame(TNBC4_scRNA@reductions$umap@cell.embeddings)
TNBC4_uamp$sample.id              <- rownames(TNBC4_uamp)
TNBC4_uamp                        <- merge(TNBC4_uamp,TNBC4_subtype,by="sample.id")
ggplot(data = TNBC4_uamp,aes(x=UMAP_1,y=UMAP_2,col=subtype))+
  geom_point()+ggtitle("UAMP plot of TNBC4 by using all genes")
table(TNBC4_subtype$subtype)

####*TNBC5####
library(readr)
TNBC5_matrix                    <- read.delim("Pro_TNBC/data/scRNASeq/Data_Gao2021_Breast/GSE148673_RAW/GSM4476490_combined_UMIcount_CellTypes_TNBC5.txt",sep = "\t")
TNBC5_matrix                    <- as.matrix(TNBC5_matrix)  %>% t() %>% as.data.frame()
table(TNBC5_matrix$copykat.pred)
TNBC5_tumor_matrix              <- TNBC5_matrix[TNBC5_matrix$copykat.pred=="T",]
TNBC5_tumor_matrix              <- TNBC5_tumor_matrix[,-c(1,2)]  %>% as.matrix() %>% t()
TNBC5_scRNA                     <- CreateSeuratObject(counts = TNBC5_tumor_matrix)
#calculate the proportion of ribosomal genes in cells
TNBC5_scRNA[["percent.mt"]]     <- PercentageFeatureSet(TNBC5_scRNA,pattern = "^MT-")
head(TNBC5_scRNA@meta.data,5)
VlnPlot(TNBC5_scRNA,features = c("nCount_RNA","nFeature_RNA","percent.mt"),ncol=3)
TNBC5_scRNA                     <- subset(TNBC5_scRNA,nFeature_RNA >0 & nFeature_RNA<8000)
#normalized
TNBC5_scRNA                     <- NormalizeData(TNBC5_scRNA)
#features choosing
TNBC5_scRNA                     <- FindVariableFeatures(TNBC5_scRNA,selection.method = "vst",nfeatures = 2000)#choose high variability between cells.
#identify the 10 most highly variable genes
top10                           <- head(VariableFeatures(TNBC5_scRNA),10)
#data scaling 
all_gene                        <- rownames(TNBC5_scRNA)
TNBC5_scRNA                     <- ScaleData(TNBC5_scRNA,features = all_gene)#scaled data存放在 pbmc[["RNA"]]@scale.data
#linear dimensionality reduction
TNBC5_scRNA                     <- RunPCA(TNBC5_scRNA,features = VariableFeatures(object=TNBC5_scRNA))
#dimension selection
TNBC5_scRNA                     <- JackStraw(TNBC5_scRNA,num.replicate = 100)
TNBC5_scRNA                     <- ScoreJackStraw(TNBC5_scRNA,dims = 1:20)
JackStrawPlot(TNBC5_scRNA,dims=1:20)
ElbowPlot(TNBC5_scRNA,ndims = 50)#choose dims of 30
#cell clustering
TNBC5_scRNA                     <- FindNeighbors(TNBC5_scRNA,dims = 1:30)
TNBC5_scRNA                     <- FindClusters(TNBC5_scRNA,resolution = 0.1)
table(Idents(TNBC5_scRNA))
#unlinear dimension reduction(UMAP or tSNE)
TNBC5_scRNA                     <- RunUMAP(TNBC5_scRNA,dims = 1:30)
TNBC5_scRNA                     <- RunTSNE(TNBC5_scRNA,dims = 1:30,check_duplicates = FALSE)
DimPlot(TNBC5_scRNA,reduction = "umap",label = T)
DimPlot(TNBC5_scRNA,reduction = "tsne",label = T)
save(TNBC5_scRNA,file = "Pro_TNBC/output/data/scRNASeq/Gao.et.al.2021/TNBC5/TNBC5_scRNA.RData")

#run our predicter
TNBC5_exp                         <- as.matrix(TNBC5_scRNA@assays$RNA@data) #normalied data
TNBC5_exp_bct93                   <- as.data.frame(TNBC5_exp)
TNBC5_exp_bct93$SYMBOL            <- rownames(TNBC5_exp_bct93)
TNBC5_exp_bct93                   <- merge(UBS93.data$UBS93.gene.df,TNBC5_exp_bct93,by="SYMBOL")
rown                              <- TNBC5_exp_bct93$ENSEMBL
TNBC5_exp_bct93                   <- TNBC5_exp_bct93[,-c(1:3)] %>% as.matrix()
rownames(TNBC5_exp_bct93)         <- rown
TNBC5_cell_subtype                <- breast.cancer.predictor(expr.of.sample = TNBC5_exp_bct93,
                                                             expr.of.centroid = UBS93.data$UBS93.centroid,
                                                             marker.gene = UBS93.data$UBS93.gene.df$ENSEMBL,
                                                             HER2.amp.signature.genes = UBS93.data$HER2.amp.signature.genes,
                                                             ssgsea.cutoff = UBS93.data$ssgsea.cutoff)

TNBC5_subtype                     <- TNBC5_cell_subtype$subtype
TNBC5_subtype$sample.id           <- rownames(TNBC5_subtype)
TNBC5_uamp                        <- as.data.frame(TNBC5_scRNA@reductions$umap@cell.embeddings)
TNBC5_uamp$sample.id              <- rownames(TNBC5_uamp)
TNBC5_uamp                        <- merge(TNBC5_uamp,TNBC5_subtype,by="sample.id")
ggplot(data = TNBC5_uamp,aes(x=UMAP_1,y=UMAP_2,col=subtype))+
  geom_point()+ggtitle("UAMP plot of TNBC5 by using all genes")
table(TNBC5_subtype$subtype)


####Qian et al. 2020####
####*sc5rJUQ033####
load("~/Pro_TNBC/output/data/scRNASeq/Qian_et_2020/sc5rJUQ033_scRNA/sc5rJUQ033_scRNA.RData")
sc5rJUQ033_tumor_scRNA                       <- subset(sc5rJUQ033_scRNA,celltype=="Epithelial_cells"&cell=="tumorcell")
#run our predicter
sc5rJUQ033_tumor_exp                         <- as.matrix(sc5rJUQ033_tumor_scRNA@assays$RNA@data) #normalied data
sc5rJUQ033_tumor_exp_bct93                   <- as.data.frame(sc5rJUQ033_tumor_exp)
sc5rJUQ033_tumor_exp_bct93$SYMBOL            <- rownames(sc5rJUQ033_tumor_exp_bct93)
sc5rJUQ033_tumor_exp_bct93                   <- merge(UBS93.data$UBS93.gene.df,sc5rJUQ033_tumor_exp_bct93,by="SYMBOL")
rown                                         <- sc5rJUQ033_tumor_exp_bct93$ENSEMBL
sc5rJUQ033_tumor_exp_bct93                   <- sc5rJUQ033_tumor_exp_bct93[,-c(1:3)] %>% as.matrix()
rownames(sc5rJUQ033_tumor_exp_bct93)         <- rown
sc5rJUQ033_tumor_cell_subtype                <- breast.cancer.predictor(expr.of.sample = sc5rJUQ033_tumor_exp_bct93,
                                                                        expr.of.centroid = UBS93.data$UBS93.centroid,
                                                                        marker.gene = UBS93.data$UBS93.gene.df$ENSEMBL,
                                                                        HER2.amp.signature.genes = UBS93.data$HER2.amp.signature.genes,
                                                                        ssgsea.cutoff = UBS93.data$ssgsea.cutoff)

sc5rJUQ033_tumor_subtype                     <- sc5rJUQ033_tumor_cell_subtype$subtype
sc5rJUQ033_tumor_subtype$sample.id           <- rownames(sc5rJUQ033_tumor_subtype)
sc5rJUQ033_tumor_uamp                        <- as.data.frame(sc5rJUQ033_tumor_scRNA@reductions$umap@cell.embeddings)
sc5rJUQ033_tumor_uamp$sample.id              <- rownames(sc5rJUQ033_tumor_uamp)
sc5rJUQ033_tumor_uamp                        <- merge(sc5rJUQ033_tumor_uamp,sc5rJUQ033_tumor_subtype,by="sample.id")
ggplot(data = sc5rJUQ033_tumor_uamp,aes(x=UMAP_1,y=UMAP_2,col=subtype))+
  geom_point()+ggtitle("UAMP plot of sc5rJUQ033_tumor by using all genes")
table(sc5rJUQ033_tumor_subtype$subtype)
sc5rJUQ033_cl_id                             <- sc5rJUQ033_tumor_subtype[sc5rJUQ033_tumor_subtype$subtype=="Claudin_low",]$sample.id
sc5rJUQ033_cl_cor                            <- sc5rJUQ033_tumor_cell_subtype$cor.matrix[ ,2,drop=F] %>% as.data.frame()
sc5rJUQ033_cl_cor                            <- sc5rJUQ033_cl_cor[rownames(sc5rJUQ033_cl_cor) %in%sc5rJUQ033_cl_id,,drop=F]
sc5rJUQ033_cl_cor$patient.id                 <- rep("sc5rJUQ033",length(sc5rJUQ033_cl_id))

####*sc5rJUQ039####
load("~/Pro_TNBC/output/data/scRNASeq/Qian_et_2020/sc5rJUQ039_scRNA/sc5rJUQ039_scRNA.RData")
sc5rJUQ039_tumor_scRNA                       <- subset(sc5rJUQ039_scRNA,celltype=="Epithelial_cells"&cell=="tumorcell")
#run our predicter
sc5rJUQ039_tumor_exp                         <- as.matrix(sc5rJUQ039_tumor_scRNA@assays$RNA@data) #normalied data
sc5rJUQ039_tumor_exp_bct93                   <- as.data.frame(sc5rJUQ039_tumor_exp)
sc5rJUQ039_tumor_exp_bct93$SYMBOL            <- rownames(sc5rJUQ039_tumor_exp_bct93)
sc5rJUQ039_tumor_exp_bct93                   <- merge(UBS93.data$UBS93.gene.df,sc5rJUQ039_tumor_exp_bct93,by="SYMBOL")
rown                                         <- sc5rJUQ039_tumor_exp_bct93$ENSEMBL
sc5rJUQ039_tumor_exp_bct93                   <- sc5rJUQ039_tumor_exp_bct93[,-c(1:3)] %>% as.matrix()
rownames(sc5rJUQ039_tumor_exp_bct93)         <- rown
sc5rJUQ039_tumor_cell_subtype                <- breast.cancer.predictor(expr.of.sample = sc5rJUQ039_tumor_exp_bct93,
                                                                        expr.of.centroid = UBS93.data$UBS93.centroid,
                                                                        marker.gene = UBS93.data$UBS93.gene.df$ENSEMBL,
                                                                        HER2.amp.signature.genes = UBS93.data$HER2.amp.signature.genes,
                                                                        ssgsea.cutoff = UBS93.data$ssgsea.cutoff)

sc5rJUQ039_tumor_subtype                     <- sc5rJUQ039_tumor_cell_subtype$subtype
sc5rJUQ039_tumor_subtype$sample.id           <- rownames(sc5rJUQ039_tumor_subtype)
sc5rJUQ039_tumor_uamp                        <- as.data.frame(sc5rJUQ039_tumor_scRNA@reductions$umap@cell.embeddings)
sc5rJUQ039_tumor_uamp$sample.id              <- rownames(sc5rJUQ039_tumor_uamp)
sc5rJUQ039_tumor_uamp                        <- merge(sc5rJUQ039_tumor_uamp,sc5rJUQ039_tumor_subtype,by="sample.id")
ggplot(data = sc5rJUQ039_tumor_uamp,aes(x=UMAP_1,y=UMAP_2,col=subtype))+
  geom_point()+ggtitle("UAMP plot of sc5rJUQ039_tumor by using all genes")
table(sc5rJUQ039_tumor_subtype$subtype)
sc5rJUQ039_cl_id                             <- sc5rJUQ039_tumor_subtype[sc5rJUQ039_tumor_subtype$subtype=="Claudin_low",]$sample.id
sc5rJUQ039_cl_cor                            <- sc5rJUQ039_tumor_cell_subtype$cor.matrix[ ,2,drop=F] %>% as.data.frame()
sc5rJUQ039_cl_cor                            <- sc5rJUQ039_cl_cor[rownames(sc5rJUQ039_cl_cor) %in%sc5rJUQ039_cl_id,,drop=F]
sc5rJUQ039_cl_cor$patient.id                 <- rep("sc5rJUQ039",length(sc5rJUQ039_cl_id))

####*sc5rJUQ042####
load("~/Pro_TNBC/output/data/scRNASeq/Qian_et_2020/sc5rJUQ042_scRNA/sc5rJUQ042_scRNA.RData")
sc5rJUQ042_tumor_scRNA                       <- subset(sc5rJUQ042_scRNA,celltype=="Epithelial_cells"&cell=="tumorcell")
#run our predicter
sc5rJUQ042_tumor_exp                         <- as.matrix(sc5rJUQ042_tumor_scRNA@assays$RNA@data) #normalied data
sc5rJUQ042_tumor_exp_bct93                   <- as.data.frame(sc5rJUQ042_tumor_exp)
sc5rJUQ042_tumor_exp_bct93$SYMBOL            <- rownames(sc5rJUQ042_tumor_exp_bct93)
sc5rJUQ042_tumor_exp_bct93                   <- merge(UBS93.data$UBS93.gene.df,sc5rJUQ042_tumor_exp_bct93,by="SYMBOL")
rown                                         <- sc5rJUQ042_tumor_exp_bct93$ENSEMBL
sc5rJUQ042_tumor_exp_bct93                   <- sc5rJUQ042_tumor_exp_bct93[,-c(1:3)] %>% as.matrix()
rownames(sc5rJUQ042_tumor_exp_bct93)         <- rown
sc5rJUQ042_tumor_cell_subtype                <- breast.cancer.predictor(expr.of.sample = sc5rJUQ042_tumor_exp_bct93,
                                                                        expr.of.centroid = UBS93.data$UBS93.centroid,
                                                                        marker.gene = UBS93.data$UBS93.gene.df$ENSEMBL,
                                                                        HER2.amp.signature.genes = UBS93.data$HER2.amp.signature.genes,
                                                                        ssgsea.cutoff = UBS93.data$ssgsea.cutoff)

sc5rJUQ042_tumor_subtype                     <- sc5rJUQ042_tumor_cell_subtype$subtype
sc5rJUQ042_tumor_subtype$sample.id           <- rownames(sc5rJUQ042_tumor_subtype)
sc5rJUQ042_tumor_uamp                        <- as.data.frame(sc5rJUQ042_tumor_scRNA@reductions$umap@cell.embeddings)
sc5rJUQ042_tumor_uamp$sample.id              <- rownames(sc5rJUQ042_tumor_uamp)
sc5rJUQ042_tumor_uamp                        <- merge(sc5rJUQ042_tumor_uamp,sc5rJUQ042_tumor_subtype,by="sample.id")
ggplot(data = sc5rJUQ042_tumor_uamp,aes(x=UMAP_1,y=UMAP_2,col=subtype))+
  geom_point()+ggtitle("UAMP plot of sc5rJUQ042_tumor by using all genes")
table(sc5rJUQ042_tumor_subtype$subtype)
sc5rJUQ042_cl_id                             <- sc5rJUQ042_tumor_subtype[sc5rJUQ042_tumor_subtype$subtype=="Claudin_low",]$sample.id
sc5rJUQ042_cl_cor                            <- sc5rJUQ042_tumor_cell_subtype$cor.matrix[ ,2,drop=F] %>% as.data.frame()
sc5rJUQ042_cl_cor                            <- sc5rJUQ042_cl_cor[rownames(sc5rJUQ042_cl_cor) %in%sc5rJUQ042_cl_id,,drop=F]
sc5rJUQ042_cl_cor$patient.id                 <- rep("sc5rJUQ042",length(sc5rJUQ042_cl_id))

####*sc5rJUQ045####
load("~/Pro_TNBC/output/data/scRNASeq/Qian_et_2020/sc5rJUQ045_scRNA/sc5rJUQ045_scRNA.RData")
sc5rJUQ045_tumor_scRNA                       <- subset(sc5rJUQ045_scRNA,celltype=="epithelial_cells"&cell=="tumorcell")
#run our predicter
sc5rJUQ045_tumor_exp                         <- as.matrix(sc5rJUQ045_tumor_scRNA@assays$RNA@data) #normalied data
sc5rJUQ045_tumor_exp_bct93                   <- as.data.frame(sc5rJUQ045_tumor_exp)
sc5rJUQ045_tumor_exp_bct93$SYMBOL            <- rownames(sc5rJUQ045_tumor_exp_bct93)
sc5rJUQ045_tumor_exp_bct93                   <- merge(UBS93.data$UBS93.gene.df,sc5rJUQ045_tumor_exp_bct93,by="SYMBOL")
rown                                         <- sc5rJUQ045_tumor_exp_bct93$ENSEMBL
sc5rJUQ045_tumor_exp_bct93                   <- sc5rJUQ045_tumor_exp_bct93[,-c(1:3)] %>% as.matrix()
rownames(sc5rJUQ045_tumor_exp_bct93)         <- rown
sc5rJUQ045_tumor_cell_subtype                <- breast.cancer.predictor(expr.of.sample = sc5rJUQ045_tumor_exp_bct93,
                                                                        expr.of.centroid = UBS93.data$UBS93.centroid,
                                                                        marker.gene = UBS93.data$UBS93.gene.df$ENSEMBL,
                                                                        HER2.amp.signature.genes = UBS93.data$HER2.amp.signature.genes,
                                                                        ssgsea.cutoff = UBS93.data$ssgsea.cutoff)

sc5rJUQ045_tumor_subtype                     <- sc5rJUQ045_tumor_cell_subtype$subtype
sc5rJUQ045_tumor_subtype$sample.id           <- rownames(sc5rJUQ045_tumor_subtype)
sc5rJUQ045_tumor_uamp                        <- as.data.frame(sc5rJUQ045_tumor_scRNA@reductions$umap@cell.embeddings)
sc5rJUQ045_tumor_uamp$sample.id              <- rownames(sc5rJUQ045_tumor_uamp)
sc5rJUQ045_tumor_uamp                        <- merge(sc5rJUQ045_tumor_uamp,sc5rJUQ045_tumor_subtype,by="sample.id")
ggplot(data = sc5rJUQ045_tumor_uamp,aes(x=UMAP_1,y=UMAP_2,col=subtype))+
  geom_point()+ggtitle("UAMP plot of sc5rJUQ045_tumor by using all genes")
table(sc5rJUQ045_tumor_subtype$subtype)
sc5rJUQ045_cl_id                             <- sc5rJUQ045_tumor_subtype[sc5rJUQ045_tumor_subtype$subtype=="Claudin_low",]$sample.id
sc5rJUQ045_cl_cor                            <- sc5rJUQ045_tumor_cell_subtype$cor.matrix[ ,2,drop=F] %>% as.data.frame()
sc5rJUQ045_cl_cor                            <- sc5rJUQ045_cl_cor[rownames(sc5rJUQ045_cl_cor) %in%sc5rJUQ045_cl_id,,drop=F]
sc5rJUQ045_cl_cor$patient.id                 <- rep("sc5rJUQ045",length(sc5rJUQ045_cl_id))

####*sc5rJUQ053####
load("~/Pro_TNBC/output/data/scRNASeq/Qian_et_2020/sc5rJUQ053_scRNA/sc5rJUQ053_scRNA.RData")
sc5rJUQ053_tumor_scRNA                       <- subset(sc5rJUQ053_scRNA,celltype=="Epithelial_cells"&cell=="tumorcell")
#run our predicter
sc5rJUQ053_tumor_exp                         <- as.matrix(sc5rJUQ053_tumor_scRNA@assays$RNA@data) #normalied data
sc5rJUQ053_tumor_exp_bct93                   <- as.data.frame(sc5rJUQ053_tumor_exp)
sc5rJUQ053_tumor_exp_bct93$SYMBOL            <- rownames(sc5rJUQ053_tumor_exp_bct93)
sc5rJUQ053_tumor_exp_bct93                   <- merge(UBS93.data$UBS93.gene.df,sc5rJUQ053_tumor_exp_bct93,by="SYMBOL")
rown                                         <- sc5rJUQ053_tumor_exp_bct93$ENSEMBL
sc5rJUQ053_tumor_exp_bct93                   <- sc5rJUQ053_tumor_exp_bct93[,-c(1:3)] %>% as.matrix()
rownames(sc5rJUQ053_tumor_exp_bct93)         <- rown
sc5rJUQ053_tumor_cell_subtype                <- breast.cancer.predictor(expr.of.sample = sc5rJUQ053_tumor_exp_bct93,
                                                                        expr.of.centroid = UBS93.data$UBS93.centroid,
                                                                        marker.gene = UBS93.data$UBS93.gene.df$ENSEMBL,
                                                                        HER2.amp.signature.genes = UBS93.data$HER2.amp.signature.genes,
                                                                        ssgsea.cutoff = UBS93.data$ssgsea.cutoff)

sc5rJUQ053_tumor_subtype                     <- sc5rJUQ053_tumor_cell_subtype$subtype
sc5rJUQ053_tumor_subtype$sample.id           <- rownames(sc5rJUQ053_tumor_subtype)
sc5rJUQ053_tumor_uamp                        <- as.data.frame(sc5rJUQ053_tumor_scRNA@reductions$umap@cell.embeddings)
sc5rJUQ053_tumor_uamp$sample.id              <- rownames(sc5rJUQ053_tumor_uamp)
sc5rJUQ053_tumor_uamp                        <- merge(sc5rJUQ053_tumor_uamp,sc5rJUQ053_tumor_subtype,by="sample.id")
ggplot(data = sc5rJUQ053_tumor_uamp,aes(x=UMAP_1,y=UMAP_2,col=subtype))+
  geom_point()+ggtitle("UAMP plot of sc5rJUQ053_tumor by using all genes")
table(sc5rJUQ053_tumor_subtype$subtype)
sc5rJUQ053_cl_id                             <- sc5rJUQ053_tumor_subtype[sc5rJUQ053_tumor_subtype$subtype=="Claudin_low",]$sample.id
sc5rJUQ053_cl_cor                            <- sc5rJUQ053_tumor_cell_subtype$cor.matrix[ ,2,drop=F] %>% as.data.frame()
sc5rJUQ053_cl_cor                            <- sc5rJUQ053_cl_cor[rownames(sc5rJUQ053_cl_cor) %in%sc5rJUQ053_cl_id,,drop=F]
sc5rJUQ053_cl_cor$patient.id                 <- rep("sc5rJUQ053",length(sc5rJUQ053_cl_id))
Qian_cl_cor                                  <- bind_rows(sc5rJUQ033_cl_cor, sc5rJUQ039_cl_cor, sc5rJUQ042_cl_cor,sc5rJUQ045_cl_cor,sc5rJUQ053_cl_cor)
table(Qian_cl_cor$patient.id)


####GSE118389####
####*PT039####
load("~/Pro_TNBC/output/data/scRNASeq/Karaayvas_et_2018/PT039/PT039_scRNA.RData")
##tumor cell
PT039_tumor_scRNA                       <- subset(PT039_scRNA,celltype=="Epithelial_cells"&cell=="tumorcell")
PT039_tumor_exp                         <- as.matrix(PT039_tumor_scRNA@assays$RNA@data) #normalied data
PT039_tumor_exp_bct93                   <- as.data.frame(PT039_tumor_exp)
PT039_tumor_exp_bct93$SYMBOL            <- rownames(PT039_tumor_exp_bct93)
PT039_tumor_exp_bct93                   <- merge(UBS93.data$UBS93.gene.df,PT039_tumor_exp_bct93,by="SYMBOL")
rown                                    <- PT039_tumor_exp_bct93$ENSEMBL
PT039_tumor_exp_bct93                   <- PT039_tumor_exp_bct93[,-c(1:3)] %>% as.matrix()
rownames(PT039_tumor_exp_bct93)         <- rown
PT039_tumor_cell_subtype                <- breast.cancer.predictor(expr.of.sample = PT039_tumor_exp_bct93,
                                                                   expr.of.centroid = UBS93.data$UBS93.centroid,
                                                                   marker.gene = UBS93.data$UBS93.gene.df$ENSEMBL,
                                                                   HER2.amp.signature.genes = UBS93.data$HER2.amp.signature.genes,
                                                                   ssgsea.cutoff = UBS93.data$ssgsea.cutoff)

PT039_tumor_subtype                     <- PT039_tumor_cell_subtype$subtype
PT039_tumor_subtype$sample.id           <- rownames(PT039_tumor_subtype)
PT039_tumor_uamp                        <- as.data.frame(PT039_tumor_scRNA@reductions$umap@cell.embeddings)
PT039_tumor_uamp$sample.id              <- rownames(PT039_tumor_uamp)
PT039_tumor_uamp                        <- merge(PT039_tumor_uamp,PT039_tumor_subtype,by="sample.id")
ggplot(data = PT039_tumor_uamp,aes(x=UMAP_1,y=UMAP_2,col=subtype))+
  geom_point()+ggtitle("UAMP plot of PT039_tumor by using all genes")
table(PT039_tumor_subtype$subtype)
PT039_cl_id                             <- PT039_tumor_subtype[PT039_tumor_subtype$subtype=="Claudin_low",]$sample.id
PT039_cl_cor                            <- PT039_tumor_cell_subtype$cor.matrix[ ,2,drop=F] %>% as.data.frame()
PT039_cl_cor                            <- PT039_cl_cor[rownames(PT039_cl_cor) %in%PT039_cl_id,,drop=F]
PT039_cl_cor$patient.id                 <- rep("PT039",length(PT039_cl_id))

####*PT084####
load("~/Pro_TNBC/output/data/scRNASeq/Karaayvas_et_2018/PT084/PT084_scRNA.RData")
##tumor cell
PT084_tumor_scRNA                       <- subset(PT084_scRNA,celltype=="Epithelial_cells"&cell=="tumorcell")
PT084_tumor_exp                         <- as.matrix(PT084_tumor_scRNA@assays$RNA@data) #normalied data
PT084_tumor_exp_bct93                   <- as.data.frame(PT084_tumor_exp)
PT084_tumor_exp_bct93$SYMBOL            <- rownames(PT084_tumor_exp_bct93)
PT084_tumor_exp_bct93                   <- merge(UBS93.data$UBS93.gene.df,PT084_tumor_exp_bct93,by="SYMBOL")
rown                                    <- PT084_tumor_exp_bct93$ENSEMBL
PT084_tumor_exp_bct93                   <- PT084_tumor_exp_bct93[,-c(1:3)] %>% as.matrix()
rownames(PT084_tumor_exp_bct93)         <- rown
PT084_tumor_cell_subtype                <- breast.cancer.predictor(expr.of.sample = PT084_tumor_exp_bct93,
                                                                   expr.of.centroid = UBS93.data$UBS93.centroid,
                                                                   marker.gene = UBS93.data$UBS93.gene.df$ENSEMBL,
                                                                   HER2.amp.signature.genes = UBS93.data$HER2.amp.signature.genes,
                                                                   ssgsea.cutoff = UBS93.data$ssgsea.cutoff)

PT084_tumor_subtype                     <- PT084_tumor_cell_subtype$subtype
PT084_tumor_subtype$sample.id           <- rownames(PT084_tumor_subtype)
PT084_tumor_uamp                        <- as.data.frame(PT084_tumor_scRNA@reductions$umap@cell.embeddings)
PT084_tumor_uamp$sample.id              <- rownames(PT084_tumor_uamp)
PT084_tumor_uamp                        <- merge(PT084_tumor_uamp,PT084_tumor_subtype,by="sample.id")
ggplot(data = PT084_tumor_uamp,aes(x=UMAP_1,y=UMAP_2,col=subtype))+
  geom_point()+ggtitle("UAMP plot of PT084_tumor by using all genes")
table(PT084_tumor_subtype$subtype)
PT084_cl_id                             <- PT084_tumor_subtype[PT084_tumor_subtype$subtype=="Claudin_low",]$sample.id
PT084_cl_cor                            <- PT084_tumor_cell_subtype$cor.matrix[ ,2,drop=F] %>% as.data.frame()
PT084_cl_cor                            <- PT084_cl_cor[rownames(PT084_cl_cor) %in%PT084_cl_id,,drop=F]
PT084_cl_cor$patient.id                 <- rep("PT084",length(PT084_cl_id))

####*PT089####
load("~/Pro_TNBC/output/data/scRNASeq/Karaayvas_et_2018/PT089/PT089_scRNA.RData")
##tumor cell
PT089_tumor_scRNA                       <- subset(PT089_scRNA,celltype=="Epithelial_cells"&cell=="tumorcell")
PT089_tumor_exp                         <- as.matrix(PT089_tumor_scRNA@assays$RNA@data) #normalied data
PT089_tumor_exp_bct93                   <- as.data.frame(PT089_tumor_exp)
PT089_tumor_exp_bct93$SYMBOL            <- rownames(PT089_tumor_exp_bct93)
PT089_tumor_exp_bct93                   <- merge(UBS93.data$UBS93.gene.df,PT089_tumor_exp_bct93,by="SYMBOL")
rown                                    <- PT089_tumor_exp_bct93$ENSEMBL
PT089_tumor_exp_bct93                   <- PT089_tumor_exp_bct93[,-c(1:3)] %>% as.matrix()
rownames(PT089_tumor_exp_bct93)         <- rown
PT089_tumor_cell_subtype                <- breast.cancer.predictor(expr.of.sample = PT089_tumor_exp_bct93,
                                                                   expr.of.centroid = UBS93.data$UBS93.centroid,
                                                                   marker.gene = UBS93.data$UBS93.gene.df$ENSEMBL,
                                                                   HER2.amp.signature.genes = UBS93.data$HER2.amp.signature.genes,
                                                                   ssgsea.cutoff = UBS93.data$ssgsea.cutoff)

PT089_tumor_subtype                     <- PT089_tumor_cell_subtype$subtype
PT089_tumor_subtype$sample.id           <- rownames(PT089_tumor_subtype)
PT089_tumor_uamp                        <- as.data.frame(PT089_tumor_scRNA@reductions$umap@cell.embeddings)
PT089_tumor_uamp$sample.id              <- rownames(PT089_tumor_uamp)
PT089_tumor_uamp                        <- merge(PT089_tumor_uamp,PT089_tumor_subtype,by="sample.id")
ggplot(data = PT089_tumor_uamp,aes(x=UMAP_1,y=UMAP_2,col=subtype))+
  geom_point()+ggtitle("UAMP plot of PT089_tumor by using all genes")
table(PT089_tumor_subtype$subtype)
PT089_cl_id                             <- PT089_tumor_subtype[PT089_tumor_subtype$subtype=="Claudin_low",]$sample.id
PT089_cl_cor                            <- PT089_tumor_cell_subtype$cor.matrix[ ,2,drop=F] %>% as.data.frame()
PT089_cl_cor                            <- PT089_cl_cor[rownames(PT089_cl_cor) %in%PT089_cl_id,,drop=F]
PT089_cl_cor$patient.id                 <- rep("PT089",length(PT089_cl_id))

####*PT126####
load("~/Pro_TNBC/output/data/scRNASeq/Karaayvas_et_2018/PT126/PT126_scRNA.RData")
##tumor cell
PT126_tumor_scRNA                       <- subset(PT126_scRNA,celltype=="Epithelial_cells"&cell=="tumorcell")
PT126_tumor_exp                         <- as.matrix(PT126_tumor_scRNA@assays$RNA@data) #normalied data
PT126_tumor_exp_bct93                   <- as.data.frame(PT126_tumor_exp)
PT126_tumor_exp_bct93$SYMBOL            <- rownames(PT126_tumor_exp_bct93)
PT126_tumor_exp_bct93                   <- merge(UBS93.data$UBS93.gene.df,PT126_tumor_exp_bct93,by="SYMBOL")
rown                                    <- PT126_tumor_exp_bct93$ENSEMBL
PT126_tumor_exp_bct93                   <- PT126_tumor_exp_bct93[,-c(1:3)] %>% as.matrix()
rownames(PT126_tumor_exp_bct93)         <- rown
PT126_tumor_cell_subtype                <- breast.cancer.predictor(expr.of.sample = PT126_tumor_exp_bct93,
                                                                   expr.of.centroid = UBS93.data$UBS93.centroid,
                                                                   marker.gene = UBS93.data$UBS93.gene.df$ENSEMBL,
                                                                   HER2.amp.signature.genes = UBS93.data$HER2.amp.signature.genes,
                                                                   ssgsea.cutoff = UBS93.data$ssgsea.cutoff)

PT126_tumor_subtype                     <- PT126_tumor_cell_subtype$subtype
PT126_tumor_subtype$sample.id           <- rownames(PT126_tumor_subtype)
PT126_tumor_uamp                        <- as.data.frame(PT126_tumor_scRNA@reductions$umap@cell.embeddings)
PT126_tumor_uamp$sample.id              <- rownames(PT126_tumor_uamp)
PT126_tumor_uamp                        <- merge(PT126_tumor_uamp,PT126_tumor_subtype,by="sample.id")
ggplot(data = PT126_tumor_uamp,aes(x=UMAP_1,y=UMAP_2,col=subtype))+
  geom_point()+ggtitle("UAMP plot of PT126_tumor by using all genes")
table(PT126_tumor_subtype$subtype)
PT126_cl_id                             <- PT126_tumor_subtype[PT126_tumor_subtype$subtype=="Claudin_low",]$sample.id
PT126_cl_cor                            <- PT126_tumor_cell_subtype$cor.matrix[ ,2,drop=F] %>% as.data.frame()
PT126_cl_cor                            <- PT126_cl_cor[rownames(PT126_cl_cor) %in%PT126_cl_id,,drop=F]
PT126_cl_cor$patient.id                 <- rep("PT126",length(PT126_cl_id))


####GSE138536####
GSE138536_HBC_metadata              <- read.delim("Pro_TNBC/data/scRNASeq/GSE138536/GSE138536_HBC_metadata.txt",sep = "\t")
GSE138536_TNBC_metadata             <- subset(GSE138536_HBC_metadata,GSE138536_HBC_metadata$Clinical_subtype=="Basal")
GSE138536_TNBC_tumor                <- subset(GSE138536_TNBC_metadata,GSE138536_TNBC_metadata$tumor_vs_normal=="T")
GSE138536_HBC_TPM                   <- read.delim("Pro_TNBC/data/scRNASeq/GSE138536/GSE138536_HBC_TranscriptMatrixSalmon_TPM.txt",sep = "\t")
GSE138536_HBC_TPM                   <- as.matrix(GSE138536_HBC_TPM)
save(GSE138536_HBC_TPM,file = "Pro_TNBC/data/scRNASeq/GSE138536/GSE138536_HBC_TPM.RData")

####*SU4####
SU4_tumor_metadata                  <- subset(GSE138536_TNBC_tumor,GSE138536_TNBC_tumor$Patient=="SU4")
SU4_TPM                             <- GSE138536_HBC_TPM[,colnames(GSE138536_HBC_TPM) %in% SU4_tumor_metadata$UniqueID]
library(Seurat)
SU4_scRNA                           <- CreateSeuratObject(counts = SU4_TPM)
SU4_tumor_scRNA                     <- SU4_scRNA
#calculate the proportion of ribosomal genes in cells
SU4_tumor_scRNA[["percent.mt"]]     <- PercentageFeatureSet(SU4_tumor_scRNA,pattern = "^MT-")
head(SU4_tumor_scRNA@meta.data,5)
VlnPlot(SU4_tumor_scRNA,features = c("nCount_RNA","nFeature_RNA","percent.mt"),ncol=3)
SU4_tumor_scRNA                       <- subset(SU4_tumor_scRNA,nFeature_RNA >0 & nFeature_RNA<10000 )
#normalized
SU4_tumor_scRNA                       <- NormalizeData(SU4_tumor_scRNA)
#run our predicter
SU4_tumor_exp                         <- as.matrix(SU4_tumor_scRNA@assays$RNA@data) #normalied data
SU4_tumor_exp_bct93                   <- as.data.frame(SU4_tumor_exp)
SU4_tumor_exp_bct93$SYMBOL            <- rownames(SU4_tumor_exp_bct93)
SU4_tumor_exp_bct93                   <- merge(UBS93.data$UBS93.gene.df,SU4_tumor_exp_bct93,by="SYMBOL")
rown                                  <- SU4_tumor_exp_bct93$ENSEMBL
SU4_tumor_exp_bct93                   <- SU4_tumor_exp_bct93[,-c(1:3)] %>% as.matrix()
rownames(SU4_tumor_exp_bct93)         <- rown
SU4_tumor_cell_subtype                <- breast.cancer.predictor(expr.of.sample = SU4_tumor_exp_bct93,
                                                                 expr.of.centroid = UBS93.data$UBS93.centroid,
                                                                 marker.gene = UBS93.data$UBS93.gene.df$ENSEMBL,
                                                                 HER2.amp.signature.genes = UBS93.data$HER2.amp.signature.genes,
                                                                 ssgsea.cutoff = UBS93.data$ssgsea.cutoff)

SU4_tumor_subtype                     <- SU4_tumor_cell_subtype$subtype
SU4_tumor_subtype$sample.id           <- rownames(SU4_tumor_subtype)
SU4_tumor_uamp                        <- as.data.frame(SU4_tumor_scRNA@reductions$umap@cell.embeddings)
SU4_tumor_uamp$sample.id <- rownames(SU4_tumor_uamp)
SU4_tumor_uamp                        <- merge(SU4_tumor_uamp,SU4_tumor_subtype,by="sample.id")
ggplot(data = SU4_tumor_uamp,aes(x=UMAP_1,y=UMAP_2,col=subtype))+
  geom_point()+ggtitle("UAMP plot of SU4_tumor by using all genes")
table(SU4_tumor_subtype$subtype)

SU4_cl_id                             <- SU4_tumor_subtype[SU4_tumor_subtype$subtype=="Claudin_low",]$sample.id
SU4_cl_cor                            <- SU4_tumor_cell_subtype$cor.matrix[ ,2,drop=F] %>% as.data.frame()
SU4_cl_cor                            <- SU4_cl_cor[rownames(SU4_cl_cor) %in%SU4_cl_id,,drop=F]
SU4_cl_cor$patient.id                 <- rep("SU4",length(SU4_cl_id))

####*SU58####
SU58_tumor_metadata                  <- subset(GSE138536_TNBC_tumor,GSE138536_TNBC_tumor$Patient=="SU58")
SU58_TPM                             <- GSE138536_HBC_TPM[,colnames(GSE138536_HBC_TPM) %in% SU58_tumor_metadata$UniqueID]
library(Seurat)
SU58_scRNA                           <- CreateSeuratObject(counts = SU58_TPM)
SU58_tumor_scRNA                     <- SU58_scRNA
#calculate the proportion of ribosomal genes in cells
SU58_tumor_scRNA[["percent.mt"]]     <- PercentageFeatureSet(SU58_tumor_scRNA,pattern = "^MT-")
head(SU58_tumor_scRNA@meta.data,5)
VlnPlot(SU58_tumor_scRNA,features = c("nCount_RNA","nFeature_RNA","percent.mt"),ncol=3)
SU58_tumor_scRNA                       <- subset(SU58_tumor_scRNA,nFeature_RNA >0 & nFeature_RNA<10000 )
#normalized
SU58_tumor_scRNA                       <- NormalizeData(SU58_tumor_scRNA)
#run our predicter
SU58_tumor_exp                         <- as.matrix(SU58_tumor_scRNA@assays$RNA@data) #normalied data
SU58_tumor_exp_bct93                   <- as.data.frame(SU58_tumor_exp)
SU58_tumor_exp_bct93$SYMBOL            <- rownames(SU58_tumor_exp_bct93)
SU58_tumor_exp_bct93                   <- merge(UBS93.data$UBS93.gene.df,SU58_tumor_exp_bct93,by="SYMBOL")
rown                                   <- SU58_tumor_exp_bct93$ENSEMBL
SU58_tumor_exp_bct93                   <- SU58_tumor_exp_bct93[,-c(1:3)] %>% as.matrix()
rownames(SU58_tumor_exp_bct93)         <- rown
SU58_tumor_cell_subtype                <- breast.cancer.predictor(expr.of.sample = SU58_tumor_exp_bct93,
                                                                  expr.of.centroid = UBS93.data$UBS93.centroid,
                                                                  marker.gene = UBS93.data$UBS93.gene.df$ENSEMBL,
                                                                  HER2.amp.signature.genes = UBS93.data$HER2.amp.signature.genes,
                                                                  ssgsea.cutoff = UBS93.data$ssgsea.cutoff)

SU58_tumor_subtype                     <- SU58_tumor_cell_subtype$subtype
SU58_tumor_subtype$sample.id           <- rownames(SU58_tumor_subtype)
SU58_tumor_uamp                        <- as.data.frame(SU58_tumor_scRNA@reductions$umap@cell.embeddings)
SU58_tumor_uamp$sample.id              <- rownames(SU58_tumor_uamp)
SU58_tumor_uamp                        <- merge(SU58_tumor_uamp,SU58_tumor_subtype,by="sample.id")
ggplot(data = SU58_tumor_uamp,aes(x=UMAP_1,y=UMAP_2,col=subtype))+
  geom_point()+ggtitle("UAMP plot of SU58_tumor by using all genes")
table(SU58_tumor_subtype$subtype)
SU58_cl_id                              <- SU58_tumor_subtype[SU58_tumor_subtype$subtype=="Claudin_low",]$sample.id
SU58_cl_cor                             <- SU58_tumor_cell_subtype$cor.matrix[ ,2,drop=F] %>% as.data.frame()
SU58_cl_cor                             <- SU58_cl_cor[rownames(SU58_cl_cor) %in%SU58_cl_id,,drop=F]
SU58_cl_cor$patient.id                  <- rep("SU58",length(SU58_cl_id))

GSE118389_cl_cor                        <- bind_rows(PT039_cl_cor,PT089_cl_cor,PT084_cl_cor,PT126_cl_cor,SU4_cl_cor,SU58_cl_cor)

####Bassez####
####*cohor1####
library(Matrix)
library(Seurat)
library(dplyr)
library(ggplot2)
library(GSVA)
library(readr)
counts_cells_cohort1                   <- readRDS("~/Pro_TNBC/data/scRNASeq/Data_Bassez2021_Breast/1863-counts_cells_cohort1.rds")
counts_cells_cohort1                   <- as.matrix(counts_cells_cohort1)
BIOKEY_metaData_cohort1_web            <- read_csv("Pro_TNBC/data/scRNASeq/Data_Bassez2021_Breast/1872-BIOKEY_metaData_cohort1_web.csv")
BIOKEY_metaData_cohort1_TNBC           <- subset(BIOKEY_metaData_cohort1_web,BC_type=="TNBC")
table(BIOKEY_metaData_cohort1_TNBC$patient_id)

patient_id                   <- names(table(BIOKEY_metaData_cohort1_TNBC$patient_id))
Bassez_subtype_result        <- NULL
for (i in 1:length(patient_id)) {
  A                   <- patient_id[i] %>% as.character()
  tumor_cells         <- subset(BIOKEY_metaData_cohort1_TNBC,patient_id== A &cellType=="Cancer_cell"& timepoint=="Pre")
  tumor_exp           <- counts_cells_cohort1[,colnames(counts_cells_cohort1)%in% tumor_cells$Cell]
  tumor_scRNA         <- CreateSeuratObject(counts = tumor_exp)
  tumor_scRNA         <- NormalizeData(tumor_scRNA)
  tumor_exp_bct93                   <- as.data.frame(tumor_scRNA@assays$RNA@data)
  tumor_exp_bct93$SYMBOL            <- rownames(tumor_exp_bct93)
  tumor_exp_bct93                   <- merge(UBS93.data$UBS93.gene.df,tumor_exp_bct93,by="SYMBOL")
  rown                               <- tumor_exp_bct93$ENSEMBL
  tumor_exp_bct93                   <- tumor_exp_bct93[,-c(1:3)] %>% as.matrix()
  rownames(tumor_exp_bct93)         <- rown
  predictor.res                     <- breast.cancer.predictor(expr.of.sample = tumor_exp_bct93,
                                                expr.of.centroid = UBS93.data$UBS93.centroid,
                                                marker.gene = UBS93.data$UBS93.gene.df$ENSEMBL,
                                                HER2.amp.signature.genes = UBS93.data$HER2.amp.signature.genes,
                                                ssgsea.cutoff = UBS93.data$ssgsea.cutoff)
  
  Subtype                    <- table(predictor.res$subtype) %>% as.data.frame() 
  if (is.null(Bassez_subtype_result)) {
    Bassez_subtype_result <- Subtype
  } else {
    Bassez_subtype_result <- merge(Bassez_subtype_result, Subtype,by="subtype",all=T)
  }
}
colnames(Bassez_subtype_result)[2:14]   <- patient_id
rown                                    <- as.character(Bassez_subtype_result$subtype)
Bassez_subtype_result                   <- Bassez_subtype_result[,-1]    
rownames(Bassez_subtype_result)         <- rown
Bassez_subtype_result                   <- as.matrix(Bassez_subtype_result)  %>% t()  %>% as.data.frame()
Bassez_subtype_result$total_numbers     <- apply(Bassez_subtype_result[,1:5],1,function(x){sum(x,na.rm =T)})
Bassez_subtype_result$Basal_ratio       <- Bassez_subtype_result$Basal / Bassez_subtype_result$total_numbers
Bassez_subtype_result$Claudin_low_ratio <- Bassez_subtype_result$Claudin_low / Bassez_subtype_result$total_numbers
Bassez_subtype_result$HER2_amp_ratio    <- Bassez_subtype_result$HER2_amp / Bassez_subtype_result$total_numbers
Bassez_subtype_result$Luminal_ratio     <- Bassez_subtype_result$Luminal / Bassez_subtype_result$total_numbers
Bassez_subtype_ratio                    <- Bassez_subtype_result[,7:10]
write.csv(Bassez_subtype_ratio,file = "Pro_TNBC/paper/data/results/section_4/TNBC/Bassez_subtype_ratio.csv")


####*cohort2####
counts_cells_cohort2                   <- readRDS("~/Pro_TNBC/data/scRNASeq/Data_Bassez2021_Breast/1867-counts_cells_cohort2.rds")
counts_cells_cohort2                   <- as.matrix(counts_cells_cohort2)
BIOKEY_metaData_cohort2_web            <- read_csv("Pro_TNBC/data/scRNASeq/Data_Bassez2021_Breast/1871-BIOKEY_metaData_cohort2_web.csv")
BIOKEY_metaData_cohort2_TNBC           <- subset(BIOKEY_metaData_cohort2_web,BC_type=="TNBC")
table(BIOKEY_metaData_cohort2_TNBC$patient_id)

patient_id                             <- c("BIOKEY_33" , "BIOKEY_35", "BIOKEY_36", "BIOKEY_39", "BIOKEY_41")
Bassez_subtype_result_2                <- NULL
for (i in 1:length(patient_id)) {
  A                   <- patient_id[i] %>% as.character()
  tumor_cells         <- subset(BIOKEY_metaData_cohort2_TNBC,patient_id== A &cellType=="Cancer_cell"& timepoint=="Pre")
  tumor_exp           <- counts_cells_cohort2[,colnames(counts_cells_cohort2)%in% tumor_cells$Cell]
  tumor_scRNA         <- CreateSeuratObject(counts = tumor_exp)
  tumor_scRNA         <- NormalizeData(tumor_scRNA)
  tumor_exp_bct93                   <- as.data.frame(tumor_scRNA@assays$RNA@data)
  tumor_exp_bct93$SYMBOL            <- rownames(tumor_exp_bct93)
  tumor_exp_bct93                   <- merge(UBS93.data$UBS93.gene.df,tumor_exp_bct93,by="SYMBOL")
  rown                                       <- tumor_exp_bct93$ENSEMBL
  tumor_exp_bct93                   <- tumor_exp_bct93[,-c(1:3)] %>% as.matrix()
  rownames(tumor_exp_bct93)         <- rown
  predictor.res       <- breast.cancer.predictor(expr.of.sample = tumor_exp_bct93,
                                                 expr.of.centroid = UBS93.data$UBS93.centroid,
                                                 marker.gene = UBS93.data$UBS93.gene.df$ENSEMBL,
                                                 HER2.amp.signature.genes = UBS93.data$HER2.amp.signature.genes,
                                                 ssgsea.cutoff = UBS93.data$ssgsea.cutoff)
  
  Subtype                    <- table(predictor.res$subtype) %>% as.data.frame() 
  if (is.null(Bassez_subtype_result_2)) {
    Bassez_subtype_result_2 <- Subtype
  } else {
    Bassez_subtype_result_2 <- merge(Bassez_subtype_result_2, Subtype,by="subtype",all=T)
  }
}
colnames(Bassez_subtype_result_2)[2:6] <- patient_id
rown                                   <- as.character(Bassez_subtype_result_2$subtype)
Bassez_subtype_result_2                <- Bassez_subtype_result_2[,-1]    
rownames(Bassez_subtype_result_2)      <- rown
Bassez_subtype_result_2                <- as.matrix(Bassez_subtype_result_2)  %>% t()  %>% as.data.frame()
Bassez_subtype_result_2$total_numbers  <- apply(Bassez_subtype_result_2[,1:4],1,function(x){sum(x,na.rm =T)})
Bassez_subtype_result_2$Basal_ratio    <- Bassez_subtype_result_2$Basal / Bassez_subtype_result_2$total_numbers
Bassez_subtype_result_2$Claudin_low_ratio <- Bassez_subtype_result_2$Claudin_low / Bassez_subtype_result_2$total_numbers
Bassez_subtype_result_2$HER2_amp_ratio    <- Bassez_subtype_result_2$HER2_amp / Bassez_subtype_result_2$total_numbers
Bassez_subtype_result_2$Luminal_ratio     <- Bassez_subtype_result_2$Luminal / Bassez_subtype_result_2$total_numbers
Bassez_subtype_ratio_2                    <- Bassez_subtype_result_2[,6:9]
Bassez_subtype_ratio                      <- rbind(Bassez_subtype_ratio,Bassez_subtype_ratio_2)
write.csv(Bassez_subtype_ratio,file = "Pro_TNBC/paper/data/results/section_4/TNBC/Bassez_subtype_ratio.csv")


write.csv(Bassez_subtype_result,file = "Pro_TNBC/paper/data/results/section_4/TNBC/Bassez_subtype_result.csv")
write.csv(Bassez_subtype_result_2,file = "Pro_TNBC/paper/data/results/section_4/TNBC/Bassez_subtype_result_2.csv")


####*the correlation coefficients between claudin-low cells and centroid####
#cohort1
for (i in 1:length(patient_id)) {
  A                   <- patient_id[i] %>% as.character()
  tumor_cells         <- subset(BIOKEY_metaData_cohort1_TNBC,patient_id== A &cellType=="Cancer_cell"& timepoint=="Pre")
  tumor_exp           <- counts_cells_cohort1[,colnames(counts_cells_cohort1)%in% tumor_cells$Cell]
  tumor_scRNA         <- CreateSeuratObject(counts = tumor_exp)
  tumor_scRNA         <- NormalizeData(tumor_scRNA)
  tumor_exp_bct93                   <- as.data.frame(tumor_scRNA@assays$RNA@data)
  tumor_exp_bct93$SYMBOL            <- rownames(tumor_exp_bct93)
  tumor_exp_bct93                   <- merge(UBS93.data$UBS93.gene.df,tumor_exp_bct93,by="SYMBOL")
  rown                                       <- tumor_exp_bct93$ENSEMBL
  tumor_exp_bct93                   <- tumor_exp_bct93[,-c(1:3)] %>% as.matrix()
  rownames(tumor_exp_bct93)         <- rown
  predictor.res                     <- breast.cancer.predictor(expr.of.sample = tumor_exp_bct93,
                                                 expr.of.centroid = UBS93.data$UBS93.centroid,
                                                 marker.gene = UBS93.data$UBS93.gene.df$ENSEMBL,
                                                 HER2.amp.signature.genes = UBS93.data$HER2.amp.signature.genes,
                                                 ssgsea.cutoff = UBS93.data$ssgsea.cutoff)
  
  correlation          <- predictor.res$cor.matrix[,2,drop=F] %>% as.data.frame()
  correlation$cell.id  <- rownames(correlation)
  subtype              <- predictor.res$subtype 
  subtype$cell.id      <- rownames(subtype)
  claudin_low_id       <- subtype[subtype=="Claudin_low",]$cell.id
  claudin_low_cor      <- correlation[correlation$cell.id %in% claudin_low_id,]
  claudin_low_cor$patient.id   <- rep(A,length(claudin_low_id)) 
  if (is.null(claudin_low_cormat)) {
    claudin_low_cormat <- claudin_low_cor
  } else {
    claudin_low_cormat <- rbind(claudin_low_cormat, claudin_low_cor)
  }
}
table(claudin_low_cormat$patient.id)

#cohort2
for (i in 1:length(patient_id)) {
  A                                 <- patient_id[i] %>% as.character()
  tumor_cells                       <- subset(BIOKEY_metaData_cohort2_TNBC,patient_id== A &cellType=="Cancer_cell"& timepoint=="Pre")
  tumor_exp                         <- counts_cells_cohort2[,colnames(counts_cells_cohort2)%in% tumor_cells$Cell]
  tumor_scRNA                       <- CreateSeuratObject(counts = tumor_exp)
  tumor_scRNA                       <- NormalizeData(tumor_scRNA)
  tumor_exp_bct93                   <- as.data.frame(tumor_scRNA@assays$RNA@data)
  tumor_exp_bct93$SYMBOL            <- rownames(tumor_exp_bct93)
  tumor_exp_bct93                   <- merge(UBS93.data$UBS93.gene.df,tumor_exp_bct93,by="SYMBOL")
  rown                              <- tumor_exp_bct93$ENSEMBL
  tumor_exp_bct93                   <- tumor_exp_bct93[,-c(1:3)] %>% as.matrix()
  rownames(tumor_exp_bct93)         <- rown
  predictor.res       <- breast.cancer.predictor(expr.of.sample = tumor_exp_bct93,
                                                 expr.of.centroid = UBS93.data$UBS93.centroid,
                                                 marker.gene = UBS93.data$UBS93.gene.df$ENSEMBL,
                                                 HER2.amp.signature.genes = UBS93.data$HER2.amp.signature.genes,
                                                 ssgsea.cutoff = UBS93.data$ssgsea.cutoff)
  
  correlation          <- predictor.res$cor.matrix[,2,drop=F] %>% as.data.frame()
  correlation$cell.id  <- rownames(correlation)
  subtype              <- predictor.res$subtype 
  subtype$cell.id      <- rownames(subtype)
  claudin_low_id       <- subtype[subtype=="Claudin_low",]$cell.id
  claudin_low_cor      <- correlation[correlation$cell.id %in% claudin_low_id,]
  claudin_low_cor$patient.id   <- rep(A,length(claudin_low_id)) 
  if (is.null(claudin_low_cormat)) {
    claudin_low_cormat <- claudin_low_cor
  } else {
    claudin_low_cormat <- rbind(claudin_low_cormat, claudin_low_cor)
  }
}

Bassez_CL_cor                    <- claudin_low_cormat
ggplot(Bassez_CL_cor,aes(x=patient.id,y=Claudin_low))+geom_boxplot()
table(Bassez_CL_cor$patient.id)
Bassez_CL_cor                    <- Bassez_CL_cor[,-2]
TNBC_CL_cor                      <- bind_rows(Bassez_CL_cor,Qian_cl_cor,GSE118389_cl_cor,GSE176078_TNBC_CLcormat,GSE161529_TNBC_CLcormat)
save(TNBC_CL_cor,file = "Pro_TNBC/paper/data/results/section_4/TNBC/TNBC_CL_cor.RData")

