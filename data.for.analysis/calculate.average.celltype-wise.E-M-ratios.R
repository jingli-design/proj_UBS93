############################################################################################################
#### R code to calculate the average celltype-wise E-M-ratios for all genes by using single-cell RNAseq data of breast cancer patients
############################################################################################################

library(Matrix)
library(Seurat)
library(dplyr)
library(readr)

####1.create seurat object####
dir="Pro_TNBC/data/scRNASeq/26/Wu_etal_2021_BRCA_scRNASeq/"
list.files(dir)
counts                            <- Read10X(
  data.dir=dir,
  gene.column = 1,
  cell.column = 2,
  unique.features = TRUE,
  strip.suffix = FALSE
)
class(counts)
GSE176078_scRNA                   <- CreateSeuratObject(counts = counts,min.cell=3,min.features=200)
save(GSE176078_scRNA,file = "Pro_TNBC/paper/data/results/section_1/GSE176078_scRNA.RData")
####2.compute the E-M-ratio ####
####*CID3586(HER2+/ER+)####
CID3586_scRNA                     <- subset(GSE176078_scRNA,orig.ident=="CID3586")
CID3586_scRNA[["percent.mt"]]     <- PercentageFeatureSet(CID3586_scRNA,pattern = "^MT-")
head(CID3586_scRNA@meta.data,5)
VlnPlot(CID3586_scRNA,features = c("nCount_RNA","nFeature_RNA","percent.mt"),ncol=3)
CID3586_scRNA                     <- subset(CID3586_scRNA,nFeature_RNA >0 & nFeature_RNA<4000&percent.mt <12 )

#normalized
CID3586_scRNA                     <- NormalizeData(CID3586_scRNA,normalization.method = "RC")
CID3586_exprmat_CPM               <- as.matrix(CID3586_scRNA@assays$RNA@data)#data after normalizing
metadata                          <- read_csv("Pro_TNBC/data/scRNASeq/26/Wu_etal_2021_BRCA_scRNASeq/metadata.csv")
CID3586_metadata                  <- metadata[metadata$orig.ident=="CID3586",]
CID3586                           <- CID3586_metadata[,c(1,9)]
colnames(CID3586)[1]              <- "cell.id"
CID3586_exprdf                    <- t(CID3586_exprmat_CPM) %>% as.data.frame()
CID3586_exprdf$cell.id            <- rownames(CID3586_exprdf)
CID3586_exprdf                    <- merge(CID3586,CID3586_exprdf,by="cell.id")
normal.epithelial                 <- apply(CID3586_exprdf[CID3586_exprdf$celltype_major=="Normal Epithelial",3:27721],2,mean)
T.cells                           <- apply(CID3586_exprdf[CID3586_exprdf$celltype_major=="T-cells",3:27721],2,mean)
CID3586_score                     <- data.frame(normal.epithelial,T.cells)
CID3586_score$B.cells             <- apply(CID3586_exprdf[CID3586_exprdf$celltype_major=="B-cells",3:27721],2,mean)
CID3586_score$CAFs                <- apply(CID3586_exprdf[CID3586_exprdf$celltype_major=="CAFs",3:27721],2,mean)
CID3586_score$PVL                 <- apply(CID3586_exprdf[CID3586_exprdf$celltype_major=="PVL",3:27721],2,mean)
CID3586_score$Endothelial         <- apply(CID3586_exprdf[CID3586_exprdf$celltype_major=="Endothelial",3:27721],2,mean)
CID3586_score$Myeloid             <- apply(CID3586_exprdf[CID3586_exprdf$celltype_major=="Myeloid",3:27721],2,mean)
CID3586_score$TME                 <- (CID3586_score$T.cells*4587 + 
                                   CID3586_score$B.cells*314 +
                                   CID3586_score$CAFs*175 + 
                                   CID3586_score$PVL*21 + 
                                   CID3586_score$Endothelial*143 + 
                                   CID3586_score$Myeloid*198)/5546

pseudocount                       <- 0.1
CID3586_score$CID3586_S           <- (CID3586_score$normal.epithelial + pseudocount) /(CID3586_score$TME + pseudocount)
write.csv(CID3586_score,file="Pro_TNBC/paper/data/results/section_1/CID3586_normalized_score.csv")

####*CID3921(HER2+)####
CID3921_scRNA                     <- subset(GSE176078_scRNA,orig.ident=="CID3921")
CID3921_scRNA[["percent.mt"]]     <- PercentageFeatureSet(CID3921_scRNA,pattern = "^MT-")
head(CID3921_scRNA@meta.data,5)
VlnPlot(CID3921_scRNA,features = c("nCount_RNA","nFeature_RNA","percent.mt"),ncol=3)
CID3921_scRNA                     <- subset(CID3921_scRNA,nFeature_RNA >0 & nFeature_RNA<8000&percent.mt <20 )

#normalized
CID3921_scRNA                     <- NormalizeData(CID3921_scRNA,normalization.method = "RC")
CID3921_exprmat_CPM               <- as.matrix(CID3921_scRNA@assays$RNA@data)#data after normalizing
metadata                          <- read_csv("Pro_TNBC/data/scRNASeq/26/Wu_etal_2021_BRCA_scRNASeq/metadata.csv")
CID3921_metadata                  <- metadata[metadata$orig.ident=="CID3921",]
CID3921                           <- CID3921_metadata[,c(1,9)]
colnames(CID3921)[1]              <- "cell.id"
CID3921_exprdf                    <- t(CID3921_exprmat_CPM) %>% as.data.frame()
CID3921_exprdf$cell.id            <- rownames(CID3921_exprdf)
CID3921_exprdf                    <- merge(CID3921,CID3921_exprdf,by="cell.id")

cancer.epithelial            <- apply(CID3921_exprdf[CID3921_exprdf$celltype_major=="Cancer Epithelial",3:27721],2,mean)
T.cells                      <- apply(CID3921_exprdf[CID3921_exprdf$celltype_major=="T-cells",3:27721],2,mean)
CID3921_score                <- data.frame(cancer.epithelial,T.cells)
CID3921_score$B.cells        <- apply(CID3921_exprdf[CID3921_exprdf$celltype_major=="B-cells",3:27721],2,mean)
CID3921_score$CAFs           <- apply(CID3921_exprdf[CID3921_exprdf$celltype_major=="CAFs",3:27721],2,mean)
CID3921_score$PVL            <- apply(CID3921_exprdf[CID3921_exprdf$celltype_major=="PVL",3:27721],2,mean)
CID3921_score$Endothelial    <- apply(CID3921_exprdf[CID3921_exprdf$celltype_major=="Endothelial",3:27721],2,mean)
CID3921_score$Myeloid        <- apply(CID3921_exprdf[CID3921_exprdf$celltype_major=="Myeloid",3:27721],2,mean)
CID3921_score$Plasmablasts   <- apply(CID3921_exprdf[CID3921_exprdf$celltype_major=="Plasmablasts",3:27721],2,mean)
CID3921_score$TME            <- (CID3921_score$T.cells*1473 + 
                                   CID3921_score$B.cells*162 +
                                   CID3921_score$CAFs*106 + 
                                   CID3921_score$PVL*72 + 
                                   CID3921_score$Endothelial*210 + 
                                   CID3921_score$Myeloid*385+
                                   CID3921_score$Plasmablasts*175)/2583

CID3921_score$CID3921_S      <- (CID3921_score$cancer.epithelial + pseudocount)/(CID3921_score$TME + pseudocount)
write.csv(CID3921_score,file = "Pro_TNBC/paper/data/results/section_1/CID3921_normalized_score.csv")


####*CID4495(TNBC)####
CID4495_scRNA                     <- subset(GSE176078_scRNA,orig.ident=="CID4495")
CID4495_scRNA[["percent.mt"]]     <- PercentageFeatureSet(CID4495_scRNA,pattern = "^MT-")
head(CID4495_scRNA@meta.data,5)
VlnPlot(CID4495_scRNA,features = c("nCount_RNA","nFeature_RNA","percent.mt"),ncol=3)
CID4495_scRNA                     <- subset(CID4495_scRNA,nFeature_RNA >0 & nFeature_RNA<7000&percent.mt <20 )

#normalized
CID4495_scRNA                     <- NormalizeData(CID4495_scRNA,normalization.method = "RC")
CID4495_exprmat_CPM               <- as.matrix(CID4495_scRNA@assays$RNA@data)#data after normalizing

##compute the score in CID4495##
metadata                <- read_csv("Pro_TNBC/data/scRNASeq/26/Wu_etal_2021_BRCA_scRNASeq/metadata.csv")
CID4495_metadata        <- metadata[metadata$orig.ident=="CID4495",]
CID4495                 <- CID4495_metadata[,c(1,9)]
colnames(CID4495)[1]    <- "cell.id"
CID4495_exprdf          <- t(CID4495_exprmat_CPM) %>% as.data.frame()
CID4495_exprdf$cell.id  <- rownames(CID4495_exprdf)
CID4495_exprdf          <- merge(CID4495,CID4495_exprdf,by="cell.id")

Epithelial              <- apply(CID4495_exprdf[CID4495_exprdf$celltype_major=="Cancer Epithelial",3:27721],2,mean)
T.cells                 <- apply(CID4495_exprdf[CID4495_exprdf$celltype_major=="T-cells",3:27721],2,mean)
CID4495_score           <- data.frame(Epithelial,T.cells)
CID4495_score$CAFs      <- apply(CID4495_exprdf[CID4495_exprdf$celltype_major=="CAFs",3:27721],2,mean)
CID4495_score$PVL       <- apply(CID4495_exprdf[CID4495_exprdf$celltype_major=="PVL",3:27721],2,mean)
CID4495_score$Endothelial    <- apply(CID4495_exprdf[CID4495_exprdf$celltype_major=="Endothelial",3:27721],2,mean)
CID4495_score$Myeloid        <- apply(CID4495_exprdf[CID4495_exprdf$celltype_major=="Myeloid",3:27721],2,mean)
CID4495_score$B_cells        <- apply(CID4495_exprdf[CID4495_exprdf$celltype_major=="B-cells",3:27721],2,mean)
CID4495_score$Plasmablasts   <- apply(CID4495_exprdf[CID4495_exprdf$celltype_major=="Plasmablasts",3:27721],2,mean)          
CID4495_score$TME            <- (CID4495_score$T.cells*3504 + 
                                   CID4495_score$CAFs*232 + 
                                   CID4495_score$PVL*191 + 
                                   CID4495_score$Endothelial*184 + 
                                   CID4495_score$Myeloid*897+
                                   CID4495_score$B_cells*773 +
                                   CID4495_score$Plasmablasts*1020 )/6801

CID4495_score$"CID4495_S"          <- (CID4495_score$Epithelial + pseudocount) /(CID4495_score$TME + pseudocount)
write.csv(CID4495_score,file = "Pro_TNBC/paper/data/results/section_1/CID4495_normalized_score.csv")

####*CID3948(ER+)####
CID3948_scRNA                     <- subset(GSE176078_scRNA,orig.ident=="CID3948")
CID3948_scRNA[["percent.mt"]]     <- PercentageFeatureSet(CID3948_scRNA,pattern = "^MT-")
head(CID3948_scRNA@meta.data,5)
VlnPlot(CID3948_scRNA,features = c("nCount_RNA","nFeature_RNA","percent.mt"),ncol=3)
CID3948_scRNA                     <- subset(CID3948_scRNA,nFeature_RNA >0 & nFeature_RNA<4000&percent.mt <20 )

#normalized
CID3948_scRNA                     <- NormalizeData(CID3948_scRNA,normalization.method = "RC")
CID3948_exprmat_CPM                   <- as.matrix(CID3948_scRNA@assays$RNA@data)#data after normalizing

##compute the score in CID3948##
metadata                <- read_csv("Pro_TNBC/data/scRNASeq/26/Wu_etal_2021_BRCA_scRNASeq/metadata.csv")
CID3948_metadata        <- metadata[metadata$orig.ident=="CID3948",]
CID3948                 <- CID3948_metadata[,c(1,9)]
colnames(CID3948)[1]    <- "cell.id"
CID3948_exprdf          <- t(CID3948_exprmat_CPM) %>% as.data.frame()
CID3948_exprdf$cell.id  <- rownames(CID3948_exprdf)
CID3948_exprdf          <- merge(CID3948,CID3948_exprdf,by="cell.id")
Epithelial              <- apply(CID3948_exprdf[CID3948_exprdf$celltype_major=="Cancer Epithelial",3:27721],2,mean)
T.cells                 <- apply(CID3948_exprdf[CID3948_exprdf$celltype_major=="T-cells",3:27721],2,mean)
CID3948_score           <- data.frame(Epithelial,T.cells)
CID3948_score$CAFs      <- apply(CID3948_exprdf[CID3948_exprdf$celltype_major=="CAFs",3:27721],2,mean)
CID3948_score$PVL       <- apply(CID3948_exprdf[CID3948_exprdf$celltype_major=="PVL",3:27721],2,mean)
CID3948_score$Endothelial    <- apply(CID3948_exprdf[CID3948_exprdf$celltype_major=="Endothelial",3:27721],2,mean)
CID3948_score$Myeloid        <- apply(CID3948_exprdf[CID3948_exprdf$celltype_major=="Myeloid",3:27721],2,mean)
CID3948_score$B_cells        <- apply(CID3948_exprdf[CID3948_exprdf$celltype_major=="B-cells",3:27721],2,mean)
CID3948_score$Plasmablasts   <- apply(CID3948_exprdf[CID3948_exprdf$celltype_major=="Plasmablasts",3:27721],2,mean)          
CID3948_score$TME            <- (CID3948_score$T.cells*1465 + 
                                   CID3948_score$CAFs*15 + 
                                   CID3948_score$PVL*62 + 
                                   CID3948_score$Endothelial*85 + 
                                   CID3948_score$Myeloid*122+
                                   CID3948_score$B_cells*85 +
                                   CID3948_score$Plasmablasts*232 )/2066

CID3948_score$"CID3948_S"     <- (CID3948_score$Epithelial + pseudocount)/(CID3948_score$TME + pseudocount)
write.csv(CID3948_score,file = "Pro_TNBC/paper/data/results/section_1/CID3948_normalized_score.csv")

####*CID4066(ER+/HER2+)####
CID4066_scRNA                     <- subset(GSE176078_scRNA,orig.ident=="CID4066")
CID4066_scRNA[["percent.mt"]]     <- PercentageFeatureSet(CID4066_scRNA,pattern = "^MT-")
head(CID4066_scRNA@meta.data,5)
VlnPlot(CID4066_scRNA,features = c("nCount_RNA","nFeature_RNA","percent.mt"),ncol=3)
CID4066_scRNA                     <- subset(CID4066_scRNA,nFeature_RNA >0 & nFeature_RNA<7000&percent.mt <20 )

#normalized
CID4066_scRNA                     <- NormalizeData(CID4066_scRNA,normalization.method = "RC")
CID4066_exprmat_CPM               <- as.matrix(CID4066_scRNA@assays$RNA@data)#data after normalizing

##compute the score in CID4066##
metadata                <- read_csv("Pro_TNBC/data/scRNASeq/26/Wu_etal_2021_BRCA_scRNASeq/metadata.csv")
CID4066_metadata        <- metadata[metadata$orig.ident=="CID4066",]
CID4066                 <- CID4066_metadata[,c(1,9)]
colnames(CID4066)[1]    <- "cell.id"
CID4066_exprdf          <- t(CID4066_exprmat_CPM) %>% as.data.frame()
CID4066_exprdf$cell.id  <- rownames(CID4066_exprdf)
CID4066_exprdf          <- merge(CID4066,CID4066_exprdf,by="cell.id")

Epithelial              <- apply(CID4066_exprdf[CID4066_exprdf$celltype_major=="Cancer Epithelial"|CID4066_exprdf$celltype_major=="Normal Epithelial",3:27721],2,mean)
T.cells                 <- apply(CID4066_exprdf[CID4066_exprdf$celltype_major=="T-cells",3:27721],2,mean)
CID4066_score           <- data.frame(Epithelial,T.cells)
CID4066_score$CAFs      <- apply(CID4066_exprdf[CID4066_exprdf$celltype_major=="CAFs",3:27721],2,mean)
CID4066_score$PVL       <- apply(CID4066_exprdf[CID4066_exprdf$celltype_major=="PVL",3:27721],2,mean)
CID4066_score$Endothelial    <- apply(CID4066_exprdf[CID4066_exprdf$celltype_major=="Endothelial",3:27721],2,mean)
CID4066_score$Myeloid        <- apply(CID4066_exprdf[CID4066_exprdf$celltype_major=="Myeloid",3:27721],2,mean)
CID4066_score$B_cells        <- apply(CID4066_exprdf[CID4066_exprdf$celltype_major=="B-cells",3:27721],2,mean)
CID4066_score$TME            <- (CID4066_score$T.cells*2171 + 
                                   CID4066_score$CAFs*923 + 
                                   CID4066_score$PVL*630 + 
                                   CID4066_score$Endothelial*535 + 
                                   CID4066_score$Myeloid*221+
                                   CID4066_score$B_cells*38 )/4518

CID4066_score$"CID4066_S"    <- (CID4066_score$Epithelial + pseudocount)/(CID4066_score$TME + pseudocount)
write.csv(CID4066_score,file = "Pro_TNBC/paper/data/results/section_1/CID4066_normalized_score.csv")

####*CID45171(HER2+)####
CID45171_scRNA                     <- subset(GSE176078_scRNA,orig.ident=="CID45171")
CID45171_scRNA[["percent.mt"]]     <- PercentageFeatureSet(CID45171_scRNA,pattern = "^MT-")
head(CID45171_scRNA@meta.data,5)
VlnPlot(CID45171_scRNA,features = c("nCount_RNA","nFeature_RNA","percent.mt"),ncol=3)
CID45171_scRNA                     <- subset(CID45171_scRNA,nFeature_RNA >0 & nFeature_RNA<7000&percent.mt <20 )

#normalized
CID45171_scRNA                     <- NormalizeData(CID45171_scRNA,normalization.method = "RC")
CID45171_exprmat_CPM               <- as.matrix(CID45171_scRNA@assays$RNA@data)#data after normalizing

##compute the score in CID45171##
metadata                 <- read_csv("Pro_TNBC/data/scRNASeq/26/Wu_etal_2021_BRCA_scRNASeq/metadata.csv")
CID45171_metadata        <- metadata[metadata$orig.ident=="CID45171",]
CID45171                 <- CID45171_metadata[,c(1,9)]
colnames(CID45171)[1]    <- "cell.id"
CID45171_exprdf          <- t(CID45171_exprmat_CPM) %>% as.data.frame()
CID45171_exprdf$cell.id  <- rownames(CID45171_exprdf)
CID45171_exprdf          <- merge(CID45171,CID45171_exprdf,by="cell.id")

Epithelial               <- apply(CID45171_exprdf[CID45171_exprdf$celltype_major=="Cancer Epithelial",3:27721],2,mean)
T.cells                  <- apply(CID45171_exprdf[CID45171_exprdf$celltype_major=="T-cells",3:27721],2,mean)
CID45171_score           <- data.frame(Epithelial,T.cells)
CID45171_score$CAFs      <- apply(CID45171_exprdf[CID45171_exprdf$celltype_major=="CAFs",3:27721],2,mean)
CID45171_score$PVL       <- apply(CID45171_exprdf[CID45171_exprdf$celltype_major=="PVL",3:27721],2,mean)
CID45171_score$Endothelial    <- apply(CID45171_exprdf[CID45171_exprdf$celltype_major=="Endothelial",3:27721],2,mean)
CID45171_score$Myeloid        <- apply(CID45171_exprdf[CID45171_exprdf$celltype_major=="Myeloid",3:27721],2,mean)
CID45171_score$B_cells        <- apply(CID45171_exprdf[CID45171_exprdf$celltype_major=="B-cells",3:27721],2,mean)
CID45171_score$TME            <- (CID45171_score$T.cells*1346 + 
                                    CID45171_score$CAFs*32 + 
                                    CID45171_score$PVL*13 + 
                                    CID45171_score$Endothelial*15 + 
                                    CID45171_score$Myeloid*172+
                                    CID45171_score$B_cells*56 )/1634

CID45171_score$"CID45171_S"         <- (CID45171_score$Epithelial + pseudocount)/(CID45171_score$TME + pseudocount)
write.csv(CID45171_score,file = "Pro_TNBC/paper/data/results/section_1/CID45171_normalized_score.csv")

####*CID44041(TNBC)####
CID44041_scRNA                     <- subset(scRNA,orig.ident=="CID44041")
CID44041_scRNA[["percent.mt"]]     <- PercentageFeatureSet(CID44041_scRNA,pattern = "^MT-")
head(CID44041_scRNA@meta.data,5)
VlnPlot(CID44041_scRNA,features = c("nCount_RNA","nFeature_RNA","percent.mt"),ncol=3)
CID44041_scRNA                     <- subset(CID44041_scRNA,nFeature_RNA >0 & nFeature_RNA<3500&percent.mt <10 )

#normalized
CID44041_scRNA                     <- NormalizeData(CID44041_scRNA,normalization.method = "RC")
CID44041_exprmat_CPM                   <- as.matrix(CID44041_scRNA@assays$RNA@data)#data after normalizing


##compute the score in CID44041##
metadata                 <- read_csv("Pro_TNBC/data/scRNASeq/26/Wu_etal_2021_BRCA_scRNASeq/metadata.csv")
CID44041_metadata        <- metadata[metadata$orig.ident=="CID44041",]
CID44041                 <- CID44041_metadata[,c(1,9)]
colnames(CID44041)[1]    <- "cell.id"
CID44041_exprdf          <- t(CID44041_exprmat_CPM) %>% as.data.frame()
CID44041_exprdf$cell.id  <- rownames(CID44041_exprdf)
CID44041_exprdf          <- merge(CID44041,CID44041_exprdf,by="cell.id")

Epithelial               <- apply(CID44041_exprdf[CID44041_exprdf$celltype_major=="Normal Epithelial",3:27721],2,mean)
T.cells                  <- apply(CID44041_exprdf[CID44041_exprdf$celltype_major=="T-cells",3:27721],2,mean)
CID44041_score           <- data.frame(Epithelial,T.cells)
CID44041_score$CAFs      <- apply(CID44041_exprdf[CID44041_exprdf$celltype_major=="CAFs",3:27721],2,mean)
CID44041_score$PVL       <- apply(CID44041_exprdf[CID44041_exprdf$celltype_major=="PVL",3:27721],2,mean)
CID44041_score$Endothelial    <- apply(CID44041_exprdf[CID44041_exprdf$celltype_major=="Endothelial",3:27721],2,mean)
CID44041_score$Myeloid        <- apply(CID44041_exprdf[CID44041_exprdf$celltype_major=="Myeloid",3:27721],2,mean)
CID44041_score$B_cells        <- apply(CID44041_exprdf[CID44041_exprdf$celltype_major=="B-cells",3:27721],2,mean)
CID44041_score$TME            <- (CID44041_score$T.cells*742 + 
                                    CID44041_score$CAFs*680 + 
                                    CID44041_score$PVL*128 + 
                                    CID44041_score$Endothelial*146 + 
                                    CID44041_score$Myeloid*105+
                                    CID44041_score$B_cells*175)/1995

CID44041_score$"CID44041_S"   <- (CID44041_score$Epithelial + pseudocount)/(CID44041_score$TME + pseudocount)
write.csv(CID44041_score,file = "Pro_TNBC/paper/data/results/section_1/CID44041_normalized_score.csv")

####*CID3963(ER+)####
CID3963_scRNA                     <- subset(GSE176078_scRNA,orig.ident=="CID3963")
CID3963_scRNA[["percent.mt"]]     <- PercentageFeatureSet(CID3963_scRNA,pattern = "^MT-")
head(CID3963_scRNA@meta.data,5)
VlnPlot(CID3963_scRNA,features = c("nCount_RNA","nFeature_RNA","percent.mt"),ncol=3)
CID3963_scRNA                     <- subset(CID3963_scRNA,nFeature_RNA >0 & nFeature_RNA<4000&percent.mt <20 )

#normalized
CID3963_scRNA                     <- NormalizeData(CID3963_scRNA,normalization.method = "RC")
CID3963_exprmat_CPM               <- as.matrix(CID3963_scRNA@assays$RNA@data)#data after normalizing

##compute the score in CID3963##
metadata                <- read_csv("Pro_TNBC/data/scRNASeq/26/Wu_etal_2021_BRCA_scRNASeq/metadata.csv")
CID3963_metadata        <- metadata[metadata$orig.ident=="CID3963",]
CID3963                 <- CID3963_metadata[,c(1,9)]
colnames(CID3963)[1]    <- "cell.id"
CID3963_exprdf          <- t(CID3963_exprmat_CPM) %>% as.data.frame()
CID3963_exprdf$cell.id  <- rownames(CID3963_exprdf)
CID3963_exprdf          <- merge(CID3963,CID3963_exprdf,by="cell.id")

Epithelial              <- apply(CID3963_exprdf[CID3963_exprdf$celltype_major=="Cancer Epithelial"|CID3963_exprdf$celltype_major=="Normal Epithelial",3:27721],2,mean)
T.cells                 <- apply(CID3963_exprdf[CID3963_exprdf$celltype_major=="T-cells",3:27721],2,mean)
CID3963_score           <- data.frame(Epithelial,T.cells)
CID3963_score$CAFs      <- apply(CID3963_exprdf[CID3963_exprdf$celltype_major=="CAFs",3:27721],2,mean)
CID3963_score$PVL       <- apply(CID3963_exprdf[CID3963_exprdf$celltype_major=="PVL",3:27721],2,mean)
CID3963_score$Endothelial    <- apply(CID3963_exprdf[CID3963_exprdf$celltype_major=="Endothelial",3:27721],2,mean)
CID3963_score$Myeloid        <- apply(CID3963_exprdf[CID3963_exprdf$celltype_major=="Myeloid",3:27721],2,mean)
CID3963_score$TME            <- (CID3963_score$T.cells*2672 + 
                                   CID3963_score$CAFs*23 + 
                                   CID3963_score$PVL*28 + 
                                   CID3963_score$Endothelial*102 + 
                                   CID3963_score$Myeloid*474)/3299

CID3963_score$"CID3963_S"    <- (CID3963_score$Epithelia + pseudocount)/(CID3963_score$TME + pseudocount)
write.csv(CID3963_score,file = "Pro_TNBC/paper/data/results/section_1//CID3963_normalized_score.csv")


####3.compute the average celltype-wise E-M-ratios####
CID44041_normalized_score <- read_csv("Pro_TNBC/paper/data/results/section_1/CID44041_normalized_score.csv")
CID4495_normalized_score <- read_csv("Pro_TNBC/paper/data/results/section_1/CID4495_normalized_score.csv")
CID3948_normalized_score <- read_csv("Pro_TNBC/paper/data/results/section_1/CID3948_normalized_score.csv")
CID3963_normalized_score <- read_csv("Pro_TNBC/paper/data/results/section_1/CID3963_normalized_score.csv")
CID4066_normalized_score <- read_csv("Pro_TNBC/paper/data/results/section_1/CID4066_normalized_score.csv")
CID45171_normalized_score <- read_csv("Pro_TNBC/paper/data/results/section_1/CID45171_normalized_score.csv")
CID3586_normalized_score <- read_csv("Pro_TNBC/paper/data/results/section_1/CID3586_normalized_score.csv")
CID3921_normalized_score <- read_csv("Pro_TNBC/paper/data/results/section_1/CID3921_normalized_score.csv")

tnbc_a                 <- CID44041_normalized_score[,c(1,10)]
colnames(tnbc_a)[1]    <- "gene"
tnbc_b                 <- CID4495_normalized_score[,c(1,11)]
colnames(tnbc_b)[1]    <- "gene"
tnbc                   <- merge(tnbc_a,tnbc_b,by="gene")
tnbc$tnbc_mean         <- (tnbc$CID44041_S+tnbc$CID4495_S)/2
nontnbc_a              <- CID3948_normalized_score[,c(1,11)]
colnames(nontnbc_a)[1] <- "gene"
nontnbc_b              <- CID3963_normalized_score[,c(1,9)]
colnames(nontnbc_b)[1] <- "gene"
nontnbc_c              <- CID4066_normalized_score[,c(1,10)]
colnames(nontnbc_c)[1] <- "gene"
nontnbc_d              <- CID45171_normalized_score[,c(1,10)]
colnames(nontnbc_d)[1] <- "gene"
nontnbc_e              <- CID3586_normalized_score[,c(1,10)]
colnames(nontnbc_e)[1] <- "gene"
nontnbc_f              <- CID3921_normalized_score[,c(1,11)]
colnames(nontnbc_f)[1] <- "gene"
nontnbc                <- merge(nontnbc_a,nontnbc_b,by="gene")
nontnbc                <- merge(nontnbc,nontnbc_c,by="gene")
nontnbc                <- merge(nontnbc,nontnbc_d,by="gene")
nontnbc                <- merge(nontnbc,nontnbc_e,by="gene")
nontnbc                <- merge(nontnbc,nontnbc_f,by="gene")
nontnbc$nontnbc_mean   <- (nontnbc$CID3948_S+
                             nontnbc$CID3963_S+
                             nontnbc$CID4066_S+
                             nontnbc$CID45171_S+
                             nontnbc$CID3586_S+
                             nontnbc$CID3921_S)/6
mean_score             <- merge(nontnbc,tnbc,by="gene")
mean_score$S_mean      <- (mean_score$nontnbc_mean+mean_score$tnbc_mean)/2
write.csv(mean_score,file="Pro_TNBC/paper/data/results/section_1/mean_score.csv")

