
library(readr)
library(ggplot2)
library(dplyr)
library(org.Hs.eg.db)
library(clusterProfiler)
library(Seurat)
source('Pro_TNBC/paper/code/ggplot.style.R')

load("Pro_TNBC/data/CCLE/CCLE.RData")
load("Pro_TNBC/paper/data/results/section_1/UBS93.gene.df.RData")
load("Pro_TNBC/output/data/scRNASeq/26_sample/GSE176078_scRNA.RData")
load("Pro_TNBC/output/data/scRNASeq/26_sample/top100.gene/top100_gene.RData")
######################################################################################################################
####  Fig 1:Development of a uniform gene panel for breast cancer subtyping.###
######################################################################################################################


#### fig 1b: Box plot of the mean values of UBS93 gene in different cell types of different clinical patients.####
####*CID3963(ER+)####
CID3963_UBS93genelist_exprmat         <- CID3963_exprmat_CPM[rownames(CID3963_exprmat_CPM) %in% UBS93.gene.df$SYMBOL,]
gene_mean                             <- apply(CID3963_UBS93genelist_exprmat,2,mean)
CID3963_UBS93genelistmean_df          <- as.data.frame(gene_mean)
CID3963_UBS93genelistmean_df$cell.id  <- rownames(CID3963_UBS93genelistmean_df)
library(readr)
metadata                              <- read_csv("Pro_TNBC/data/scRNASeq/26/Wu_etal_2021_BRCA_scRNASeq/metadata.csv")
CID3963_metadata                      <- metadata[metadata$orig.ident=="CID3963",]
CID3963                               <- CID3963_metadata[,c(1,9)]
colnames(CID3963)[1]                  <- "cell.id"
CID3963_UBS93genelistmean_df          <- merge(CID3963_UBS93genelistmean_df,CID3963,by="cell.id")
write.csv(CID3963_UBS93genelistmean_df,file = "Pro_TNBC/output/data/scRNASeq/26_sample/UBS93.gene/CID3963_UBS93genelistmean_df.csv")
CID3963_UBS93genelistmean_df          <- read_csv("Pro_TNBC/output/data/scRNASeq/26_sample/UBS93.gene/CID3963_UBS93genelistmean_df.csv")
CID3963_UBS93genelistmean_df          <- subset(CID3963_UBS93genelistmean_df,celltype_major!="Normal Epithelial")
CID3963_UBS93genelistmean_df$celltype_major    <- factor(CID3963_UBS93genelistmean_df$celltype_major,levels = c("Cancer Epithelial","Endothelial","CAFs","PVL","Myeloid","T-cells"))
CID3963_UBS93genelistmean_df$subtype  <- rep("ER+",length(CID3963_UBS93genelistmean_df$cell.id)) 
write.csv(CID3963_UBS93genelistmean_df,file = "Pro_TNBC/paper/data/results/section_1/CID3963_UBS93genelistmean_df.csv")

####*CID4066(ER+/HER2+)####
CID4066_UBS93genelist_exprmat         <- CID4066_exprmat_CPM[rownames(CID4066_exprmat_CPM) %in% UBS93.gene.df$SYMBOL,]
gene_mean                             <- apply(CID4066_UBS93genelist_exprmat,2,mean)
CID4066_UBS93genelistmean_df          <- as.data.frame(gene_mean)
CID4066_UBS93genelistmean_df$cell.id  <- rownames(CID4066_UBS93genelistmean_df)
library(readr)
metadata                <- read_csv("Pro_TNBC/data/scRNASeq/26/Wu_etal_2021_BRCA_scRNASeq/metadata.csv")
CID4066_metadata        <- metadata[metadata$orig.ident=="CID4066",]
CID4066                 <- CID4066_metadata[,c(1,9)]
colnames(CID4066)[1]    <- "cell.id"
CID4066_UBS93genelistmean_df    <- merge(CID4066_UBS93genelistmean_df,CID4066,by="cell.id")
write.csv(CID4066_UBS93genelistmean_df,file = "Pro_TNBC/output/data/scRNASeq/26_sample/UBS93.gene/CID4066_UBS93genelistmean_df.csv")

CID4066_UBS93genelistmean_df     <- read_csv("Pro_TNBC/output/data/scRNASeq/26_sample/UBS93.gene/CID4066_UBS93genelistmean_df.csv")
CID4066_UBS93genelistmean_df     <- CID4066_UBS93genelistmean_df %>% mutate(gene_panel="UBS93")
CID4066_UBS93genelistmean_df     <- subset(CID4066_UBS93genelistmean_df,celltype_major!="Normal Epithelial")
CID4066_UBS93genelistmean_df$celltype_major     <- factor(CID4066_UBS93genelistmean_df$celltype_major,levels = c("Cancer Epithelial","CAFs","Endothelial","PVL","B-cells","Myeloid","T-cells"))
CID4066_UBS93genelistmean_df$subtype  <- rep("ER+/HER2+",length(CID4066_UBS93genelistmean_df$cell.id))
write.csv(CID4066_UBS93genelistmean_df,file = "Pro_TNBC/paper/data/results/section_1/CID4066_UBS93genelistmean_df.csv")


####*CID4495(TNBC)####
CID4495_UBS93genelist_exprmat         <- CID4495_exprmat_CPM[rownames(CID4495_exprmat_CPM) %in% UBS93.gene.df$SYMBOL,]
gene_mean                             <- apply(CID4495_UBS93genelist_exprmat,2,mean)
CID4495_UBS93genelistmean_df          <- as.data.frame(gene_mean)
CID4495_UBS93genelistmean_df$cell.id  <- rownames(CID4495_UBS93genelistmean_df)
library(readr)
metadata                <- read_csv("Pro_TNBC/data/scRNASeq/26/Wu_etal_2021_BRCA_scRNASeq/metadata.csv")
CID4495_metadata        <- metadata[metadata$orig.ident=="CID4495",]
CID4495                 <- CID4495_metadata[,c(1,9)]
colnames(CID4495)[1]    <- "cell.id"
CID4495_UBS93genelistmean_df           <- merge(CID4495_UBS93genelistmean_df,CID4495,by="cell.id")
write.csv(CID4495_UBS93genelistmean_df,file = "Pro_TNBC/output/data/scRNASeq/26_sample/UBS93.gene/CID4495_UBS93genelistmean_df.csv")

CID4495_UBS93genelistmean_df   <- read_csv("Pro_TNBC/output/data/scRNASeq/26_sample/UBS93.gene/CID4495_UBS93genelistmean_df.csv")
CID4495_UBS93genelistmean_df   <- CID4495_UBS93genelistmean_df %>% mutate(gene_panel="UBS93")
CID4495_UBS93genelistmean_df$subtype  <- rep("TNBC",length(CID4495_UBS93genelistmean_df$cell.id))
write.csv(CID4495_UBS93genelistmean_df,file = "Pro_TNBC/paper/data/results/section_1/CID4495_UBS93genelistmean_df.csv")

####*CID45171(HER2+)####
CID45171_UBS93genelist_exprmat         <- CID45171_exprmat_CPM[rownames(CID45171_exprmat_CPM) %in% UBS93.gene.df$SYMBOL,]
gene_mean                              <- apply(CID45171_UBS93genelist_exprmat,2,mean)
CID45171_UBS93genelistmean_df          <- as.data.frame(gene_mean)
CID45171_UBS93genelistmean_df$cell.id  <- rownames(CID45171_UBS93genelistmean_df)
library(readr)
metadata                 <- read_csv("Pro_TNBC/data/scRNASeq/26/Wu_etal_2021_BRCA_scRNASeq/metadata.csv")
CID45171_metadata        <- metadata[metadata$orig.ident=="CID45171",]
CID45171                 <- CID45171_metadata[,c(1,9)]
colnames(CID45171)[1]    <- "cell.id"
CID45171_UBS93genelistmean_df           <- merge(CID45171_UBS93genelistmean_df,CID45171,by="cell.id")
write.csv(CID45171_UBS93genelistmean_df,file = "Pro_TNBC/output/data/scRNASeq/26_sample/UBS93.gene/CID45171_UBS93genelistmean_df.csv")

CID45171_UBS93genelistmean_df        <- read_csv("Pro_TNBC/output/data/scRNASeq/26_sample/UBS93.gene/CID45171_UBS93genelistmean_df.csv")
CID45171_UBS93genelistmean_df        <- CID45171_UBS93genelistmean_df %>% mutate(gene_panel="UBS93")
CID45171_UBS93genelistmean_df$subtype <- rep("HER2+",length(CID45171_UBS93genelistmean_df$cell.id))

patient_UBS93gene_df                    <- bind_rows(CID3963_UBS93genelistmean_df,CID4066_UBS93genelistmean_df,CID4495_UBS93genelistmean_df,CID45171_UBS93genelistmean_df)
patient_UBS93gene_df[patient_UBS93gene_df$celltype_major=="Cancer Epithelial",4]  <- "Cancer cell"
patient_UBS93gene_df[patient_UBS93gene_df$celltype_major=="B-cells",4]  <- "B-cell"
patient_UBS93gene_df[patient_UBS93gene_df$celltype_major=="T-cells",4]  <- "T-cell"
patient_UBS93gene_df[patient_UBS93gene_df$celltype_major=="CAFs",4]  <- "Fibroblast"
patient_UBS93gene_df[patient_UBS93gene_df$celltype_major=="B-cells",4]  <- "B-cell"
patient_UBS93gene_df                    <-  subset(patient_UBS93gene_df,celltype_major!= "Plasmablasts")
patient_UBS93gene_df$celltype_major     <- factor(patient_UBS93gene_df$celltype_major,levels = c("Cancer cell","Fibroblast","PVL","Endothelial","T-cell","Myeloid","B-cell"))
patient_UBS93gene_df$subtype            <- factor(patient_UBS93gene_df$subtype,levels = c("HER2+","ER+/HER2+","ER+","TNBC"))

library(ggpubr)
fig_1b                        <- ggboxplot(patient_UBS93gene_df,x = "celltype_major",y="gene_mean",color="subtype",add = "point",size= 1,add.params=list(size=1))+ 
  ggplot.style + 
  theme(axis.text  = element_text( size=25, face="bold"),plot.title = element_text(size = 55,hjust=0.5))+
  theme(panel.grid.major=element_line(colour=NA),
        panel.grid.minor = element_blank())+xlab("")+ylab("Expression")+
   scale_color_manual(values = c("Red","#1E90FF","Purple","#FFA500"))
ggsave(fig_1b,filename = "Pro_TNBC/paper/plot/section_1/Box.plot.of.the.mean.values.of.gene.panel.in.different.celltypes.of.patients.pdf",height = 15,width = 30)


####fig 1c: t-SNE plot of CCLE using all the genes#####
ccle_info                  <- read.csv("Pro_TNBC/data/CCLE/ccle_info.csv")
brca_meta                  <- CCLE.sample.meta[CCLE.sample.meta$cell.line.tumor.site=="BREAST",]
brca_info                  <- ccle_info[ccle_info$CCLE_Name  %in% brca_meta$cell.line.name,]
brca_info                  <- brca_info[brca_info$primary_disease=="Breast Cancer",]
brca.ccle.log2.rpkm        <- CCLE.log2.rpkm.matrix[,colnames(CCLE.log2.rpkm.matrix) %in% brca_info$CCLE_Name]

library(Rtsne)
cor_matrix                      <- cor(brca.ccle.log2.rpkm,method = "spearman")
dissimilarity_matrix            <- 1 - cor_matrix
normalized_dissimilarity_matrix <- dissimilarity_matrix / max(dissimilarity_matrix)
tsne_result                     <- Rtsne(normalized_dissimilarity_matrix,is_distance = T,perplexity = 15)
plot(tsne_result$Y, col = "blue", pch = 16)
df1                             <- as.data.frame(tsne_result$Y)
df1$CCLE_Name                   <- rownames(normalized_dissimilarity_matrix)
df1                             <- merge(df1,brca_info[,c(4,21)],by="CCLE_Name")
fig_1c                          <- ggplot(df1, aes(x=V1, y=V2, color=lineage_molecular_subtype)) + geom_point(size=6)+
  xlab("t-SNE1")+ 
  ylab("t-SNE2")+ggplot.style+
  theme(legend.text=element_text(size=40,face='bold'),
        legend.title=element_text(size=40,face = 'bold'),
        legend.key.size = unit(1.5, 'lines'),
        legend.position = "bottom")+scale_color_manual(values = c("Red","#1E90FF","Purple","#FFA500","#FF69B4"))+
  theme(panel.grid.major=element_line(colour=NA),
        panel.grid.minor = element_blank())
ggsave(fig_1c,filename = "Pro_TNBC/paper/plot/section_1/TSNE.of.breast.cancer.cell.lines.in.CCLE.using.all.gene.pdf",width=20,height=15)
save(tsne_result,file="Pro_TNBC/paper/data/results/section_1/tsne.of.CLLE.RData")



####Fig S1: Run TSNE for all breast cancer cell lines by UBS93 gene lists using the correlation distance####
ccle_info                  <- read.csv("Pro_TNBC/data/CCLE/ccle_info.csv")
brca_meta                  <- CCLE.sample.meta[CCLE.sample.meta$cell.line.tumor.site=="BREAST",]
brca_info                  <- ccle_info[ccle_info$CCLE_Name  %in% brca_meta$cell.line.name,]
brca_info                  <- brca_info[brca_info$primary_disease=="Breast Cancer",]
brca_info                  <- subset(brca_info,lineage_molecular_subtype!="luminal_HER2_amp")
brca.ccle.log2.rpkm        <- CCLE.log2.rpkm.matrix[,colnames(CCLE.log2.rpkm.matrix) %in% brca_info$CCLE_Name]

brca.ccle.log2.rpkm.genecov             <- as.data.frame(brca.ccle.log2.rpkm)
brca.ccle.log2.rpkm.genecov$ENSEMBL     <- rownames(brca.ccle.log2.rpkm.genecov)
symbol                                  <- bitr(brca.ccle.log2.rpkm.genecov$ENSEMBL, fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Hs.eg.db")
top100_gene_conv                        <- symbol[symbol$SYMBOL %in% top100_gene,]
brca.ccle.log2.rpkm.genecov.100gene     <- merge(top100_gene_conv,brca.ccle.log2.rpkm.genecov,by="ENSEMBL")
rown                                    <- brca.ccle.log2.rpkm.genecov.100gene$SYMBOL
brca.ccle.log2.rpkm.genecov.100gene     <- brca.ccle.log2.rpkm.genecov.100gene[,-c(1,2)]
rownames(brca.ccle.log2.rpkm.genecov.100gene)      <- rown
brca.ccle.log2.rpkm.genecov.100gene                <- as.matrix(brca.ccle.log2.rpkm.genecov.100gene) %>% t()
brca.ccle.log2.rpkm.genecov.100gene                <- t(brca.ccle.log2.rpkm.genecov.100gene)
library(Rtsne)
cor_matrix                      <- cor(brca.ccle.log2.rpkm.genecov.100gene,method = "spearman")
dissimilarity_matrix            <- 1 - cor_matrix
normalized_dissimilarity_matrix <- dissimilarity_matrix / max(dissimilarity_matrix)
tsne_result                     <- Rtsne(normalized_dissimilarity_matrix,is_distance = T,perplexity = 15)
plot(tsne_result$Y, col = "blue", pch = 16)
df1                                               <- as.data.frame(tsne_result$Y)
df1$CCLE_Name                                     <- rownames(normalized_dissimilarity_matrix)
df1                                               <- merge(df1,brca_info[,c(4,21)],by="CCLE_Name")
df1                                               <- within(df1,{
  subtype <- NA
  subtype[lineage_molecular_subtype=="basal_A"]   <- "Basal"
  subtype[lineage_molecular_subtype=="basal_B"]   <- "Claudin_low"
  subtype[lineage_molecular_subtype=="luminal"]   <- "Luminal"
  subtype[lineage_molecular_subtype=="HER2_amp"]  <- "HER2_amp"
})
df1 <- df1[!is.na(df1$subtype),]
supfig_3                                          <- ggplot(df1, aes(V1, y=V2, color=subtype)) + geom_point(size=6)+
  xlab("t-SNE1")+ 
  ylab("t-SNE2")+ggplot.style+
  theme(legend.text=element_text(size=40,face='bold'),
        legend.title=element_text(size=40,face = 'bold'),
        legend.key.size = unit(1.5, 'lines'),
        legend.position = "bottom")+scale_color_manual(values = c("Red","#1E90FF","Purple","#FFA500","#FF69B4"))+
  theme(panel.grid.major=element_line(colour=NA),
        panel.grid.minor = element_blank())
ggsave(supfig_3,filename = "Pro_TNBC/paper/plot/section_1/supfig.3.TSNE.of.breast.cancer.cell.lines.in.CCLE.using.UBS93.gene.panel.pdf",width=20,height=15)
save(tsne_result,file="Pro_TNBC/paper/data/results/section_1/UBS93.gene.tsne.of.CLLE.RData")

