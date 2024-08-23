
library(readr)
library(ggplot2)
library(dplyr)
library(org.Hs.eg.db)
library(clusterProfiler)
library(Seurat)
source('Pro_TNBC/paper/code/ggplot.style.R')

load("Pro_TNBC/data/CCLE/CCLE.RData")
load("Pro_TNBC/output/data/CCLE/UBS93.Rpackage/UBS93.gene.df.RData")
load("Pro_TNBC/output/data/scRNASeq/26_sample/GSE176078_scRNA.RData")
load("Pro_TNBC/output/data/scRNASeq/26_sample/top100.gene/top100_gene.RData")
######################################################################################################################
#  Fig 1:Development of UBS93 gene panel.
######################################################################################################################

######################################################################################################################
#  Fig S1: Scatter plot illustrating the relationship between the average E-M ratio and the rank of 27,719 genes (a), the top 2,720 genes (b), and the top 500 genes (c)
######################################################################################################################
source('Pro_TNBC/paper/code/ggplot.style.R')
library(ggplot2)
library(quantreg)
library(readr)
mean_score                           <- read_csv("Pro_TNBC/paper/data/results/section_1/mean_score.csv")
mean_score                           <- mean_score[order(mean_score$S_mean,decreasing = F),]
mean_score$rank                      <- (1: length(mean_score$S_mean))
figS1a                               <- ggplot(mean_score,aes(x=rank,y=S_mean))+
  geom_point(size=6,show.legend = F)+ggplot.style+xlab(paste("rank",sep = "\n","(n=27719)"))+ylab("score")+
  scale_x_continuous(limits = c(0, 28000),breaks = seq(0,28000,5000))+
  geom_quantile(data = mean_score, mapping = aes(
    x = rank, y = S_mean,color="red",linewidth = 1),quantiles = 0.5,show.legend = F)+
  theme(panel.grid.major=element_line(colour=NA),
        panel.grid.minor = element_blank())+geom_vline(aes(xintercept=25000,color="blue",linewidth=1))
ggsave(figS1a,filename = "Pro_TNBC/paper/plot/method/the.rank.charts.of.S_mean.pdf",width=20,height=15)
mean_score_2720                       <- mean_score[25000:27719,]
mean_score_2720$rank                  <- (1: length(mean_score_2720$S_mean))
figS1b                               <- ggplot(mean_score_2720,aes(x=rank,y=S_mean))+
  geom_point(size=6,show.legend = F)+ggplot.style+
  xlab(paste("rank",sep = "\n","(n=2720)"))+ylab("score")+
  scale_x_continuous(limits = c(0, 2800),breaks = seq(0,2800,500))+ 
  geom_quantile(data = mean_score_2720, 
                mapping = aes(x = rank, y = S_mean,color="red",linewidth = 1),
                quantiles = 0.5,show.legend = F)+
  theme(panel.grid.major=element_line(colour=NA),
        panel.grid.minor = element_blank())+geom_vline(aes(xintercept=2220,color="blue",linewidth=1))
ggsave(figS1b,filename = "Pro_TNBC/paper/plot/method/the.second.rank.charts.of.S_mean.pdf",width=20,height=15)
mean_score_500                       <- mean_score_2720[2221:2720,]
mean_score_500$rank                  <- (1: length(mean_score_500$S_mean))
figS1c                              <- ggplot(mean_score_500,aes(x=rank,y=S_mean))+
  geom_point(size=6,show.legend = F)+ggplot.style+
  xlab(paste("rank",sep = "\n","(n=500)"))+ylab("score")+
  scale_x_continuous(limits = c(0, 500),breaks = seq(0,500,100))+ 
  geom_quantile(data = mean_score_500, 
                mapping = aes(x = rank, y = S_mean,color="red",linewidth = 1),
                quantiles = 0.5,show.legend = F)+
  theme(panel.grid.major=element_line(colour=NA),
        panel.grid.minor = element_blank())+geom_vline(aes(xintercept=400,color="blue",linewidth = 1))
ggsave(figS1c,filename = "Pro_TNBC/paper/plot/method/the.third.rank.charts.of.S_mean.pdf",width=20,height=15)
top100_gene           <- mean_score_500[401:500,]$gene
ensembl               <- bitr(top100_gene,fromType = "SYMBOL",toType = "ENSEMBL",OrgDb = "org.Hs.eg.db")
top100_gene_ensembl   <- intersect(ensembl$ENSEMBL,rownames(CCLE.log2.rpkm.matrix))
ENTREZID              <- bitr(top100_gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
A                     <- merge(ensembl,ENTREZID,by="SYMBOL")
UBS93.gene.df         <- A[A$ENSEMBL %in% top100_gene_ensembl,]
save(UBS93.gene.df,file="Pro_TNBC/paper/data/results/section_1/UBS93.gene.df.RData")

#### fig 1b: Box plot of the mean values of UBS93 gene in different cell types of CID3963.####
#CID3963(ER+)#
library(ggplot2)
library(gridExtra)
library(cowplot)
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
CID3963_UBS93genelistmean_df          <- subset(CID3963_UBS93genelistmean_df,!celltype_major=="Normal Epithelial")
CID3963_UBS93genelistmean_df$celltype_major    <- factor(CID3963_UBS93genelistmean_df$celltype_major,levels = c("Cancer Epithelial","Endothelial","CAFs","PVL","Myeloid","T-cells"))
fig_1b                               <- ggboxplot(CID3963_UBS93genelistmean_df,x = "celltype_major",y="gene_mean",add = "point",size= 1,add.params=list(size=1))+ 
  ggplot.style+ 
  theme(axis.text  = element_text( size=35, face="bold"),plot.title = element_text(size = 55,hjust=0.5))+
  theme(panel.grid.major=element_line(colour=NA),
        panel.grid.minor = element_blank())+xlab("Celltype")+ylab("Expression")+
  scale_x_discrete("Celltype", labels = c("Cancer Epithelial" = paste("Cancer cell",sep = "\n","(n=195)"),"CAFs" = paste("CAF",sep = "\n","(n=23)"),"Myeloid" = paste("Myeloid",sep = "\n","(n=474)"),"PVL" = paste("PVL",sep = "\n","(n=28)"),"T-cells" = paste("T-cells",sep = "\n","(n=2672)"),"Endothelial" = paste("Endothelial",sep = "\n","(n=102)")))
ggsave(fig_1b,filename = "Pro_TNBC/paper/plot/section_1/Box.plot.of.the.mean.values.of.gene.panel.in.different.celltypes.of.CID3963.pdf",height = 15,width = 20)


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
        legend.position = "bottom")+scale_color_manual(values = c("Red","Blue","Purple","#FFCC00","#E95790"))+
  theme(panel.grid.major=element_line(colour=NA),
        panel.grid.minor = element_blank())
ggsave(fig_1c,filename = "Pro_TNBC/paper/plot/section_1/TSNE.of.breast.cancer.cell.lines.in.CCLE.using.all.gene.pdf",width=20,height=15)
save(tsne_result,file="Pro_TNBC/paper/data/results/section_1/tsne.of.CLLE.RData")



######################################################################################################################
#  Fig S2: Box plot of the mean values of TOP100 gene expression in different cell types
######################################################################################################################

####*Figure S2a: CID4066(ER+/HER2+)####
library(ggplot2)
library(gridExtra)
library(cowplot)
CID4066_UBS93genelist_exprmat         <- CID4066_exprmat_CPM[rownames(CID4066_exprmat_CPM) %in% UBS93.gene.df$SYMBOL,]
gene_mean                             <- apply(CID4066_UBS93genelist_exprmat,2,mean)
CID4066_UBS93genelistmean_df          <- as.data.frame(gene_mean)
CID4066_UBS93genelistmean_df$cell.id  <- rownames(CID4066_UBS93genelistmean_df)
library(readr)
metadata                <- read_csv("Pro_TNBC/data/scRNASeq/26/Wu_etal_2021_BRCA_scRNASeq/metadata.csv")
CID4066_metadata        <- metadata[metadata$orig.ident=="CID4066",]
CID4066                 <- CID4066_metadata[,c(1,9)]
colnames(CID4066)[1]    <- "cell.id"
CID4066_UBS93genelistmean_df <- merge(CID4066_UBS93genelistmean_df,CID4066,by="cell.id")
write.csv(CID4066_UBS93genelistmean_df,file = "Pro_TNBC/output/data/scRNASeq/26_sample/UBS93.gene/CID4066_UBS93genelistmean_df.csv")

CID4066_UBS93genelistmean_df   <- read_csv("Pro_TNBC/output/data/scRNASeq/26_sample/UBS93.gene/CID4066_UBS93genelistmean_df.csv")
CID4066_UBS93genelistmean_df   <- CID4066_UBS93genelistmean_df %>% mutate(gene_panel="UBS93")
CID4066_UBS93genelistmean_df     <- subset(CID4066_UBS93genelistmean_df,celltype_major!="Normal Epithelial")
CID4066_UBS93genelistmean_df$celltype_major     <- factor(CID4066_UBS93genelistmean_df$celltype_major,levels = c("Cancer Epithelial","CAFs","Endothelial","PVL","B-cells","Myeloid","T-cells"))
supfig_2a                      <- ggboxplot(CID4066_UBS93genelistmean_df,x = "celltype_major",y="gene_mean",add = "point",size= 1,add.params=list(size=1))+ 
  ggplot.style + 
  theme(axis.text  = element_text( size=25, face="bold"),plot.title = element_text(size = 55,hjust=0.5))+
  theme(panel.grid.major=element_line(colour=NA),
        panel.grid.minor = element_blank())+xlab("Celltype")+ylab("Expression")+
  scale_x_discrete("Celltype", labels = c("Cancer Epithelial" = paste("Cancer cell",sep = "\n","(n=514)"),"CAFs" = paste("CAF",sep = "\n","(n=923)"),"Myeloid" = paste("Myeloid",sep = "\n","(n=221)"),"PVL" = paste("PVL",sep = "\n","(n=630)"),"T-cells" = paste("T-cells",sep = "\n","(n=2171)"),"Endothelial" = paste("Endothelial",sep = "\n","(n=535)"),"B-cells"=paste("B-cells",sep = "\n","(n=38)")))
ggsave(supfig_2a,filename = "Pro_TNBC/paper/plot/section_1/Box.plot.of.the.mean.values.of.gene.panel.in.different.celltypes.of.CID4066.pdf",height = 15,width = 20)

####*Figure S2b: CID4495(TNBC)####
library(ggplot2)
library(gridExtra)
library(cowplot)
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
supfig_2b                        <- ggplot(CID4495_UBS93genelistmean_df,aes(x=reorder(celltype_major,gene_mean,decreasing=T),y=gene_mean))+
  geom_boxplot(width=0.6,alpha=0.8)+ggplot.style + 
  theme(axis.text  = element_text( size=25, face="bold"),plot.title = element_text(size = 55,hjust=0.5))+
  geom_point(size=1,position="jitter")+
  xlab("Celltype")+ylab("Expression")+
  theme(panel.grid.major=element_line(colour=NA),
        panel.grid.minor = element_blank())+
  scale_x_discrete("Celltype", labels = c("Cancer Epithelial" = paste("Cancer cell",sep = "\n","(n=1177)"),"CAFs" = paste("CAF",sep = "\n","(n=232)"),"Myeloid" = paste("Myeloid",sep = "\n","(n=897)"),"PVL" = paste("PVL",sep = "\n","(n=191)"),"T-cells" = paste("T-cells",sep = "\n","(n=3504)"),"Endothelial" = paste("Endothelial",sep = "\n","(n=184)"),"B-cells"=paste("B-cells",sep = "\n","(n=773)"),"Plasmablasts"=paste("Plasmablasts",sep = "\n","(n=1020)")))
CID4495_UBS93genelistmean_df$celltype_major      <- factor(CID4495_UBS93genelistmean_df$celltype_major,levels = c("Cancer Epithelial","CAFs","Endothelial","PVL","Plasmablasts","Myeloid","B-cells","T-cells"))
supfig_2b                      <- ggboxplot(CID4495_UBS93genelistmean_df,x = "celltype_major",y="gene_mean",add = "point",size= 1,add.params=list(size=1))+ 
  ggplot.style+ 
  theme(axis.text  = element_text( size=25, face="bold"),plot.title = element_text(size = 55,hjust=0.5))+
  scale_color_manual(values =c(PAM50="blue",UBS93="red"))+
  theme(panel.grid.major=element_line(colour=NA),
        panel.grid.minor = element_blank())+xlab("Celltype")+ylab("Expression")+
  scale_x_discrete("Celltype", labels = c("Cancer Epithelial" = paste("Cancer cell",sep = "\n","(n=1177)"),"CAFs" = paste("CAF",sep = "\n","(n=232)"),"Myeloid" = paste("Myeloid",sep = "\n","(n=897)"),"PVL" = paste("PVL",sep = "\n","(n=191)"),"T-cells" = paste("T-cells",sep = "\n","(n=3504)"),"Endothelial" = paste("Endothelial",sep = "\n","(n=184)"),"B-cells"=paste("B-cells",sep = "\n","(n=773)"),"Plasmablasts"=paste("Plasmablasts",sep = "\n","(n=1020)")))
ggsave(supfig_2b,filename = "Pro_TNBC/paper/plot/section_1/Box.plot.of.the.mean.values.of.gene.panel.in.different.celltypes.of.CID4495.pdf",height = 15,width = 20)

####*Figure S2c: CID45171(HER2+)####
library(ggplot2)
library(gridExtra)
library(cowplot)
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

CID45171_UBS93genelistmean_df   <- read_csv("Pro_TNBC/output/data/scRNASeq/26_sample/UBS93.gene/CID45171_UBS93genelistmean_df.csv")
CID45171_UBS93genelistmean_df   <- CID45171_UBS93genelistmean_df %>% mutate(gene_panel="UBS93")
supfig_2c                        <- ggplot(CID45171_UBS93genelistmean_df,aes(x=reorder(celltype_major,gene_mean,decreasing=T),y=gene_mean))+
  geom_boxplot(width=0.6,alpha=0.8)+ggplot.style + 
  theme(axis.text  = element_text( size=25, face="bold"),plot.title = element_text(size = 55,hjust=0.5))+
  geom_point(size=1,position = "jitter")+
  xlab("Celltype")+ylab("Expression")+
  theme(panel.grid.major=element_line(colour=NA),
        panel.grid.minor = element_blank())+
  scale_x_discrete("Celltype", labels = c("Cancer Epithelial" = paste("Cancer cell",sep = "\n","(n=805)"),"CAFs" = paste("CAF",sep = "\n","(n=32)"),"Myeloid" = paste("Myeloid",sep = "\n","(n=172)"),"PVL" = paste("PVL",sep = "\n","(n=13)"),"T-cells" = paste("T-cells",sep = "\n","(n=1346)"),"Endothelial" = paste("Endothelial",sep = "\n","(n=15)"),"B-cells"=paste("B-cells",sep = "\n","(n=56)")))
CID45171_UBS93genelistmean_df$celltype_major     <- factor(CID45171_UBS93genelistmean_df$celltype_major,levels = c("Cancer Epithelial","CAFs","PVL","Endothelial","T-cells","Myeloid","B-cells"))
supfig_2c                       <- ggboxplot(CID45171_UBS93genelistmean_df,x = "celltype_major",y="gene_mean",add = "point",size= 1,add.params=list(size=1))+ 
  ggplot.style + 
  theme(axis.text  = element_text( size=25, face="bold"),plot.title = element_text(size = 55,hjust=0.5))+
  theme(panel.grid.major=element_line(colour=NA),
        panel.grid.minor = element_blank())+xlab("Celltype")+ylab("Expression")+
  scale_x_discrete("Celltype", labels = c("Cancer Epithelial" = paste("Cancer cell",sep = "\n","(n=805)"),"CAFs" = paste("CAF",sep = "\n","(n=32)"),"Myeloid" = paste("Myeloid",sep = "\n","(n=172)"),"PVL" = paste("PVL",sep = "\n","(n=13)"),"T-cells" = paste("T-cells",sep = "\n","(n=1346)"),"Endothelial" = paste("Endothelial",sep = "\n","(n=15)"),"B-cells"=paste("B-cells",sep = "\n","(n=56)")))
ggsave(supfig_2c,filename = "Pro_TNBC/paper/plot/section_1/Box.plot.of.the.mean.values.of.gene.panel.in.different.celltypes.of.CID45171.pdf",height = 15,width = 20)




######################################################################################################################
#  Fig S3: Run TSNE for all breast cancer cell lines by UBS93 gene lists using the correlation distance
######################################################################################################################

####Figure S3#####
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
supfig_3                                          <- ggplot(df1, aes(x=V1, y=V2, color=subtype)) + geom_point(size=6)+
  xlab("t-SNE1")+ 
  ylab("t-SNE2")+ggplot.style+
  theme(legend.text=element_text(size=40,face='bold'),
        legend.title=element_text(size=40,face = 'bold'),
        legend.key.size = unit(1.5, 'lines'),
        legend.position = "bottom")+scale_color_manual(values = c("Red","Blue","Purple","#FFCC00","Black"))+
  theme(panel.grid.major=element_line(colour=NA),
        panel.grid.minor = element_blank())
ggsave(supfig_3,filename = "Pro_TNBC/paper/plot/section_1/supfig.3.TSNE.of.breast.cancer.cell.lines.in.CCLE.using.UBS93.gene.panel.pdf",width=20,height=15)
save(tsne_result,file="Pro_TNBC/paper/data/results/section_1/UBS93.gene.tsne.of.CLLE.RData")
