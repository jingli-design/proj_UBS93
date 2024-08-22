################################################################
#### 3.1: score and mean expression in different cell types
################################################################


#### fig 3a: Box plot of the mean score of UBS93 gene panel and PAM50 gene panel.####
library(readr)
library(ggplot2)
load("Pro_TNBC/paper/data/results/section_1/UBS93.gene.df.RData")
load("~/Pro_TNBC/output/data/scRNASeq/26_sample/pam50.gene/PAM50.gene.name.RData")
mean_score         <- read_csv("Pro_TNBC/paper/data/results/section_1/mean_score.csv")
UBS93_score        <- mean_score[mean_score$gene %in% UBS93.gene.df$SYMBOL,]
PAM50_score        <- mean_score[mean_score$gene %in% pam50.gene,]
UBS93_pam50_score  <- data.frame(s_mean=c(UBS93_score$S_mean,PAM50_score$S_mean),group=c(rep("UBS93_score",93),rep("pam50_gene",50)))
fig_3a             <- ggplot(UBS93_pam50_score,aes(x=group,y=s_mean))+
  geom_boxplot(outlier.shape=NA,lwd=1.2)+
  geom_point(size=4,position="jitter") + ggplot.style + 
  xlab('') + ylab('score')+
  theme(panel.grid.major=element_line(colour=NA),
        panel.grid.minor = element_blank())+
  scale_x_discrete("", labels = c("UBS93_score" = paste("UBS93",sep = "\n","(n=93)"),"pam50_gene" = paste("PAM50",sep = "\n","(n=50)")))
ggsave(fig_3a,filename = "Pro_TNBC/paper/plot/section_3/the.score.of.UBS93.gene.panel.and.pam50.gene.panel.pdf",width=20,height=15)
wilcox.test(s_mean~group,data = UBS93_pam50_score)#p-value < 2.2e-16
save(UBS93_pam50_score,file="Pro_TNBC/paper/data/section_3/fig3a.UBS93.pam50.score.RData")

#### fig 3b: Box plot of the mean values of UBS93 gene and PAM50 gene expression in different cell types of CID3963.####
#CID3963(ER+)#
library(ggplot2)
library(gridExtra)
library(cowplot)
CID3963_UBS93genelist_exprmat         <- CID3963_exprmat_CPM[rownames(CID3963_exprmat_CPM) %in% BCT93.gene.df$SYMBOL,]
gene_mean                             <- apply(CID3963_UBS93genelist_exprmat,2,mean)
CID3963_UBS93genelistmean_df          <- as.data.frame(gene_mean)
CID3963_UBS93genelistmean_df$cell.id  <- rownames(CID3963_UBS93genelistmean_df)
library(readr)
metadata                              <- read_csv("Pro_TNBC/data/scRNASeq/26/Wu_etal_2021_BRCA_scRNASeq/metadata.csv")
CID3963_metadata                      <- metadata[metadata$orig.ident=="CID3963",]
CID3963                               <- CID3963_metadata[,c(1,9)]
colnames(CID3963)[1]                  <- "cell.id"
CID3963_UBS93genelistmean_df          <- merge(CID3963_UBS93genelistmean_df,CID3963,by="cell.id")
library(ggplot2)
ggplot(CID3963_UBS93genelistmean_df,aes(x=celltype_major,y=gene_mean))+geom_boxplot()
write.csv(CID3963_UBS93genelistmean_df,file = "Pro_TNBC/output/data/scRNASeq/26_sample/UBS93.gene/CID3963_UBS93genelistmean_df.csv")

CID3963_UBS93genelistmean_df   <- read_csv("Pro_TNBC/output/data/scRNASeq/26_sample/UBS93.gene/CID3963_UBS93genelistmean_df.csv")
CID3963_UBS93genelistmean_df   <- CID3963_UBS93genelistmean_df %>% mutate(gene_panel="UBS93")
CID3963_pam50genelistmean_df   <- read_csv(file = "Pro_TNBC/output/data/scRNASeq/26_sample/pam50.gene/CID3963_pam50genelistmean_df.csv")
CID3963_pam50genelistmean_df   <- CID3963_pam50genelistmean_df %>% mutate(gene_panel="PAM50")
CID3963_genelistmean_df        <- rbind(CID3963_UBS93genelistmean_df,CID3963_pam50genelistmean_df)
CID3963_genelistmean_df        <- subset(CID3963_genelistmean_df,!celltype_major=="Normal Epithelial")
fig_3b                         <- ggplot(CID3963_genelistmean_df,aes(x=reorder(celltype_major,gene_mean,decreasing =T),y=gene_mean,fill=gene_panel))+
  geom_boxplot(width=0.6,alpha=0.8,position=position_dodge(0.8))+ggplot.style + 
  theme(axis.text  = element_text( size=27, face="bold"),plot.title = element_text(size = 55,hjust=0.5))+
  scale_fill_manual(values =c(PAM50="#FF66FF",UBS93="#FF9900"))+
  geom_point(size=1,position = position_jitterdodge(dodge.width=0.8))+
  theme(panel.grid.major=element_line(colour=NA),
        panel.grid.minor = element_blank())+xlab("Celltype")+ylab("Expression")+
  scale_x_discrete("Celltype", labels = c("Cancer Epithelial" = paste("Cancer cell",sep = "\n","(n=195)"),"CAFs" = paste("CAF",sep = "\n","(n=23)"),"Myeloid" = paste("Myeloid",sep = "\n","(n=474)"),"PVL" = paste("PVL",sep = "\n","(n=28)"),"T-cells" = paste("T-cells",sep = "\n","(n=2672)"),"Endothelial" = paste("Endothelial",sep = "\n","(n=102)")))
ggsave(fig_3b,filename = "Pro_TNBC/paper/plot/section_3/Box.plot.of.the.mean.values.of.gene.panel.in.different.celltypes.of.CID3963.pdf",height = 15,width = 25)
save(CID3963_genelistmean_df,file = "Pro_TNBC/paper/data/section_3/fig1b2.CID3963.genelistmean.RData")

#Wilcoxon Signed Rank Test
cell_groups               <- names(table(CID3963_genelistmean_df$celltype_major))
results                   <- data.frame(CellType = character(), p.value = numeric(), stringsAsFactors = FALSE)
for (i in 1:length(cell_groups)) {
  CID3963_genelistmean    <- subset(CID3963_genelistmean_df,celltype_major==cell_groups[i])
  p_value                 <- wilcox.test(CID3963_genelistmean$gene_mean ~ CID3963_genelistmean$gene_panel, paired = TRUE)
  result_row              <- data.frame(CellType = cell_groups[i], p.value = p_value$p.value, stringsAsFactors = FALSE)
  results                 <- rbind(results, result_row)
}

CID3963_PAM50  <- subset(CID3963_genelistmean_df,gene_panel=="PAM50_gene")
CID3963_PAM50  <- subset(CID3963_PAM50,celltype_major=="CAFs"|celltype_major=="Cancer Epithelial")
wilcox.test(gene_mean~celltype_major,data = CID3963_PAM50)
#p = 0.1274

#### supplymentary figure 8:Box plot of the mean values of UBS93 gene and PAM50 gene expression in different cell types.####
#### Figure S8a: CID4495(TNBC)#
library(ggplot2)
library(gridExtra)
library(cowplot)
CID4495_UBS93genelist_exprmat         <- CID4495_exprmat_CPM[rownames(CID4495_exprmat_CPM) %in% BCT93.gene.df$SYMBOL,]
gene_mean                             <- apply(CID4495_UBS93genelist_exprmat,2,mean)
CID4495_UBS93genelistmean_df          <- as.data.frame(gene_mean)
CID4495_UBS93genelistmean_df$cell.id  <- rownames(CID4495_UBS93genelistmean_df)
library(readr)
metadata                <- read_csv("Pro_TNBC/data/scRNASeq/26/Wu_etal_2021_BRCA_scRNASeq/metadata.csv")
CID4495_metadata        <- metadata[metadata$orig.ident=="CID4495",]
CID4495                 <- CID4495_metadata[,c(1,9)]
colnames(CID4495)[1]    <- "cell.id"
CID4495_UBS93genelistmean_df <- merge(CID4495_UBS93genelistmean_df,CID4495,by="cell.id")
library(ggplot2)
ggplot(CID4495_UBS93genelistmean_df,aes(x=celltype_major,y=gene_mean))+geom_boxplot()
write.csv(CID4495_UBS93genelistmean_df,file = "Pro_TNBC/output/data/scRNASeq/26_sample/UBS93.gene/CID4495_UBS93genelistmean_df.csv")

CID4495_UBS93genelistmean_df   <- read_csv("Pro_TNBC/output/data/scRNASeq/26_sample/UBS93.gene/CID4495_UBS93genelistmean_df.csv")
CID4495_UBS93genelistmean_df   <- CID4495_UBS93genelistmean_df %>% mutate(gene_panel="UBS93")
CID4495_pam50genelistmean_df   <- read_csv(file = "Pro_TNBC/output/data/scRNASeq/26_sample/pam50.gene/CID4495_pam50genelistmean_df.csv")
CID4495_pam50genelistmean_df   <- CID4495_pam50genelistmean_df %>% mutate(gene_panel="PAM50")
CID4495_genelistmean_df        <- rbind(CID4495_UBS93genelistmean_df,CID4495_pam50genelistmean_df)
CID4495_genelistmean_df        <- subset(CID4495_genelistmean_df,!celltype_major=="Normal Epithelial")
supfig_8a                        <- ggplot(CID4495_genelistmean_df,aes(x=reorder(celltype_major,gene_mean,decreasing =T),y=gene_mean,fill=gene_panel))+
  geom_boxplot(width=0.6,alpha=0.8,position=position_dodge(0.8))+ggplot.style + 
  theme(axis.text  = element_text( size=27, face="bold"),plot.title = element_text(size = 55,hjust=0.5))+
  scale_fill_manual(values =c(PAM50="#FF66FF",UBS93="#FF9900"))+
  geom_point(size=0.5,position = position_jitterdodge(dodge.width=0.8))+
  theme(panel.grid.major=element_line(colour=NA),
        panel.grid.minor = element_blank())+xlab("Celltype")+ylab("Expression")+
  scale_x_discrete("Celltype", labels = c("Cancer Epithelial" = paste("Cancer cell",sep = "\n","(n=1177)"),"CAFs" = paste("CAF",sep = "\n","(n=232)"),"Myeloid" = paste("Myeloid",sep = "\n","(n=897)"),"PVL" = paste("PVL",sep = "\n","(n=191)"),"T-cells" = paste("T-cells",sep = "\n","(n=3504)"),"Endothelial" = paste("Endothelial",sep = "\n","(n=184)"),"B-cells"=paste("B-cells",sep = "\n","(n=773)"),"Plasmablasts"=paste("Plasmablasts",sep = "\n","(n=1020)")))
ggsave(supfig_8a,filename = "Pro_TNBC/paper/plot/section_3/Box.plot.of.the.mean.values.of.gene.panel.in.different.celltypes.of.CID4495.pdf",height = 15,width = 25)
save(CID4495_genelistmean_df,file = "Pro_TNBC/paper/data/section_3/fig1b2.CID4495.genelistmean.RData")

#Wilcoxon Signed Rank Test
cell_groups               <- names(table(CID4495_genelistmean_df$celltype_major))
results                   <- data.frame(CellType = character(), p.value = numeric(), stringsAsFactors = FALSE)
for (i in 1:length(cell_groups)) {
  CID4495_genelistmean    <- subset(CID4495_genelistmean_df,celltype_major==cell_groups[i])
  p_value                 <- wilcox.test(CID4495_genelistmean$gene_mean ~ CID4495_genelistmean$gene_panel, paired = TRUE)
  result_row              <- data.frame(CellType = cell_groups[i], p.value = p_value$p.value, stringsAsFactors = FALSE)
  results                 <- rbind(results, result_row)
}

CID4495_PAM50 <- subset(CID4495_genelistmean_df,gene_panel=="PAM50_gene")
CID4495_PAM50 <- subset(CID4495_PAM50,celltype_major=="CAFs"|celltype_major=="Cancer Epithelial")
wilcox.test(gene_mean~celltype_major,data = CID4495_PAM50)
#p-value = 6.536e-10

####Figure S8b: CID4066(ER+/HER2+)#
library(ggplot2)
library(gridExtra)
library(cowplot)
CID4066_UBS93genelist_exprmat         <- CID4066_exprmat_CPM[rownames(CID4066_exprmat_CPM) %in% BCT93.gene.df$SYMBOL,]
gene_mean                             <- apply(CID4066_UBS93genelist_exprmat,2,mean)
CID4066_UBS93genelistmean_df          <- as.data.frame(gene_mean)
CID4066_UBS93genelistmean_df$cell.id  <- rownames(CID4066_UBS93genelistmean_df)
library(readr)
metadata                <- read_csv("Pro_TNBC/data/scRNASeq/26/Wu_etal_2021_BRCA_scRNASeq/metadata.csv")
CID4066_metadata        <- metadata[metadata$orig.ident=="CID4066",]
CID4066                 <- CID4066_metadata[,c(1,9)]
colnames(CID4066)[1]    <- "cell.id"
CID4066_UBS93genelistmean_df <- merge(CID4066_UBS93genelistmean_df,CID4066,by="cell.id")
library(ggplot2)
ggplot(CID4066_UBS93genelistmean_df,aes(x=celltype_major,y=gene_mean))+geom_boxplot()
write.csv(CID4066_UBS93genelistmean_df,file = "Pro_TNBC/output/data/scRNASeq/26_sample/UBS93.gene/CID4066_UBS93genelistmean_df.csv")

CID4066_UBS93genelistmean_df   <- read_csv("Pro_TNBC/output/data/scRNASeq/26_sample/UBS93.gene/CID4066_UBS93genelistmean_df.csv")
CID4066_UBS93genelistmean_df   <- CID4066_UBS93genelistmean_df %>% mutate(gene_panel="UBS93")
CID4066_pam50genelistmean_df   <- read_csv(file = "Pro_TNBC/output/data/scRNASeq/26_sample/pam50.gene/CID4066_pam50genelistmean_df.csv")
CID4066_pam50genelistmean_df   <- CID4066_pam50genelistmean_df %>% mutate(gene_panel="PAM50")
CID4066_genelistmean_df        <- rbind(CID4066_UBS93genelistmean_df,CID4066_pam50genelistmean_df)
CID4066_genelistmean_df        <- subset(CID4066_genelistmean_df,!celltype_major=="Normal Epithelial")
supfig_8b                      <- ggplot(CID4066_genelistmean_df,aes(x=reorder(celltype_major,gene_mean,decreasing =T),y=gene_mean,fill=gene_panel))+
  geom_boxplot(width=0.6,alpha=0.8,position=position_dodge(0.8))+ggplot.style + 
  theme(axis.text  = element_text( size=27, face="bold"),plot.title = element_text(size = 55,hjust=0.5))+
  scale_fill_manual(values =c(PAM50="#FF66FF",UBS93="#FF9900"))+
  geom_point(size=0.5,position = position_jitterdodge(dodge.width=0.8))+
  theme(panel.grid.major=element_line(colour=NA),
        panel.grid.minor = element_blank())+xlab("Celltype")+ylab("Expression")+
  scale_x_discrete("Celltype", labels = c("Cancer Epithelial" = paste("Cancer cell",sep = "\n","(n=514)"),"CAFs" = paste("CAF",sep = "\n","(n=923)"),"Myeloid" = paste("Myeloid",sep = "\n","(n=221)"),"PVL" = paste("PVL",sep = "\n","(n=630)"),"T-cells" = paste("T-cells",sep = "\n","(n=2171)"),"Endothelial" = paste("Endothelial",sep = "\n","(n=535)"),"B-cells"=paste("B-cells",sep = "\n","(n=38)")))
ggsave(supfig_8b,filename = "Pro_TNBC/paper/plot/section_3/Box.plot.of.the.mean.values.of.gene.panel.in.different.celltypes.of.CID4066.pdf",height = 15,width = 25)
save(CID4066_genelistmean_df,file = "Pro_TNBC/paper/data/section_3/fig1b2.CID4066.genelistmean.RData")

#Wilcoxon Signed Rank Test
cell_groups               <- names(table(CID4066_genelistmean_df$celltype_major))
results                   <- data.frame(CellType = character(), p.value = numeric(), stringsAsFactors = FALSE)
for (i in 1:length(cell_groups)) {
  CID4066_genelistmean    <- subset(CID4066_genelistmean_df,celltype_major==cell_groups[i])
  p_value                 <- wilcox.test(CID4066_genelistmean$gene_mean ~ CID4066_genelistmean$gene_panel, paired = TRUE)
  result_row              <- data.frame(CellType = cell_groups[i], p.value = p_value$p.value, stringsAsFactors = FALSE)
  results                 <- rbind(results, result_row)}

CID4066_PAM50 <- subset(CID4066_genelistmean_df,gene_panel=="PAM50_gene")
CID4066_PAM50 <- subset(CID4066_PAM50,celltype_major=="CAFs"|celltype_major=="Cancer Epithelial")
wilcox.test(gene_mean~celltype_major,data = CID4066_PAM50)
#p-value < 2.2e-16


####Figure S8c:CID45171(HER2+)#
library(ggplot2)
library(gridExtra)
library(cowplot)
CID45171_UBS93genelist_exprmat         <- CID45171_exprmat_CPM[rownames(CID45171_exprmat_CPM) %in% BCT93.gene.df$SYMBOL,]
gene_mean                              <- apply(CID45171_UBS93genelist_exprmat,2,mean)
CID45171_UBS93genelistmean_df          <- as.data.frame(gene_mean)
CID45171_UBS93genelistmean_df$cell.id  <- rownames(CID45171_UBS93genelistmean_df)
library(readr)
metadata                 <- read_csv("Pro_TNBC/data/scRNASeq/26/Wu_etal_2021_BRCA_scRNASeq/metadata.csv")
CID45171_metadata        <- metadata[metadata$orig.ident=="CID45171",]
CID45171                 <- CID45171_metadata[,c(1,9)]
colnames(CID45171)[1]    <- "cell.id"
CID45171_UBS93genelistmean_df           <- merge(CID45171_UBS93genelistmean_df,CID45171,by="cell.id")
library(ggplot2)
ggplot(CID45171_UBS93genelistmean_df,aes(x=celltype_major,y=gene_mean))+geom_boxplot()
write.csv(CID45171_UBS93genelistmean_df,file = "Pro_TNBC/output/data/scRNASeq/26_sample/UBS93.gene/CID45171_UBS93genelistmean_df.csv")

CID45171_UBS93genelistmean_df   <- read_csv("Pro_TNBC/output/data/scRNASeq/26_sample/UBS93.gene/CID45171_UBS93genelistmean_df.csv")
CID45171_UBS93genelistmean_df   <- CID45171_UBS93genelistmean_df %>% mutate(gene_panel="UBS93")
CID45171_pam50genelistmean_df   <- read_csv(file = "Pro_TNBC/output/data/scRNASeq/26_sample/pam50.gene/CID45171_pam50genelistmean_df.csv")
CID45171_pam50genelistmean_df   <- CID45171_pam50genelistmean_df %>% mutate(gene_panel="PAM50")
CID45171_genelistmean_df        <- rbind(CID45171_UBS93genelistmean_df,CID45171_pam50genelistmean_df)
CID45171_genelistmean_df        <- subset(CID45171_genelistmean_df,!celltype_major=="Normal Epithelial")
supfig_8c                       <- ggplot(CID45171_genelistmean_df,aes(x=reorder(celltype_major,gene_mean,decreasing =T),y=gene_mean,fill=gene_panel))+
  geom_boxplot(width=0.6,alpha=0.8,position=position_dodge(0.8))+ggplot.style + 
  theme(axis.text  = element_text( size=27, face="bold"),plot.title = element_text(size = 55,hjust=0.5))+
  scale_fill_manual(values =c(PAM50="#FF66FF",UBS93="#FF9900"))+
  geom_point(size=0.5,position = position_jitterdodge(dodge.width=0.8))+
  theme(panel.grid.major=element_line(colour=NA),
        panel.grid.minor = element_blank())+xlab("Celltype")+ylab("Expression")+
  scale_x_discrete("Celltype", labels = c("Cancer Epithelial" = paste("Cancer cell",sep = "\n","(n=805)"),"CAFs" = paste("CAF",sep = "\n","(n=32)"),"Myeloid" = paste("Myeloid",sep = "\n","(n=172)"),"PVL" = paste("PVL",sep = "\n","(n=13)"),"T-cells" = paste("T-cells",sep = "\n","(n=1346)"),"Endothelial" = paste("Endothelial",sep = "\n","(n=15)"),"B-cells"=paste("B-cells",sep = "\n","(n=56)")))
ggsave(supfig_8c,filename = "Pro_TNBC/paper/plot/section_3/Box.plot.of.the.mean.values.of.gene.panel.in.different.celltypes.of.CID45171.pdf",height = 15,width = 25)
save(CID45171_genelistmean_df,file = "Pro_TNBC/paper/data/section_3/fig1b2.CID45171.genelistmean.RData")

#Wilcoxon Signed Rank Test
cell_groups               <- names(table(CID45171_genelistmean_df$celltype_major))
results                   <- data.frame(CellType = character(), p.value = numeric(), stringsAsFactors = FALSE)
for (i in 1:length(cell_groups)) {
  CID45171_genelistmean    <- subset(CID45171_genelistmean_df,celltype_major==cell_groups[i])
  p_value                 <- wilcox.test(CID45171_genelistmean$gene_mean ~ CID45171_genelistmean$gene_panel, paired = TRUE)
  result_row              <- data.frame(CellType = cell_groups[i], p.value = p_value$p.value, stringsAsFactors = FALSE)
  results                 <- rbind(results, result_row)
}

CID45171_PAM50 <- subset(CID45171_genelistmean_df,gene_panel=="PAM50_gene")
CID45171_PAM50 <- subset(CID45171_PAM50,celltype_major=="CAFs"|celltype_major=="Cancer Epithelial")
wilcox.test(gene_mean~celltype_major,data = CID45171_PAM50)
p-value = 0.0372



##########################################################################
##### 3.2 bulk RNAseq data of mouse model and human cell lines
##########################################################################


####3.2.1: compare the genefu predicting and UBS93 predicting with bulk RNAseq data of humnan cell lines(fig_3c)####
####*GSE212143####
####pam50 subtyping
library(genefu)
library(dplyr)
data("pam50")
pam50.gene.df                  <- read.delim("Pro_TNBC/data/TCGA/TCGA/pam50_entrez_gene_id_and_ensemble_gene_id.txt", stringsAsFactors=FALSE,comment.char = '#') 
pam50.gene                     <- pam50.gene.df$ensemble.gene.id %>% as.character
colnames(pam50.gene.df)[1]     <- 'EntrezGene.ID'
colnames(pam50.gene.df)[2]     <- 'probe' 
pam50.gene.df$EntrezGene.ID    <- as.character(pam50.gene.df$EntrezGene.ID)

sampleID                       <- colnames(GSE212143.log2tpm)
pam50.gene.expr                <- GSE212143.log2tpm[pam50.gene.df$probe %>% as.character,sampleID] %>% t 
annot.matrix                   <- pam50.gene.df[,1:2] %>% as.matrix()
rownames(annot.matrix)         <- annot.matrix[,'probe']
pam50.subtype                  <- intrinsic.cluster.predict(sbt.model = pam50,data = pam50.gene.expr,annot = annot.matrix,mapping=annot.matrix,do.mapping = TRUE )
GSE212143.pam50.subtype                 <- pam50.subtype["subtype"] %>% as.data.frame()
GSE212143.pam50.subtype$run_accession   <- rownames(GSE212143.pam50.subtype)

####*Claudin low typing#
library(genefu)
data(claudinLowData)
test                      <- GSE212143.log2tpm
train                     <- claudinLowData
train$xd                  <-medianCtr(train$xd)
train_exp                 <- train$xd  %>% as.data.frame()
train_exp$ENTREZID        <- rownames(train_exp)
library(tidyverse)
library(clusterProfiler)
name                      <- bitr(rownames(train_exp),fromType = 'ENTREZID',toType = 'ENSEMBL',OrgDb = 'org.Hs.eg.db')
table(rownames(test) %in% name$ENSEMBL)
gene.name                 <- name[name$ENSEMBL %in% rownames(test),]
train_exp                 <- merge(gene.name,train_exp,by="ENTREZID")
rown                      <- train_exp$ENSEMBL
train_exp                 <- train_exp[,-c(1:2)]
rownames(train_exp)       <- rown
test                      <- test[rownames(test)%in% rownames(train_exp),]
predout                   <- claudinLow(x=train_exp, classes=as.matrix(train$classes$Group,ncol=1), y=test)
results                   <- cbind(predout$predictions, predout$distances)
colnames(results)[1]      <- "run_accession"
GSE212143.all.subtype     <- merge(GSE212143.pam50.subtype,results,by="run_accession")
GSE212143.all.subtype     <- merge(GSE212143.all.subtype,GSE212143.subtype,by="run_accession")
GSE212143.all.subtype     <- GSE212143.all.subtype[,c(1:3,6:8,4:5,9:11)]

####**GSE48213####
#####*pam50  subtyping#
library(genefu)
library(dplyr)
data("pam50")
pam50.gene.df                  <- read.delim("Pro_TNBC/data/TCGA/TCGA/pam50_entrez_gene_id_and_ensemble_gene_id.txt", stringsAsFactors=FALSE,comment.char = '#') 
pam50.gene                     <- pam50.gene.df$ensemble.gene.id %>% as.character
colnames(pam50.gene.df)[1]     <- 'EntrezGene.ID'
colnames(pam50.gene.df)[2]     <- 'probe' 
pam50.gene.df$EntrezGene.ID    <- as.character(pam50.gene.df$EntrezGene.ID)

sampleID                       <- colnames(GSE48213_log2tpm)
pam50.gene.expr                <- GSE48213_log2tpm[pam50.gene.df$probe %>% as.character,sampleID] %>% t 
annot.matrix                   <- pam50.gene.df[,1:2] %>% as.matrix()
rownames(annot.matrix)         <- annot.matrix[,'probe']
pam50.subtype                  <- intrinsic.cluster.predict(sbt.model = pam50,data = pam50.gene.expr,annot = annot.matrix,mapping=annot.matrix,do.mapping = TRUE )
GSE48213.pam50.subtype             <- pam50.subtype["subtype"] %>% as.data.frame()
GSE48213.pam50.subtype$run_accession   <- rownames(GSE48213.pam50.subtype)

#####*Claudin low typing
library(genefu)
data(claudinLowData)
test                <- GSE48213_log2tpm
train               <- claudinLowData
train$xd            <- medianCtr(train$xd)
train_exp           <- train$xd  %>% as.data.frame()
train_exp$ENTREZID  <- rownames(train_exp)
library(tidyverse)
library(clusterProfiler)
name                  <- bitr(rownames(train_exp),fromType = 'ENTREZID',toType = 'ENSEMBL',OrgDb = 'org.Hs.eg.db')
table(rownames(test) %in% name$ENSEMBL)
gene.name             <- name[name$ENSEMBL %in% rownames(test),]
train_exp             <- merge(gene.name,train_exp,by="ENTREZID")
rown                  <- train_exp$ENSEMBL
train_exp             <- train_exp[,-c(1:2)]
rownames(train_exp)   <- rown
test                  <- test[rownames(test)%in% rownames(train_exp),]
predout               <- claudinLow(x=train_exp, classes=as.matrix(train$classes$Group,ncol=1), y=test)
results               <- cbind(predout$predictions, predout$distances)
colnames(results)[1]  <- "run_accession"
GSE48213.all.subtype     <- merge(GSE48213.pam50.subtype,results,by="run_accession")
GSE48213.all.subtype     <- merge(GSE48213.all.subtype,GSE48213.subtype,by="run_accession")
GSE48213.all.subtype     <- GSE48213.all.subtype[,c(1:3,6:8,4:5,9:12)]

####*GSE73526####
#####*pam50  subtyping#
library(genefu)
library(dplyr)
data("pam50")
pam50.gene.df                  <- read.delim("Pro_TNBC/data/TCGA/TCGA/pam50_entrez_gene_id_and_ensemble_gene_id.txt", stringsAsFactors=FALSE,comment.char = '#') 
pam50.gene                     <- pam50.gene.df$ensemble.gene.id %>% as.character
colnames(pam50.gene.df)[1]     <- 'EntrezGene.ID'
colnames(pam50.gene.df)[2]     <- 'probe' 
pam50.gene.df$EntrezGene.ID    <- as.character(pam50.gene.df$EntrezGene.ID)

sampleID                       <- colnames(GSE73526_log2fpkm)
pam50.gene.expr                <- GSE73526_log2fpkm[pam50.gene.df$probe %>% as.character,sampleID] %>% t 
annot.matrix                   <- pam50.gene.df[,1:2] %>% as.matrix()
rownames(annot.matrix)         <- annot.matrix[,'probe']
pam50.subtype                  <- intrinsic.cluster.predict(sbt.model = pam50,data = pam50.gene.expr,annot = annot.matrix,mapping=annot.matrix,do.mapping = TRUE )
GSE73526.pam50.subtype             <- pam50.subtype["subtype"] %>% as.data.frame()
GSE73526.pam50.subtype$sample.id   <- rownames(GSE73526.pam50.subtype)

####*Claudin low typing#
library(genefu)
data(claudinLowData)
test                <- GSE73526_log2fpkm
train               <- claudinLowData
train$xd            <-  medianCtr(train$xd)
train_exp           <- train$xd  %>% as.data.frame()
train_exp$ENTREZID  <- rownames(train_exp)
library(tidyverse)
library(clusterProfiler)
name                  <- bitr(rownames(train_exp),fromType = 'ENTREZID',toType = 'ENSEMBL',OrgDb = 'org.Hs.eg.db')
table(rownames(test) %in% name$ENSEMBL)
gene.name             <- name[name$ENSEMBL %in% rownames(test),]
train_exp             <- merge(gene.name,train_exp,by="ENTREZID")
rown                  <- train_exp$ENSEMBL
train_exp             <- train_exp[,-c(1:2)]
rownames(train_exp)   <- rown
test                  <- test[rownames(test)%in% rownames(train_exp),]
predout               <- claudinLow(x=train_exp, classes=as.matrix(train$classes$Group,ncol=1), y=test)
results               <- cbind(predout$predictions, predout$distances)
colnames(results)[1]  <- "sample.id"
GSE73526.all.subtype     <- merge(GSE73526.pam50.subtype,results,by="sample.id")
GSE73526.all.subtype     <- merge(GSE73526.all.subtype,GSE73526.subtype,by="sample.id")
GSE73526.all.subtype     <- GSE73526.all.subtype[,c(1:3,6:8,4:5,9:11)]


####*brca.ccle####
#####*pam50  subtyping#
library(genefu)
library(dplyr)
data("pam50")
pam50.gene.df                  <- read.delim("Pro_TNBC/data/TCGA/TCGA/pam50_entrez_gene_id_and_ensemble_gene_id.txt", stringsAsFactors=FALSE,comment.char = '#') 
pam50.gene                     <- pam50.gene.df$ensemble.gene.id %>% as.character
colnames(pam50.gene.df)[1]     <- 'EntrezGene.ID'
colnames(pam50.gene.df)[2]     <- 'probe' 
pam50.gene.df$EntrezGene.ID    <- as.character(pam50.gene.df$EntrezGene.ID)

sampleID                       <- colnames(brca_rpkm)
pam50.gene.expr                <- brca_rpkm[pam50.gene.df$probe %>% as.character,sampleID] %>% t 
annot.matrix                   <- pam50.gene.df[,1:2] %>% as.matrix()
rownames(annot.matrix)         <- annot.matrix[,'probe']
pam50.subtype                  <- intrinsic.cluster.predict(sbt.model = pam50,data = pam50.gene.expr,annot = annot.matrix,mapping=annot.matrix,do.mapping = TRUE )
brca.ccle.pam50.subtype             <- pam50.subtype["subtype"] %>% as.data.frame()
brca.ccle.pam50.subtype$CCLE_Name   <- rownames(brca.ccle.pam50.subtype)

#####*Claudin low typing#
library(genefu)
data(claudinLowData)
test                <- brca_rpkm
train               <- claudinLowData
train$xd            <-  medianCtr(train$xd)
train_exp           <- train$xd  %>% as.data.frame()
train_exp$ENTREZID  <- rownames(train_exp)
library(tidyverse)
library(clusterProfiler)
name                  <- bitr(rownames(train_exp),fromType = 'ENTREZID',toType = 'ENSEMBL',OrgDb = 'org.Hs.eg.db')
table(rownames(test) %in% name$ENSEMBL)
gene.name             <- name[name$ENSEMBL %in% rownames(test),]
train_exp             <- merge(gene.name,train_exp,by="ENTREZID")
rown                  <- train_exp$ENSEMBL
train_exp             <- train_exp[,-c(1:2)]
rownames(train_exp)   <- rown
test                  <- test[rownames(test)%in% rownames(train_exp),]
predout               <- claudinLow(x=train_exp, classes=as.matrix(train$classes$Group,ncol=1), y=test)
results                   <- cbind(predout$predictions, predout$distances)
colnames(results)[1]      <- "CCLE_Name"
brca.ccle.all.subtype     <- merge(brca.ccle.pam50.subtype,results,by="CCLE_Name")
brca.ccle.all.subtype     <- merge(brca.ccle.all.subtype,brcaccle.subtype,by="CCLE_Name")
brca.ccle.all.subtype     <- brca.ccle.all.subtype[,c(1:3,6:8,4:5,9:10)]


library(ggpubr)
library(readxl)
bulk_RNAseq_predicting          <- read_excel("Pro_TNBC/paper/data/results/section_3/bulk_RNAseq_predicting.xls")
basal_predicting                <- subset(bulk_RNAseq_predicting,subtype=="Basal")
fig_S9a_basal                   <- ggplot(basal_predicting, aes(x = predictor, y = Frequency)) +
  geom_point(size = 6,color='#6600FF') + 
  geom_point(size = 6,shape=21) +  #绘制散点
  geom_line(aes(group = dataset), color = 'black', lwd = 0.5) +
  labs(x = '', y = 'Frequency')+ggplot.style
ggsave(fig_S9a_basal,filename = "Pro_TNBC/paper/plot/section_3/fig_S9a_basal_compear.pdf",height = 15,width = 20)

cl_predicting               <- subset(bulk_RNAseq_predicting,subtype=="Claudin-low")
fig_S9b_cl                   <- ggplot(cl_predicting, aes(x = predictor, y = Frequency)) +
  geom_point(size = 6,color='#6600FF') + 
  geom_point(size = 6,shape=21) +  #绘制散点
  geom_line(aes(group = dataset), color = 'black', lwd = 0.5) +
  labs(x = '', y = 'Frequency')+ggplot.style
ggsave(fig_S9b_cl,filename = "Pro_TNBC/paper/plot/section_3/fig_S9b_cl_compear.pdf",height = 15,width = 20)

Luminal_predicting <- subset(bulk_RNAseq_predicting,subtype=="Luminal")
fig_S9c_Luminal                   <- ggplot(Luminal_predicting, aes(x = predictor, y = Frequency)) +
  geom_point(size = 6,color='#6600FF') + 
  geom_point(size = 6,shape=21) +  #绘制散点
  geom_line(aes(group = dataset), color = 'black', lwd = 0.5) +
  labs(x = '', y = 'Frequency')+ggplot.style
ggsave(fig_S9c_Luminal,filename = "Pro_TNBC/paper/plot/section_3/fig_S9c_Luminal_compear.pdf",height = 15,width = 20)

HER2_predicting <- subset(bulk_RNAseq_predicting,subtype=="HER2-amp")
fig_S9d_HER2                   <- ggplot(HER2_predicting, aes(x = predictor, y = Frequency)) +
  geom_point(size = 6,color='#6600FF') + 
  geom_point(size = 6,shape=21) +  #绘制散点
  geom_line(aes(group = dataset), color = 'black', lwd = 0.5) +
  labs(x = '', y = 'Frequency')+ggplot.style
ggsave(fig_S9d_HER2,filename = "Pro_TNBC/paper/plot/section_3/fig_S9d_HER2_compear.pdf",height = 15,width = 20)

prediction               <-  data.frame(dataset=c("GSE212143","GSE48213","GSE73526","CCLE","GSE212143","GSE48213","GSE73526","CCLE"),
                          Frequency=c(0.60,0.741,0.774,0.708,0.933, 0.850, 0.849, 0.958),
                          predictor=c(rep("PAM50",4),rep("UBS93",4)))
fig_3c                   <- ggplot(prediction, aes(x = predictor, y = Frequency))+
  geom_point(size = 6,color='#6600FF') + 
  geom_point(size = 6,shape=21) +  
  geom_line(aes(group = dataset), color = 'black', lwd = 0.5) +
  labs(x = '', y = 'Frequency')+ggplot.style
ggsave(fig_3c,filename = "Pro_TNBC/paper/plot/section_3/fig_3c_total_compear.pdf",height = 15,width = 20)




##################################################################################
##### 3.3 compare UBS93 with SCsubtype using the scRNAseq data of cell lines
##################################################################################
library(readxl)
library(readr)
SCsubtype_signature_genes    <- read_excel("Pro_TNBC/data/scRNASeq/26/SCsubtype_signature_genes.xlsx")
metadata                     <- read_csv("Pro_TNBC/data/scRNASeq/26/Wu_etal_2021_BRCA_scRNASeq/metadata.csv")
tumor_metadata               <- subset(metadata,celltype_major=="Cancer Epithelial")
sampleID                     <- names(table(tumor_metadata$orig.ident))
load("Pro_TNBC/output/data/scRNASeq/26_sample/GSE176078_scRNA.RData")
scRNA$cell.id                <- rownames(scRNA@meta.data)
GSE176068_tumor_scRNA        <- subset(scRNA,cell.id %in% tumor_metadata$...1)
Tumour_ID                    <- c("CID3948","CID4290A","CID4530N","CID4535","CID3921","CID45171","CID4495","CID44971","CID44991","CID4515","CID3921")
GSE176068_tumor_train_scRNA  <- subset(GSE176068_tumor_scRNA,orig.ident %in% Tumour_ID)
save(GSE176068_tumor_train_scRNA,file = "Pro_TNBC/paper/data/results/section_3/GSE176068_tumor_train_scRNA.RData")
######*single cell RNAseq of cell lines(GSE173634)####
library(Seurat)
library(foreach)
load("Pro_TNBC/output/data/CCLE/GSE173634_scRNA.RData")
load("~/Pro_TNBC/output/data/CCLE/GSE173634_info.RData")
load("Pro_TNBC/paper/data/results/section_3/GSE176068_tumor_train_scRNA.RData")
GSE173634_scRNA$cell.id       <- rownames(GSE173634_scRNA@meta.data)
GSE173634_exp                 <- as.matrix(GSE173634_scRNA@assays$RNA@counts) %>% as.data.frame()
gene_row                      <- rownames(GSE173634_exp)
symbol                        <- bitr(gene_row,fromType = "ENSEMBL",toType = "SYMBOL",OrgDb = "org.Hs.eg.db")
GSE173634_exp$ENSEMBL         <- rownames(GSE173634_exp)
GSE173634_df                  <- merge(symbol,GSE173634_exp,by="ENSEMBL")
index                         <- order(rowMeans(GSE173634_df[,-c(1,2)]),decreasing=T)
GSE173634_df_order            <- GSE173634_df[index,]
keep                          <- !duplicated(GSE173634_df_order$SYMBOL)
GSE173634_df_order            <- GSE173634_df_order[keep,]
table(duplicated(GSE173634_df_order$SYMBOL))
rown                                       <- GSE173634_df_order$SYMBOL
GSE173634_counts                           <- GSE173634_df_order[,-c(1,2)] %>% as.matrix()
rownames(GSE173634_counts)                 <- rown
GSE176068_tumor_train_exp                  <- as.matrix(GSE176068_tumor_train_scRNA@assays$RNA@counts)
common_genes                               <- intersect(rownames(GSE173634_counts),rownames(GSE176068_tumor_train_exp))
GSE173634_common_counts                    <- GSE173634_counts[common_genes,]
GSE176068_tumor_train_common_counts        <- GSE176068_tumor_train_exp[common_genes,]
identical(rownames(GSE173634_common_counts),rownames(GSE176068_tumor_train_common_counts))
GSE173634_common_scRNA                     <- CreateSeuratObject(counts = GSE173634_common_counts)
GSE176068_tumor_train_common_scRNA         <- CreateSeuratObject(counts = GSE176068_tumor_train_common_counts)
GSE173634_info_pam50 <- subset(GSE173634_info,lineage_molecular_subtype=="basal_A"|lineage_molecular_subtype=="luminal"|lineage_molecular_subtype=="HER2_amp")
temp_allgenes  <- c(as.vector(SCsubtype_signature_genes$Basal_SC),
                    as.vector(SCsubtype_signature_genes$Her2E_SC),
                    as.vector(SCsubtype_signature_genes$LumA_SC),
                    as.vector(SCsubtype_signature_genes$LumB_SC))
temp_allgenes  <- unique(temp_allgenes[!temp_allgenes == ""])
center_sweep   <- function(x, row.w = rep(1, nrow(x))/nrow(x)) {
  get_average  <- function(v) sum(v * row.w)/sum(row.w)
  average      <- apply(x, 2, get_average)
  sweep(x, 2, average)
}
GSE173634_subtype_result        <- NULL
for (i in 1:nrow(GSE173634_info_pam50)) {
  A                          <- GSE173634_info_pam50[i,3] %>% as.character()
  scRNA                      <- subset(GSE173634_common_scRNA, orig.ident == A)
  scRNA.combined             <- merge(scRNA, y = GSE176068_tumor_train_common_scRNA)
  scRNA.combined             <- NormalizeData(object = scRNA.combined)
  scRNA.combined             <- ScaleData(scRNA.combined, features=temp_allgenes)
  tocalc                     <- as.data.frame(scRNA.combined@assays$RNA@scale.data)
  #calculate mean scsubtype scores
  outdat <- matrix(0,
                   nrow=length(SCsubtype_signature_genes),
                   ncol=ncol(tocalc),
                   dimnames=list(names(SCsubtype_signature_genes),
                                 colnames(tocalc)))
  for(j in 1:length(SCsubtype_signature_genes)){
    row    <- as.character(SCsubtype_signature_genes[[j]])
    row    <- unique(row[row != ""])
    genes  <- which(rownames(tocalc) %in% row)
    temp   <- apply(tocalc[genes,],2,function(x){mean(as.numeric(x),na.rm=TRUE)})
    outdat[j,] <- as.numeric(temp)
  }
  final          <- outdat[which(rowSums(outdat,na.rm=TRUE)!=0),]
  final          <- as.data.frame(final)
  is.num         <- sapply(final, is.numeric)
  final[is.num]  <- lapply(final[is.num], round, 4)
  finalm         <- as.matrix(final)
  ##Obtaining the highest call
  finalmt                      <- as.data.frame(t(finalm))
  scale_scscore                <- center_sweep(finalmt)
  scsubtype                    <- colnames(scale_scscore)[max.col(scale_scscore,ties.method="first")]
  scale_scscore$SCSubtypeCall  <- scsubtype
  scale_scscore$org.ident      <- scRNA.combined$orig.ident
  test_scscore                 <- subset(scale_scscore,org.ident==A)
  test_Subtype                 <- table(test_scscore$SCSubtypeCall) %>% as.data.frame() 
  colnames(test_Subtype)[2]    <- GSE173634_info_pam50[i,3] 
  colnames(test_Subtype)[1]    <- "subtype"
  if (is.null(GSE173634_subtype_result)) {
    GSE173634_subtype_result   <- test_Subtype
  } else {
    GSE173634_subtype_result   <- merge(GSE173634_subtype_result, test_Subtype,by="subtype",all=T)
  }
}

rown                                    <- as.character(GSE173634_subtype_result$subtype)
GSE173634_subtype_result                <- GSE173634_subtype_result[,-1]    
rownames(GSE173634_subtype_result)      <- rown
GSE173634_subtype_result                <- as.matrix(GSE173634_subtype_result)  %>% t()  %>% as.data.frame()
GSE173634_subtype_result                <- cbind(GSE173634_info_pam50[,2],GSE173634_subtype_result)
GSE173634_subtype_result$total_numbers  <- apply(GSE173634_subtype_result[,2:5],1,function(x){sum(x,na.rm =T)})
GSE173634_subtype_result$Lum            <- GSE173634_subtype_result$LumA + GSE173634_subtype_result$LumB
GSE173634_subtype_basal                 <- subset(GSE173634_subtype_result,lineage_molecular_subtype=="basal_A")
GSE173634_subtype_basal$ratio           <- GSE173634_subtype_basal$Basal_SC/GSE173634_subtype_basal$total_numbers
GSE173634_subtype_her2                  <- subset(GSE173634_subtype_result,lineage_molecular_subtype=="HER2_amp")
GSE173634_subtype_her2$ratio            <- GSE173634_subtype_her2$Her2E_SC/GSE173634_subtype_her2$total_numbers
GSE173634_subtype_luminal               <- subset(GSE173634_subtype_result,lineage_molecular_subtype=="luminal")
GSE173634_subtype_luminal$ratio         <- GSE173634_subtype_luminal$Lum/GSE173634_subtype_luminal$total_numbers
GSE173634_subtype                       <- rbind(GSE173634_subtype_basal,GSE173634_subtype_her2)
GSE173634_subtype                       <- rbind(GSE173634_subtype,GSE173634_subtype_luminal)
GSE173634_subtype$celllines_name        <- rownames(GSE173634_subtype)
GSE173634_subtype                       <- arrange(GSE173634_subtype,lineage_molecular_subtype)
GSE173634_subtype[GSE173634_subtype$lineage_molecular_subtype=="basal_A",1]  <- c(rep("Basal",11))
GSE173634_subtype[GSE173634_subtype$lineage_molecular_subtype=="luminal",1]  <- c(rep("Luminal",8))
GSE173634_SCsubtype_subtype                         <- GSE173634_subtype
GSE173634_SCsubtype_subtype[is.na(GSE173634_SCsubtype_subtype$ratio),8] <- 0
save(GSE173634_SCsubtype_subtype,file = "Pro_TNBC/paper/data/results/section_3/GSE173634_SCsubtype_subtype.RData")
load("~/Pro_TNBC/paper/data/results/section_2/fig2b.GSE17634.SUBTYPE.RData")
load("Pro_TNBC/paper/data/results/section_3/GSE173634_SCsubtype_subtype.RData")
GSE173634_UBS93_subtype      <- GSE173634_subtype[-c(17:20),c(1,7,8)] %>% mutate(classify="UBS93")
GSE173634_SCsubtype_subtype  <- GSE173634_SCsubtype_subtype[,c(1,8,9)] %>% mutate(classify="SCSubtype")
identical(GSE173634_UBS93_subtype$celllines_name,GSE173634_SCsubtype_subtype$celllines_name)
GSE173634_subtype              <- rbind(GSE173634_SCsubtype_subtype,GSE173634_UBS93_subtype)
fig_3c                         <- ggplot(GSE173634_subtype,aes(x=GSE173634_subtype$lineage_molecular_subtype,y=ratio,fill=classify))+
  geom_boxplot(width=0.6,alpha=1#the transparency 
               ,position=position_dodge(0.6)#determines the position of the boxes in relation to each other
               ,outlier.shape = NA)+ggplot.style + scale_fill_manual(values =c(SCSubtype="Blue",UBS93="Red"))+
  geom_point(size=2,position = position_jitterdodge(dodge.width=0.6))+
  theme(panel.grid.major=element_line(colour=NA),
        panel.grid.minor = element_blank())+xlab("Subtype")+ylab("Frequency")+
  scale_x_discrete("Subtype", labels = c("HER2_amp" = paste("HER2_amp",sep = "\n","(n=5)"),"Basal" = paste("Basal",sep = "\n","(n=11)"),"Luminal" = paste("Luminal",sep = "\n","(n=8)")))
ggsave(fig_3c,filename = "Pro_TNBC/paper/plot/section_3/fig_3c.GSE173634.compare.subtype.pdf",width = 20,height = 15)

####*single cell RNAseq of cell lines(GSE202771)#####
load("Pro_TNBC/output/data/scRNASeq/GSE202771/GSE202771_scRNA.RData")
load("Pro_TNBC/data/CCLE/scRNAseq/GSE202771_infor.RData")
common_genes                            <- intersect(rownames(GSE202771_scRNA),rownames(GSE176068_tumor_train_scRNA))
GSE202771_counts                        <- as.matrix(GSE202771_scRNA@assays$RNA@counts)
GSE202771_common_counts                 <- GSE202771_counts[common_genes,]
GSE176068_tumor_train_counts            <- as.matrix(GSE176068_tumor_train_scRNA@assays$RNA@counts)
GSE176068_tumor_train_common_counts     <- GSE176068_tumor_train_counts[common_genes,]
identical(rownames(GSE202771_common_counts),rownames(GSE176068_tumor_train_common_counts))
GSE202771_common_scRNA                  <- CreateSeuratObject(counts = GSE202771_common_counts)
GSE202771_common_scRNA$orig.ident       <- GSE202771_scRNA$orig.ident
GSE176068_tumor_train_common_scRNA      <- CreateSeuratObject(counts = GSE176068_tumor_train_common_counts)
temp_allgenes           <- c(as.vector(SCsubtype_signature_genes$Basal_SC),
                             as.vector(SCsubtype_signature_genes$Her2E_SC),
                             as.vector(SCsubtype_signature_genes$LumA_SC),
                             as.vector(SCsubtype_signature_genes$LumB_SC))
temp_allgenes                   <- unique(temp_allgenes[!temp_allgenes == ""])
GSE202771_infor_pam50           <- subset(GSE202771_infor,lineage_molecular_subtype=="basal_A"|lineage_molecular_subtype=="luminal")
GSE202771_subtype_result        <- NULL
for (i in 1:nrow(GSE202771_infor_pam50)) {
  A                          <- GSE202771_infor_pam50[i,1] %>% as.character()
  scRNA                      <- subset(GSE202771_common_scRNA, orig.ident == A)
  scRNA.combined             <- merge(scRNA,GSE176068_tumor_train_common_scRNA)
  scRNA.combined             <- NormalizeData(object = scRNA.combined)
  scRNA.combined             <- ScaleData(scRNA.combined, features=temp_allgenes)
  tocalc                     <- as.data.frame(scRNA.combined@assays$RNA@scale.data)
  #calculate mean scsubtype scores
  outdat <- matrix(0,
                   nrow=length(SCsubtype_signature_genes),
                   ncol=ncol(tocalc),
                   dimnames=list(names(SCsubtype_signature_genes),
                                 colnames(tocalc)))
  for(j in 1:length(SCsubtype_signature_genes)){
    row    <- as.character(SCsubtype_signature_genes[[j]])
    row    <- unique(row[row != ""])
    genes  <- which(rownames(tocalc) %in% row)
    temp   <- apply(tocalc[genes,],2,function(x){mean(as.numeric(x),na.rm=TRUE)})
    outdat[j,] <- as.numeric(temp)
  }
  final          <- outdat[which(rowSums(outdat,na.rm=TRUE)!=0),]
  final          <- as.data.frame(final)
  is.num         <- sapply(final, is.numeric)
  final[is.num]  <- lapply(final[is.num], round, 4)
  finalm         <- as.matrix(final)
  ##Obtaining the highest call
  finalmt                      <- as.data.frame(t(finalm))
  scale_scscore                <- center_sweep(finalmt)
  scsubtype                    <- colnames(scale_scscore)[max.col(scale_scscore,ties.method="first")]
  scale_scscore$SCSubtypeCall  <- scsubtype
  scale_scscore$org.ident      <- scRNA.combined$orig.ident
  test_scscore                 <- subset(scale_scscore,org.ident==A)
  test_Subtype                 <- table(test_scscore$SCSubtypeCall) %>% as.data.frame() 
  colnames(test_Subtype)[2]    <- GSE202771_infor_pam50[i,1] 
  colnames(test_Subtype)[1]    <- "subtype"
  if (is.null(GSE202771_subtype_result)) {
    GSE202771_subtype_result <- test_Subtype
  } else {
    GSE202771_subtype_result <- merge(GSE202771_subtype_result, test_Subtype,by="subtype",all=T)
  }
}

rown                                    <- as.character(GSE202771_subtype_result$subtype)
GSE202771_subtype_result                <- GSE202771_subtype_result[,-1]    
rownames(GSE202771_subtype_result)      <- rown
GSE202771_subtype_result                <- as.matrix(GSE202771_subtype_result)  %>% t()  %>% as.data.frame()
GSE202771_subtype_result                <- cbind(GSE202771_infor_pam50[,4],GSE202771_subtype_result)
GSE202771_subtype_result$total_numbers  <- apply(GSE202771_subtype_result[,2:5],1,function(x){sum(x,na.rm =T)})
colnames(GSE202771_subtype_result)[1]   <- "Subtype"
GSE202771_subtype_result$Lum            <- GSE202771_subtype_result$LumA + GSE202771_subtype_result$LumB
GSE202771_subtype_basal                 <- subset(GSE202771_subtype_result,Subtype=="basal_A")
GSE202771_subtype_basal$ratio           <- GSE202771_subtype_basal$Basal_SC/GSE202771_subtype_basal$total_numbers
GSE202771_subtype_luminal                   <- subset(GSE202771_subtype_result,Subtype=="luminal")
GSE202771_subtype_luminal$ratio             <- GSE202771_subtype_luminal$Lum/GSE202771_subtype_luminal$total_numbers
GSE202771_subtype                           <- rbind(GSE202771_subtype_basal,GSE202771_subtype_luminal)
GSE202771_subtype$celllines_name            <- rownames(GSE202771_subtype)
GSE202771_subtype[GSE202771_subtype$Subtype=="basal_A",1]  <- c(rep("Basal",7))
GSE202771_subtype[GSE202771_subtype$Subtype=="luminal",1]  <- c(rep("Luminal",2))
GSE202771_SCsubtype_subtype                                <- GSE202771_subtype
save(GSE202771_SCsubtype_subtype,file = "Pro_TNBC/paper/data/results/section_3/GSE202771_SCsubtype_subtype.RData")
load("~/Pro_TNBC/paper/data/results/section_2/supfig_4.GSE202771.SUBTYPE.RData")
load("Pro_TNBC/paper/data/results/section_3/GSE173634_SCsubtype_subtype.RData")
GSE202771_UBS93_subtype              <- GSE202771_subtype[-c(8:13),c(1,8,9)] %>% mutate(classify="UBS93")
GSE202771_SCsubtype_subtype          <- GSE202771_SCsubtype_subtype[,c(1,8,9)] %>% mutate(classify="SCSubtype")
identical(GSE202771_UBS93_subtype$celllines_name,GSE202771_SCsubtype_subtype$celllines_name)
colnames(GSE202771_UBS93_subtype)[1] <- "Subtype"
GSE202771_subtype                    <- rbind(GSE202771_SCsubtype_subtype,GSE202771_UBS93_subtype)
colnames(GSE173634_subtype)[1]       <- "Subtype"
scRNAseq_celllines_compare           <- rbind(GSE173634_subtype,GSE202771_subtype)
fig_3d                               <- ggplot(scRNAseq_celllines_compare,aes(x=Subtype,y=ratio,color=classify))+
  geom_boxplot(width=0.6,alpha=1#the transparency 
               ,position=position_dodge(0.8)#determines the position of the boxes in relation to each other
               ,outlier.shape = NA)+ggplot.style +
  scale_color_manual(values =c(SCSubtype="Blue",UBS93="Red"))+
  geom_point(size=5,position = position_jitterdodge(dodge.width=0.8))+
  theme(panel.grid.major=element_line(colour=NA),
        panel.grid.minor = element_blank())+xlab("Subtype")+ylab("Frequency")+
  scale_x_discrete("Subtype", labels = c("HER2_amp" = paste("HER2_amp",sep = "\n","(n=5)"),"Basal" = paste("Basal",sep = "\n","(n=18)"),"Luminal" = paste("Luminal",sep = "\n","(n=10)")))
ggsave(fig_3d,filename = "Pro_TNBC/paper/plot/section_3/fig_3d.compare.subtype.pdf",width = 20,height = 15)
save(scRNAseq_celllines_compare,file = "Pro_TNBC/paper/data/results/section_3/scRNAseq_celllines_compare.RData")

scRNAseq_celllines_compare_basal     <- subset(scRNAseq_celllines_compare,Subtype=="Basal")
wilcox.test( ratio~classify,scRNAseq_celllines_compare_basal,paired=T)#p-value = 0.000145
scRNAseq_celllines_compare_HER2      <- subset(scRNAseq_celllines_compare,Subtype=="HER2_amp")
wilcox.test( ratio~classify,scRNAseq_celllines_compare_HER2,paired=T)#p-value = 0.125
scRNAseq_celllines_compare_luminal   <- subset(scRNAseq_celllines_compare,Subtype=="Luminal")
wilcox.test( ratio~classify,scRNAseq_celllines_compare_luminal,paired=T)#p-value = 0.8457



