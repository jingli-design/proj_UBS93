########################################
####the plots for the method section####
########################################

#### Supplementary Figure 10: Scatter plot illustrating the relationship between the average E-M ratio and the rank of 27,719 genes (a), the top 2,720 genes (b), and the top 500 genes (c)####
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



####Supplementary Figure 11: ROC curve plotted based on ssGSEA scores of the combined dataset (GSE212143 and GSE48213).#####
library(readr)
library(GSVA)
load("Pro_TNBC/data/CCLE/bulk_RNAseq/cell_lines/processed_by_myself/GSE212143.gene.expression.RData")
GSE212143.log2tpm                        <- log2.tpm.matrix
GSE212143_info                           <- read_csv("Pro_TNBC/data/CCLE/bulk_RNAseq/cell_lines/processed_by_myself/GSE212143_info.csv")
load("Pro_TNBC/data/CCLE/bulk_RNAseq/cell_lines/processed_by_myself/GSE48213.gene.expression.RData")
GSE48213_info                            <- read_csv("Pro_TNBC/data/CCLE/bulk_RNAseq/cell_lines/processed_by_myself/GSE48213_info.csv")
colnames(ccle_info)[3]                   <- "cell_lines"
GSE48213_info                            <- merge(GSE48213_info,ccle_info[,c(3,21)],by="cell_lines")
write.csv(GSE48213_info,file = "Pro_TNBC/data/CCLE/bulk_RNAseq/cell_lines/processed_by_myself/GSE48213_info.csv")
GSE48213_log2tpm                         <- log2.tpm.matrix
GSE212143.UBS93.log2tpm                  <- GSE212143.log2tpm[rownames(GSE212143.log2tpm) %in% BCT93.gene.df$ENSEMBL,]
GSE48213.UBS93.log2tpm                   <- GSE48213_log2tpm[rownames(GSE48213_log2tpm) %in% BCT93.gene.df$ENSEMBL,]
identical(rownames(GSE48213.UBS93.log2tpm),rownames(GSE212143.UBS93.log2tpm))
GSE48213_horl                            <- subset(GSE48213_info,lineage_molecular_subtype=="HER2_amp"|lineage_molecular_subtype=="luminal")
GSE48213_horl                            <- GSE48213_horl[,c(4,6)]
GSE48213_horl_log2tpm                    <- GSE48213.UBS93.log2tpm[,colnames(GSE48213.UBS93.log2tpm) %in% GSE48213_horl$run_accession]
GSE212143_horl                           <- subset(GSE212143_info,molecular_subtype=="HER2_amp"|molecular_subtype=="luminal")
GSE212143_horl                           <- GSE212143_horl[,c(2,4)]
colnames(GSE212143_horl)[2]              <- "lineage_molecular_subtype"
GSE212143_horl_log2tpm                   <- GSE212143.UBS93.log2tpm[,colnames(GSE212143.UBS93.log2tpm) %in% GSE212143_horl$run_accession]
identical(rownames(GSE48213_horl_log2tpm),rownames(GSE212143_horl_log2tpm))
cell_lines_horl                          <- cbind(GSE48213_horl_log2tpm,GSE212143_horl_log2tpm)
cell_lines_horl_info                     <- rbind(GSE48213_horl,GSE212143_horl)
geneset                                <- list(h_up_l_gene=h_up_l_gene)
res                                    <- gsva(cell_lines_horl,
                                               geneset, method="ssgsea",
                                               mx.diff=T, verbose=FALSE,ssgsea.norm=FALSE) 
res                                    <- t(res)%>% as.data.frame()
res$run_accession                      <- rownames(res)
ccle_horl_res                          <- merge(res,cell_lines_horl_info,by="run_accession")
library(pROC)
ccle_horl_res$subtype                 <- ifelse(ccle_horl_res$lineage_molecular_subtype=="HER2_amp",1,0)
ccle_up_roc_obj                       <- roc(ccle_horl_res$subtype,ccle_horl_res$h_up_l_gene)
plot(ccle_up_roc_obj,col="red",
     legacy.axes=T,#change Y-axis format 
     print.auc=TRUE,#Display AUC area
     print.thres=TRUE,#Add breakpoints and 95% CI
     main="The ROC curve of HER2 up-regulated genes in CCLE")

p1  <- ggroc(ccle_up_roc_obj,legacy.axes=T,color="red")+ 
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), 
               color="darkgrey", linetype="dashed")+
  ggplot.style+
  theme(panel.grid.major=element_line(colour=NA),
        panel.grid.minor = element_blank())+scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))

ggsave(p1,filename = "Pro_TNBC/paper/plot/method/The.ROC.cure.plot.of.ssgsea.score.of.her2-amp.sigature.gene.pdf",height = 15,width = 20)
