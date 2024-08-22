
library(readr)
library(GSVA)
load("Pro_TNBC/paper/data/results/section_1/UBS93.gene.df.RData")
load("Pro_TNBC/paper/data/method/HER2_amp_signature.RData")

######draw a ROC curve and get cut-off#####
load("Pro_TNBC/data/CCLE/bulk_RNAseq/cell_lines/processed_by_myself/GSE212143.gene.expression.RData")
GSE212143.log2tpm                        <- log2.tpm.matrix
GSE212143_info                           <- read_csv("Pro_TNBC/data/CCLE/bulk_RNAseq/cell_lines/processed_by_myself/GSE212143_info.csv")
GSE48213_info                            <- read_csv("Pro_TNBC/data/CCLE/bulk_RNAseq/cell_lines/processed_by_myself/GSE48213_info.csv")
load("Pro_TNBC/data/CCLE/bulk_RNAseq/cell_lines/processed_by_myself/GSE48213.gene.expression.RData")
GSE48213_log2tpm                         <- log2.tpm.matrix
GSE212143.UBS93.log2tpm                  <- GSE212143.log2tpm[rownames(GSE212143.log2tpm) %in% UBS93.gene.df$ENSEMBL,]
GSE48213.UBS93.log2tpm                   <- GSE48213_log2tpm[rownames(GSE48213_log2tpm) %in% UBS93.gene.df$ENSEMBL,]
identical(rownames(GSE48213.UBS93.log2tpm),rownames(GSE212143.UBS93.log2tpm))
GSE48213_horl                            <- subset(GSE48213_info,lineage_molecular_subtype=="HER2_amp"|lineage_molecular_subtype=="luminal")
GSE48213_horl                            <- GSE48213_horl[,c(5,7)]
GSE48213_horl_log2tpm                    <- GSE48213.UBS93.log2tpm[,colnames(GSE48213.UBS93.log2tpm) %in% GSE48213_horl$run_accession]
GSE212143_horl                           <- subset(GSE212143_info,molecular_subtype=="HER2_amp"|molecular_subtype=="luminal")
GSE212143_horl                           <- GSE212143_horl[,c(2,4)]
colnames(GSE212143_horl)[2]              <- "lineage_molecular_subtype"
GSE212143_horl_log2tpm                   <- GSE212143.UBS93.log2tpm[,colnames(GSE212143.UBS93.log2tpm) %in% GSE212143_horl$run_accession]
identical(rownames(GSE48213_horl_log2tpm),rownames(GSE212143_horl_log2tpm))
cell_lines_horl                          <- cbind(GSE48213_horl_log2tpm,GSE212143_horl_log2tpm)
cell_lines_horl_info                     <- rbind(GSE48213_horl,GSE212143_horl)
geneset                                  <- list(h_up_l_gene=HER2_amp_signature)
res                                      <- gsva(cell_lines_horl,
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

library(readr)
library(GSVA)
load("Pro_TNBC/paper/data/results/section_1/UBS93.gene.df.RData")
load("Pro_TNBC/paper/data/method/HER2_amp_signature.RData")
load("Pro_TNBC/paper/data/method/UBS93_median_centriod.RData")
######draw a ROC curve and get cut-off#####
load("Pro_TNBC/data/CCLE/bulk_RNAseq/cell_lines/processed_by_myself/GSE212143.gene.expression.RData")
GSE212143.log2tpm                        <- log2.tpm.matrix
GSE212143_info                           <- read_csv("Pro_TNBC/data/CCLE/bulk_RNAseq/cell_lines/processed_by_myself/GSE212143_info.csv")
GSE48213_info                            <- read_csv("Pro_TNBC/data/CCLE/bulk_RNAseq/cell_lines/processed_by_myself/GSE48213_info.csv")
load("Pro_TNBC/data/CCLE/bulk_RNAseq/cell_lines/processed_by_myself/GSE48213.gene.expression.RData")
GSE48213_log2tpm                         <- log2.tpm.matrix
GSE212143.UBS93.log2tpm                  <- GSE212143.log2tpm[rownames(GSE212143.log2tpm) %in% UBS93.gene.df$ENSEMBL,]
GSE48213.UBS93.log2tpm                   <- GSE48213_log2tpm[rownames(GSE48213_log2tpm) %in% UBS93.gene.df$ENSEMBL,]
identical(rownames(GSE48213.UBS93.log2tpm),rownames(GSE212143.UBS93.log2tpm))
GSE48213_horl                            <- subset(GSE48213_info,lineage_molecular_subtype=="HER2_amp"|lineage_molecular_subtype=="luminal")
GSE48213_horl                            <- GSE48213_horl[,c(5,7)]
GSE48213_horl_log2tpm                    <- GSE48213.UBS93.log2tpm[,colnames(GSE48213.UBS93.log2tpm) %in% GSE48213_horl$run_accession]
GSE212143_horl                           <- subset(GSE212143_info,molecular_subtype=="HER2_amp"|molecular_subtype=="luminal")
GSE212143_horl                           <- GSE212143_horl[,c(2,4)]
colnames(GSE212143_horl)[2]              <- "lineage_molecular_subtype"
GSE212143_horl_log2tpm                   <- GSE212143.UBS93.log2tpm[,colnames(GSE212143.UBS93.log2tpm) %in% GSE212143_horl$run_accession]
identical(rownames(GSE48213_horl_log2tpm),rownames(GSE212143_horl_log2tpm))
cell_lines_horl                          <- cbind(GSE48213_horl_log2tpm,GSE212143_horl_log2tpm)
cell_lines_horl_info                     <- rbind(GSE48213_horl,GSE212143_horl)
geneset                                  <- list(h_up_l_gene=HER2_amp_signature)
res                                      <- gsva(cell_lines_horl,
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
UBS93.data      <- list(UBS93.centroid=UBS93.centroid,UBS93.gene.df=UBS93.gene.df,ssgsea.cutoff=ssgsea.cutoff,HER2.amp.signature.genes=HER2_amp_signature)
save(UBS93.data,file = "Pro_TNBC/paper/data/method/UBS93.data.RData")
