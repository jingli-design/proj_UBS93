

library(readr)
library(DESeq2)
load("Pro_TNBC/data/CCLE/CCLE.RData")

####Differential gene expression analysis between HER2-amp and Luminal cell lines in CCLE####
#gene expression
ccle_info                               <- read_csv("Pro_TNBC/data/CCLE/ccle_info.csv")
brca_info                               <- ccle_info[ccle_info$primary_disease=="Breast Cancer",]
brca.ccle.meta                          <- subset(CCLE.sample.meta,cell.line.tumor.site=="BREAST")
brca_info                               <- brca_info[brca_info$CCLE_Name %in% brca.ccle.meta$cell.line.name,]
brca.ccle.log2.readcount                <- CCLE.log2.read.count.matrix[,colnames(CCLE.log2.read.count.matrix) %in% brca_info$CCLE_Name]
luminal.cell.id                         <- brca_info[brca_info$lineage_molecular_subtype=="luminal",]$CCLE_Name
HER2amp.ccle.id                         <- brca_info[brca_info$lineage_molecular_subtype=="HER2_amp",]$CCLE_Name
horl.ccle.id                            <- c(HER2amp.ccle.id,luminal.cell.id)
HorL_log2_readcount                     <- brca.ccle.log2.readcount[,colnames(brca.ccle.log2.readcount) %in% horl.ccle.id]
HorL_readcount                          <- 2^HorL_log2_readcount -1
HorL_readcount                          <- round(HorL_readcount,digits = 0)#Data rounding

#subgroup 
horl_subgroup                              <- data.frame(CCLE_Name=horl.ccle.id)
horl_subgroup                              <- within(horl_subgroup,{
  group                                    <- NA
  group[CCLE_Name %in% luminal.cell.id]    <- "luminal"
  group[CCLE_Name %in% HER2amp.ccle.id]    <- "her2amp"
})
rown                                       <- horl_subgroup$CCLE_Name
horl_subgroup                              <- data.frame(group=horl_subgroup[,-1])
rownames(horl_subgroup)                    <- rown

HorL_readcount                             <- HorL_readcount[,match(rownames(horl_subgroup),colnames(HorL_readcount))]
all(rownames(horl_subgroup) == colnames(HorL_readcount))

#Specify sample grouping information
group                         <- as.factor(horl_subgroup$group)
group = relevel(group, "luminal")

#strcture DEseqDataSet object
any(is.na(HorL_readcount))
dds                           <- DESeqDataSetFromMatrix(HorL_readcount, horl_subgroup, design= ~ group)

#Calculate the difference multiple
dds                           <- DESeq(dds)
res                           <- results(dds, contrast = c('group', 'her2amp','luminal'))
res = res[order(res$pvalue),]
head(res)
#
res[which(res$log2FoldChange >= 1 & res$padj < 0.05),'sig']          <- 'up'
res[which(res$log2FoldChange <= -1 & res$padj < 0.05),'sig']         <- 'down'
res [which(abs(res$log2FoldChange) <= 1 | res$padj >= 0.05),'sig']   <- 'none'
#Output differential gene summary table
diff_gene_deseq2_horl_subgroup          <- subset(res, sig %in% c('up', 'down')) %>% as.data.frame()
load("~/Pro_TNBC/output/data/CCLE/UBS93.Rpackage/UBS93.gene.df.RData")
diff_gene_deseq2_horl_subgroup$ENSEMBL  <- rownames(diff_gene_deseq2_horl_subgroup)
diff_gene_deseq2_horl_subgroup          <- merge(diff_gene_deseq2_horl_subgroup,UBS93.gene.df,by="ENSEMBL")
write.csv(diff_gene_deseq2_horl_subgroup,file="Pro_TNBC/output/data/CCLE/diff_gene_deseq2_horl_subgroup.csv")
h_up_l_gene                             <- subset(diff_gene_deseq2_horl_subgroup,sig=="up")
HER2_amp_signature                      <- h_up_l_gene$ENSEMBL
save(HER2_amp_signature,file="Pro_TNBC/paper/data/method/HER2_amp_signature.RData")
