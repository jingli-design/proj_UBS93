library(readr)
library(org.Hs.eg.db)
library(clusterProfiler)
load("Pro_TNBC/data/CCLE/CCLE.RData")
load("Pro_TNBC/paper/data/results/section_1/UBS93.gene.df")
####calculate the centroid matrix of the UBS93 gene set using median####
ccle_info         <- read.csv("Pro_TNBC/data/CCLE/ccle_info.csv")
brca_meta         <- CCLE.sample.meta[CCLE.sample.meta$cell.line.tumor.site=="BREAST",]
brca_info         <- ccle_info[ccle_info$CCLE_Name  %in% brca_meta$cell.line.name,]
brca_info         <- brca_info[brca_info$primary_disease=="Breast Cancer",]
brca_info         <- brca_info[-2,]
basalA.ccle       <- subset(brca_info,brca_info$lineage_molecular_subtype == "basal_A") 
basalA.ccle.id    <- basalA.ccle$CCLE_Name
basalB.ccle       <- subset(brca_info,brca_info$lineage_molecular_subtype == "basal_B") 
basalB.ccle.id    <- basalB.ccle$CCLE_Name
nonbasal.ccle.id  <- setdiff(brca_info$CCLE_Name,basalA.ccle.id)
nonbasal.ccle.id  <- setdiff(nonbasal.ccle.id,basalB.ccle.id)


ccle_UBS93_rpkm                 <- CCLE.log2.rpkm.matrix[rownames(CCLE.log2.rpkm.matrix) %in% UBS93.gene.df$ENSEMBL,]
brcaccle_UBS93_rpkm             <- ccle_UBS93_rpkm[,colnames(ccle_UBS93_rpkm) %in% brca_info$CCLE_Name]
brcaccle_UBS93_rpkm_rank        <- apply(brcaccle_UBS93_rpkm , 1,rank )  
Basal                           <- apply(brcaccle_UBS93_rpkm_rank[rownames(brcaccle_UBS93_rpkm_rank) %in% basalA.ccle.id,],2,median) %>% as.data.frame()
Basal$Claudin_low               <- apply(brcaccle_UBS93_rpkm_rank[rownames(brcaccle_UBS93_rpkm_rank) %in% basalB.ccle.id,],2,median)
Basal$H_or_L                    <- apply(brcaccle_UBS93_rpkm_rank[rownames(brcaccle_UBS93_rpkm_rank) %in% nonbasal.ccle.id,],2,median)
UBS93_median_centriod                  <- Basal
colnames(UBS93_median_centriod)[1]     <- "Basal"
save(UBS93_median_centriod,file="Pro_TNBC/paper/data/method/UBS93_median_centriod.RData")





