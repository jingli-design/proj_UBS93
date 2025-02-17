

################################################################################
####Find the differentiated genes between Basal and Claudin-low cells
################################################################################

####*DE analysis of basal and claudin-low cells in TNBC 0554####
library(Seurat)
Idents(TNBC_BRCA1_0554_tumor_scRNA)      <- "subtype"
de_results_TNBC0554                      <- FindMarkers(TNBC_BRCA1_0554_tumor_scRNA, slot="counts",ident.1 = "Claudin_low", ident.2 = "Basal",test.use = "DESeq2")
de_results_TNBC0554                      <- within(de_results_TNBC0554,{
  sig                                    <- NA
  sig[p_val_adj < 0.05 &avg_log2FC >0.5]   <- "up"
  sig[p_val_adj < 0.05 &avg_log2FC < -0.5] <- "down"
})
de_results_TNBC0554[is.na(de_results_TNBC0554$sig),6]      <- "none"
table(de_results_TNBC0554$sig)
de_results_TNBC0554$ENSEMBL               <- rownames(de_results_TNBC0554)
symbol                                    <- bitr(de_results_TNBC0554$ENSEMBL,fromType = "ENSEMBL",toType = "SYMBOL",OrgDb = "org.Hs.eg.db")
de_results_TNBC0554                       <- merge(symbol,de_results_TNBC0554,by="ENSEMBL")


####*HDQP1 in GSE173634####
Idents(HDQP1_scRNA)                      <- "subtype"
de_results_HDQP1_GSE173634               <- FindMarkers(HDQP1_scRNA,slot="counts",ident.1 = "Claudin_low",ident.2 = "Basal",test.use = "DESeq2")
de_results_HDQP1_GSE173634               <- within(de_results_HDQP1_GSE173634,{
  sig                                    <- NA
  sig[p_val_adj < 0.05 &avg_log2FC >0.5]   <- "up"
  sig[p_val_adj < 0.05 &avg_log2FC < -0.5] <- "down"
})
de_results_HDQP1_GSE173634[is.na(de_results_HDQP1_GSE173634$sig),6] <- "none"
table(de_results_HDQP1_GSE173634$sig)
de_results_HDQP1_GSE173634$ENSEMBL               <- rownames(de_results_HDQP1_GSE173634)
symbol                                           <- bitr(de_results_HDQP1_GSE173634$ENSEMBL,fromType = "ENSEMBL",toType = "SYMBOL",OrgDb = "org.Hs.eg.db")
de_results_HDQP1_GSE173634                       <- merge(symbol,de_results_HDQP1_GSE173634,by="ENSEMBL")



####*HDQP1 in GSE202271####
Idents(HDQP1_scRNA)                      <- "subtype"
de_results_HDQP1_GSE202771               <- FindMarkers(HDQP1_scRNA,slot="counts",ident.1 = "Claudin_low",ident.2 = "Basal",test.use = "DESeq2")
de_results_HDQP1_GSE202771               <- within(de_results_HDQP1_GSE202771,{
  sig                                    <- NA
  sig[p_val_adj < 0.05 &avg_log2FC >0.5]   <- "up"
  sig[p_val_adj < 0.05 &avg_log2FC < -0.5] <- "down"
})
de_results_HDQP1_GSE202771[is.na(de_results_HDQP1_GSE202771$sig),6] <- "none"
table(de_results_HDQP1_GSE202771$sig)
de_results_HDQP1_GSE202771$SYMBOL                <- rownames(de_results_HDQP1_GSE202771)
library(clusterProfiler)
symbol                                           <- bitr(de_results_HDQP1_GSE202771$SYMBOL,fromType = "SYMBOL",toType = "ENSEMBL",OrgDb = "org.Hs.eg.db")
de_results_HDQP1_GSE202771                       <- merge(de_results_HDQP1_GSE202771,symbol,by="SYMBOL")
symbol                                           <- bitr(de_results_HDQP1_GSE202771$SYMBOL,fromType = "SYMBOL",toType = "ENSEMBL",OrgDb = "org.Hs.eg.db")
de_results_HDQP1_GSE202771                       <- merge(symbol,de_results_HDQP1_GSE202771,by="SYMBOL")


####*intersected de genes####
de_HDQP1_claudinlow_basal_intersect                 <- merge(de_results_HDQP1_GSE173634,de_results_HDQP1_GSE202771,by="ENSEMBL")
de_HDQP1_claudinlow_basal_intersect                 <- de_HDQP1_claudinlow_basal_intersect[,c(1,2,4,7,8,11,14,15)]
colnames(de_HDQP1_claudinlow_basal_intersect)[2:8]  <- c("SYMBOL","log2FC_HDQP1_GSE173634","padj_HDQP1_GSE173634","sig_1","log2FC_HDQP1_GSE202771","padj_HDQP1_GSE202771","sig_2")
de_intersect_gene                     <- merge(de_HDQP1_claudinlow_basal_intersect,de_results_TNBC0554,by="ENSEMBL")
de_intersect_gene                     <- de_intersect_gene[,c(1:8,11,14,15)]
colnames(de_intersect_gene)[c(9:11)]  <- c("log2FC_TNBC0554","padj_TNBC0554","sig_0554")
de_intersect_gene                     <- within(de_intersect_gene,{
  sig                                                <- NA
  sig[sig_1 =="up"&sig_2=="up"&sig_0554=="up"]       <- "up"
  sig[sig_1 =="down"&sig_2=="down"&sig_0554=="down"] <- "down"
})
table(de_intersect_gene$sig)
save(de_intersect_gene,file = "Pro_TNBC/paper/data/results/section_5/de_intersect_gene_DESeq2.RData")
save(de_results_TNBC0554,file = "Pro_TNBC/paper/data/results/section_4/TNBC_0554/de_results_TNBC0554_deseq2.RData")
save(de_results_HDQP1_GSE173634,file = "Pro_TNBC/paper/data/results/section_4/GSE173634_HDQP1/de_results_HDQP1_GSE173634_deseq2.RData")
save(de_results_HDQP1_GSE202771,file = "Pro_TNBC/paper/data/results/section_4/GSE202771_HDQP1/de_results_HDQP1_GSE202771_deseq2.RData")

de_results_HDQP1_GSE173634   <- subset(de_results_HDQP1_GSE173634,sig=="up"|sig=="down")
de_results_HDQP1_GSE202771   <- subset(de_results_HDQP1_GSE202771,sig=="up"|sig=="down")
de_results_TNBC0554          <- subset(de_results_TNBC0554,sig=="up"|sig=="down")
de_intersect_gene            <- de_intersect_gene[!is.na(de_intersect_gene$sig),]
write.csv(de_results_HDQP1_GSE202771,file = "Pro_TNBC/paper/data/results/section_4/GSE202771_HDQP1/de_results_HDQP1_GSE202771.csv")
write.csv(de_results_HDQP1_GSE173634,file = "Pro_TNBC/paper/data/results/section_4/GSE173634_HDQP1/de_results_HDQP1_GSE173634.csv")
write.csv(de_results_TNBC0554,file = "Pro_TNBC/paper/data/results/section_4/TNBC_0554/de_results_TNBC0554.csv")
write.csv(de_intersect_gene,file = "Pro_TNBC/paper/data/results/section_5/de_intersect_gene.csv")



