################################################################################
####Loss of ELF3 expression drives Basal to Claudin-low trans-differentiation
################################################################################


################################################################################
####Find the differentiated genes between Basal and Claudin-low cells
################################################################################

####*DE analysis of basal and claudin-low cells in TNBC 0554####
library(Seurat)
Idents(TNBC_BRCA1_0554_tumor_scRNA)      <- "subtype"
de_results_TNBC0554                      <- FindMarkers(TNBC_BRCA1_0554_tumor_scRNA, slot="counts",ident.1 = "Claudin_low", ident.2 = "Basal",test.use = "DESeq2")
de_results_TNBC0554                      <- within(de_results_TNBC0554,{
  sig                                    <- NA
  sig[p_val_adj < 0.05 &avg_log2FC >1]   <- "up"
  sig[p_val_adj < 0.05 &avg_log2FC < -1] <- "down"
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
  sig[p_val_adj < 0.05 &avg_log2FC >1]   <- "up"
  sig[p_val_adj < 0.05 &avg_log2FC < -1] <- "down"
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
  sig[p_val_adj < 0.05 &avg_log2FC >1]   <- "up"
  sig[p_val_adj < 0.05 &avg_log2FC < -1] <- "down"
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
de_HDQP1_claudinlow_basal_intersect                 <- de_HDQP1_claudinlow_basal_intersect[,c(1,3,6,7,8,10,13,14)]
colnames(de_HDQP1_claudinlow_basal_intersect)[2:8]  <- c("log2FoldChange_HDQP1_1","padj_HDQP1_1","sig_1","SYMBOL","log2FoldChange_HDQP1_2","padj_HDQP1_2","sig_2")
de_intersect_gene                     <- merge(de_HDQP1_claudinlow_basal_intersect,de_results_TNBC0554,by="ENSEMBL")
de_intersect_gene                     <- de_intersect_gene[,c(1:8,10,13,14)]
colnames(de_intersect_gene)[c(9:11)]  <- c("log2FoldChange_0554","padj_0554","sig_0554")
de_intersect_gene                     <- within(de_intersect_gene,{
  sig                                                <- NA
  sig[sig_1 =="up"&sig_2=="up"&sig_0554=="up"]       <- "up"
  sig[sig_1 =="down"&sig_2=="down"&sig_0554=="down"] <- "down"
})
table(de_intersect_gene$sig)
save(de_intersect_gene,file = "Pro_TNBC/paper/data/results/section_5/de_intersect_gene_DESeq2.RData")
save(de_results_TNBC0554,file = "Pro_TNBC/paper/data/results/section_5/de_results_TNBC0554_deseq2.RData")
save(de_results_HDQP1_GSE173634,file = "Pro_TNBC/paper/data/results/section_5/de_results_HDQP1_GSE173634_deseq2.RData")
save(de_results_HDQP1_GSE202771,file = "Pro_TNBC/paper/data/results/section_5/de_results_HDQP1_GSE202771_deseq2.RData")


####*Figure 5a: the valcano plot of different expression genes in TNBC 0554####
load("Pro_TNBC/paper/data/results/section_5/de_results_TNBC0554_deseq2.RData")
load("~/Pro_TNBC/paper/data/results/section_5/de_intersect_gene_DESeq2.RData")
de_gene                                   <- de_intersect_gene[de_intersect_gene$sig=="up"|de_intersect_gene$sig=="down",]$SYMBOL
de_gene                                   <- de_gene[!is.na(de_gene)]
library(ggrepel)
fig5a                                     <- ggplot(de_results_TNBC0554,aes(avg_log2FC, -log10(p_val_adj),color=sig))+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black")+
  geom_vline(xintercept = c(-1,1), linetype = "dashed", color = "black")+
  geom_point()+
  scale_color_manual(values=c("#CC0000","#BBBBBB","#2f5688"))+
  scale_size_continuous(range = c(0,1))+
  theme(legend.title = element_blank() )+
  geom_text_repel(data = filter(de_results_TNBC0554, SYMBOL %in% de_gene),
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 20),
                  aes(label = SYMBOL),
                  size = 8, 
                  color = 'black') +
  xlab("Log2FC")+
  ylab("-Log10(FDR q-value)")+ggplot.style

ggsave(fig5a,filename = "Pro_TNBC/paper/plot/section_5/Volcano.map.of.0554.pdf",width = 20,height = 15)


####*Figure 5b: the valcano plot of different expression genes in HDQP1 cell line (GSE173634)####
load("Pro_TNBC/paper/data/results/section_5/de_results_HDQP1_GSE173634_deseq2.RData")
fig5b                                            <- ggplot(de_results_HDQP1_GSE173634,aes(avg_log2FC, -log10(p_val_adj),color=sig))+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black")+
  geom_vline(xintercept = c(-1,1), linetype = "dashed", color = "black")+
  geom_point()+
  scale_color_manual(values=c("#CC0000","#BBBBBB","#2f5688"))+
  scale_size_continuous(range = c(0,1))+
  theme(legend.title = element_blank() )+
  geom_text_repel(data = filter(de_results_HDQP1_GSE173634, SYMBOL %in% de_gene),
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 20),
                  aes(label = SYMBOL),
                  size = 8, 
                  color = 'black') +
  xlab("Log2FC")+
  ylab("-Log10(FDR q-value)")+ggplot.style

ggsave(fig5b,filename = "Pro_TNBC/paper/plot/section_5/Volcano.map.of.GSE173634.HDQP1.pdf",width = 20,height = 15)


####*Figure 5c: the valcano plot of different expression genes in HDQP1 cell line (GSE202771)####
load("Pro_TNBC/paper/data/results/section_5/de_results_HDQP1_GSE202771_deseq2.RData")
fig5c        <- ggplot(de_results_HDQP1_GSE202771,aes(avg_log2FC, -log10(p_val_adj),color=sig))+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black")+
  geom_vline(xintercept = c(-1,1), linetype = "dashed", color = "black")+
  geom_point()+
  scale_color_manual(values=c("#CC0000","#BBBBBB","#2f5688"))+
  scale_size_continuous(range = c(0,1))+
  theme(legend.title = element_blank() )+
  geom_text_repel(data = filter(de_results_HDQP1_GSE202771, SYMBOL %in% de_gene),
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 20),
                  aes(label = SYMBOL),
                  size = 8, 
                  color = 'black') +
  xlab("Log2FC")+
  ylab("-Log10(FDR q-value)")+ggplot.style

ggsave(fig5c,filename = "Pro_TNBC/paper/plot/section_5/Volcano.map.of.GSE202271.HDQP1.pdf",width = 20,height = 15)
