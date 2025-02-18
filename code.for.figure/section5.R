#############################################################################################
####ELF3 knock-down induces partial Claudin-low phenotype in Basal-like cancer.####
#############################################################################################



####Figure 6: ELF3 knock-down induces partial Claudin-low phenotype in Basal-like cancer.####
####*Figure 6a: a heatmap for showing the average cnv in every cell type####
load("~/Pro_TNBC/output/data/scRNASeq/38sample/TNBC_BRCA1_0554/TNBC_BRCA1_0554_scRNA.RData")
run.final.infercnv_obj = readRDS("Pro_TNBC/output/data/scRNASeq/38sample/TNBC_BRCA1_0554/infercnv_basalandCL/run.final.infercnv_obj")
expr              <- run.final.infercnv_obj@expr.data
dim(expr)
expr[1:4,1:4]
gn                <- rownames(expr)
library(AnnoProbe)  
geneInfor=annoGene(gn,"ENSEMBL",'human')  
colnames(geneInfor)
geneInfor=geneInfor[with(geneInfor, order(chr, start)),c(3,4:6)]      
geneInfor=geneInfor[!duplicated(geneInfor[,1]),]
length(unique(geneInfor[,1]))
geneFile          <- geneInfor
geneFile$chr      <- factor(geneFile$chr,levels = paste0("chr", seq(1, 22)))
geneFile          <- arrange(geneFile,chr)
identical(rownames(expr),geneFile$ENSEMBL)
expr              <- expr[match(geneFile$ENSEMBL,rownames(expr)),]
geneFile$new_cluster    <- sapply(as.character(geneFile$chr), function(x) substr(x, 4, nchar(x)))
geneFile$new_cluster    <- factor(geneFile$new_cluster,levels = seq(1, 22))
library(ComplexHeatmap)
metadata          <- TNBC_BRCA1_0554_scRNA@meta.data[,11,drop=F]
identical(colnames(expr),rownames(metadata))
metadata          <- metadata[match(colnames(expr),rownames(metadata)),,drop=F]
library(dplyr)
library(tibble)
expr_df            <- as.data.frame(t(expr))
identical(rownames(expr_df),rownames(metadata))
expr_df            <- cbind(expr_df,metadata)
expr_df[,1:6281]   <- apply(expr_df[,1:6281],2,as.numeric)
averages           <- expr_df %>% group_by(cell_subtype) %>% summarise_all(mean)
averages           <- column_to_rownames(averages,"cell_subtype")
averages                 <- t(averages)
identical(geneFile$ENSEMBL,rownames(averages))
colnames(averages)
averages                  <- averages[,-c(7,8)]
colnames(averages)[1:10]  <- c("B-cell","Basal-like","CMP","Claudin-low","DC","Fibroblasts","Monocyte","Pericyte","T-cell","Normal-epithlial")
new_cluster              <- geneFile$new_cluster
top_labels               <- HeatmapAnnotation(
  cluster = anno_block( labels = levels(new_cluster), 
                        labels_gp = gpar(cex = 0.9, col = "black"))) 
pdf("Pro_TNBC/paper/plot/section_5/heatmap.of.cnv.in.TNBC.0554.pdf",width = 15,height = 10)
averages_heatmap  <- Heatmap(t(averages),
                             cluster_rows = T,
                             cluster_columns = F,
                             show_column_names = F,
                             show_row_names = T,
                             column_split = new_cluster,
                             show_heatmap_legend=T,
                             top_annotation = top_labels,
                             column_title = "Genomic Region",
                             column_title_side = c("bottom"),
                             column_title_gp = gpar(fontsize = 25),
                             heatmap_legend_param = list(
                               title = "Modified Expression",
                               title_position = "leftcenter-rot", 
                               at=c(0.8,1.2), 
                               legend_height = unit(3, "cm") 
                             ))
draw(averages_heatmap, heatmap_legend_side = "right") # 图例位置
dev.off()
range(averages)


save(averages,file = "Pro_TNBC/paper/data/results/section_5/averages.cnv.in.TNBC0554.RData")




####*Figure 6b: Differential expression between sh-ELF3 and control group. ####
HCC70_SHELF3_low_VS_high_deres_symbol_highCount  <- subset(HCC70_SHELF3_low_VS_high_deres_symbol_highCount,avg_log2FC <10)
HCC70_SHELF3_low_VS_high_deres_symbol_highCount             <- within(HCC70_SHELF3_low_VS_high_deres_symbol_highCount,{
  sig                                   <- NA
  sig[p_val_adj < 0.05 & avg_log2FC > 0.5]   <- "up"
  sig[p_val_adj < 0.05 & avg_log2FC < -0.5]  <- "down"
})
HCC70_SHELF3_low_VS_high_deres_symbol_highCount[is.na(HCC70_SHELF3_low_VS_high_deres_symbol_highCount$sig),8] <- "none"
#
de_genes                                  <- c("ELF3","CLDN4","CDH3","SDC4","VTCN1","MPZL1")
library(ggrepel)
fig5f                                     <- ggplot(HCC70_SHELF3_low_VS_high_deres_symbol_highCount,aes(avg_log2FC, -log10(p_val_adj),color=sig))+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black")+
  geom_vline(xintercept = c(-0.5,0.5), linetype = "dashed", color = "black")+
  geom_point(size=5)+
  scale_color_manual(values=c("#2f5688","#BBBBBB","#CC0000"))+
  scale_size_continuous(range = c(0,1))+
  theme(legend.title = element_blank() )+
  geom_label_repel(data = subset(HCC70_SHELF3_low_VS_high_deres_symbol_highCount, SYMBOL %in% de_genes),
                   max.overlaps = getOption("ggrepel.max.overlaps", default = 20),
                   aes(label = SYMBOL),
                   size = 8, 
                   color = 'black',fontface = "italic",
                   segment.color = "black",   
                   segment.linetype = "solid", 
                   segment.size = 0.5 ,
                   box.padding = 3,        
                   point.padding = 0.2 ) +
  xlab("log2fc")+
  ylab("-log10padjust")+ggplot.style

ggsave(fig5f,filename = "Pro_TNBC/paper/plot/section_5/Volcano.map.of.ELF3-low.pdf",width = 20,height = 15)





####*Figure 6c: heatmap of ELF3-low down genes in CCLE cells ####
load("~/Pro_TNBC/paper/data/results/section_5/BL_up_genes.RData")
load("~/Pro_TNBC/paper/data/results/section_5/cl_low_genes.RData")
##down genes in CCLE #
load("Pro_TNBC/data/CCLE/CCLE.RData")
brca_info                  <- read_csv("Pro_TNBC/output/data/CCLE/brca_info.csv")
brca_info                  <- brca_info[!is.na(brca_info$lineage_molecular_subtype),]
down_gene_id               <- subset(HCC70_SHELF3_low_VS_high_deres_symbol_highCount,sig=="down")
brca_log2rpkm              <- CCLE.log2.rpkm.matrix[,colnames(CCLE.log2.rpkm.matrix) %in% brca_info$CCLE_Name]
brca_log2rpkm_down         <- brca_log2rpkm[rownames(brca_log2rpkm) %in% down_gene_id$ENSEMBL,] %>% as.data.frame()
brca_log2rpkm_down$ENSEMBL <- rownames(brca_log2rpkm_down)
brca_log2rpkm_down         <- merge(down_gene_id[,c(2,9)],brca_log2rpkm_down,by="ENSEMBL")
brca_log2rpkm_down         <- column_to_rownames(brca_log2rpkm_down,var = "SYMBOL")
brca_log2rpkm_down         <- brca_log2rpkm_down[,-1] %>% as.matrix() %>% t()
brca_info                  <- arrange(brca_info,lineage_molecular_subtype)
brca_log2rpkm_down                   <- brca_log2rpkm_down[match(brca_info$CCLE_Name,rownames(brca_log2rpkm_down)),]
brca_log2rpkm_down                   <- t(brca_log2rpkm_down)
brca_info                 <- brca_info[-c(49,50),]
brca_info$subtype        <- as.factor(brca_info$lineage_molecular_subtype)   
brca_info                <- arrange(brca_info,subtype)
brca_log2rpkm_down       <- brca_log2rpkm_down[,-c(49,50)]
brca_log2rpkm_down                   <- brca_log2rpkm_down[,match(brca_info$CCLE_Name,colnames(brca_log2rpkm_down))]
identical(colnames(brca_log2rpkm_down),brca_info$CCLE_Name)
cl_low_genes_symbol       <- bitr(cl_low_genes,fromType = "ENSEMBL",toType = "SYMBOL",OrgDb = "org.Hs.eg.db")       
brca_log2rpkm_down_cldown <- brca_log2rpkm_down[rownames(brca_log2rpkm_down)%in% cl_low_genes_symbol$SYMBOL,]
brca_log2rpkm_down_cldown_scale <- apply(brca_log2rpkm_down_cldown,1,scale)
rownames(brca_log2rpkm_down_cldown_scale) <- colnames(brca_log2rpkm_down_cldown)

BL_up_genes_symbol        <- bitr(BL_up_genes,fromType = "ENSEMBL",toType = "SYMBOL",OrgDb = "org.Hs.eg.db")       
brca_log2rpkm_down_BLup   <- brca_log2rpkm_down[rownames(brca_log2rpkm_down)%in% BL_up_genes_symbol$SYMBOL,]
brca_log2rpkm_down_BLup_scale <- apply(brca_log2rpkm_down_BLup,1,scale)
rownames(brca_log2rpkm_down_BLup_scale) <- colnames(brca_log2rpkm_down_BLup)

lest_genes <- setdiff(rownames(brca_log2rpkm_down),intersect(cl_low_genes_symbol$SYMBOL,rownames(brca_log2rpkm_down)))
lest_genes  <- setdiff(lest_genes,intersect(BL_up_genes_symbol$SYMBOL,rownames(brca_log2rpkm_down)))
brca_log2rpkm_down_lest   <- brca_log2rpkm_down[rownames(brca_log2rpkm_down)%in%lest_genes,]
brca_log2rpkm_down_lest_scale <- apply(brca_log2rpkm_down_lest,1,scale)
rownames(brca_log2rpkm_down_lest_scale) <- colnames(brca_log2rpkm_down_lest)

col.annotation <- HeatmapAnnotation(
  subtype = brca_info$subtype,
  col = list(subtype=c("basal_A"="Red","basal_B"= "Blue","HER2_amp"= "Purple","luminal"= "#FFCC00")))

pdf(file = "Pro_TNBC/paper/plot/section_5/heatmap.of.ELF3-low.down.genes.pdf",width = 20,height = 20)
#set.seed(123)
ht1 <- Heatmap(
  t(brca_log2rpkm_down_cldown_scale),
  name = "Expression",
  col =  colorRamp2(c(-5, 0.1, 5), c("navy","white","firebrick3")),
  cluster_rows = T,
  cluster_columns = F,
  show_column_names = F,
  show_row_names = T,
  show_row_dend = T,
  top_annotation = col.annotation,
  width = ncol(brca_log2rpkm_down_cldown)*unit(7, "mm"), 
  height = nrow(brca_log2rpkm_down_cldown)*unit(5, "mm"),
  rect_gp = gpar(col = "black", lwd = 0.1)
)
ht2 <- Heatmap(
  t(brca_log2rpkm_down_BLup_scale),
  col =  colorRamp2(c(-5, 0.1, 5), c("navy","white","firebrick3")),
  cluster_rows = T,
  cluster_columns = F,
  show_column_names = F,
  show_row_names = T,
  show_row_dend = T,
  width = ncol(brca_log2rpkm_down_BLup)*unit(7, "mm"), 
  height = nrow(brca_log2rpkm_down_BLup)*unit(5, "mm"),
  rect_gp = gpar(col = "black", lwd = 0.1)
)
ht3 <- Heatmap(
  t(brca_log2rpkm_down_lest_scale),
  col =  colorRamp2(c(-5, 0.1, 5), c("navy","white","firebrick3")),
  cluster_rows = T,
  cluster_columns = F,
  show_column_names = F,
  show_row_names = T,
  show_row_dend = T,
  width = ncol(brca_log2rpkm_down_lest)*unit(7, "mm"), 
  height = nrow(brca_log2rpkm_down_lest)*unit(5, "mm"),
  rect_gp = gpar(col = "black", lwd = 0.1)
)
ht_list = ht1 %v% ht2 %v% ht3
draw(ht_list)
dev.off()


####supplementary Figure 7: ELF3 expression of Basal-like and Claudin-low cancer cells in HDQP1 cell line (from GSE173634 and GSE202771 datasets) and TNBC patient 0554. ####
load("Pro_TNBC/paper/data/results/section_4/GSE202771_HDQP1_scRNA.RData")
GSE202771.HDQP1.cpm              <- as.matrix(GSE202771_HDQP1_scRNA@assays$RNA@data) %>% t() %>% as.data.frame()
GSE202771.HDQP1.cpm.ELF3         <- GSE202771.HDQP1.cpm[,"ELF3",drop=F]
identical(rownames(GSE202771.HDQP1.cpm.ELF3),rownames(GSE202771_HDQP1_scRNA@meta.data))
GSE202771.HDQP1.cpm.ELF3$subtype <- GSE202771_HDQP1_scRNA$subtype
GSE202771.HDQP1.cpm.ELF3$sample  <- rep("GSE202771.HDQP1",nrow(GSE202771.HDQP1.cpm.ELF3))

load("Pro_TNBC/paper/data/results/section_4/GSE173634_HDQP1_scRNA.RData")
GSE173634.HDQP1.cpm                     <- as.matrix(GSE173634.HDQP1.scRNA@assays$RNA@data) %>% t() %>% as.data.frame()
GSE173634.HDQP1.cpm.ELF3         <- GSE173634.HDQP1.cpm[,"ENSG00000163435",drop=F]
identical(rownames(GSE173634.HDQP1.cpm.ELF3),rownames(GSE173634.HDQP1.scRNA@meta.data))
GSE173634.HDQP1.cpm.ELF3$subtype <- GSE173634.HDQP1.scRNA$subtype
GSE173634.HDQP1.cpm.ELF3$sample  <- rep("GSE173634.HDQP1",nrow(GSE173634.HDQP1.cpm.ELF3))

load("Pro_TNBC/paper/data/results/section_4/TNBC_BRCA1_0554_tumor_scRNA(tsne_by_cor).RData")
TNBC.0554.cpm                     <- as.matrix(TNBC_BRCA1_0554_tumor_scRNA@assays$RNA@data) %>% t() %>% as.data.frame()
TNBC.0554.cpm.ELF3         <- TNBC.0554.cpm[,"ENSG00000163435",drop=F]
identical(rownames(TNBC.0554.cpm.ELF3),rownames(TNBC_BRCA1_0554_tumor_scRNA@meta.data))
TNBC.0554.cpm.ELF3$subtype <- TNBC_BRCA1_0554_tumor_scRNA$subtype
TNBC.0554.cpm.ELF3$sample  <- rep("TNBC.0554",nrow(TNBC.0554.cpm.ELF3))

colnames(GSE173634.HDQP1.cpm.ELF3)[1] <- "ELF3"
colnames(TNBC.0554.cpm.ELF3)[1] <- "ELF3"
cpm.ELF3 <- bind_rows(GSE173634.HDQP1.cpm.ELF3,GSE202771.HDQP1.cpm.ELF3,TNBC.0554.cpm.ELF3)
cpm.ELF3  <- subset(cpm.ELF3,subtype=="Basal"|subtype=="Claudin_low")

figS7     <- ggboxplot(cpm.ELF3,x = "sample",y = "ELF3",color = "subtype",add = "point",size= 1,add.params=list(size=3))+ 
  ggplot.style+
  scale_color_manual(values =c(Claudin_low="Blue",Basal="Red"))+
  theme(panel.grid.major=element_line(colour=NA),
        panel.grid.minor = element_blank())+xlab("")+ylab("ELF3")
ggsave(figS7,filename = "Pro_TNBC/paper/plot/section_5/the.expression.of.ELF3.pdf",width = 20,height = 15)

cpm.ELF3.GSE173634     <- subset(cpm.ELF3,sample=="GSE173634.HDQP1")
wilcox.test(ELF3~subtype,cpm.ELF3.GSE173634)#p-value = 2.886e-08
cpm.ELF3.GSE202771     <- subset(cpm.ELF3,sample=="GSE202771.HDQP1")
wilcox.test(ELF3~subtype,cpm.ELF3.GSE202771)# p-value < 2.2e-16
cpm.ELF3.TNBC          <- subset(cpm.ELF3,sample=="TNBC.0554")
wilcox.test(ELF3~subtype,cpm.ELF3.TNBC)# p-value < 2.2e-16


####Supplementary Figure 8a :cor with centroid####
HCC70_SHELF3_high_low_scRNA <- readRDS("Pro_TNBC/paper/data/results/section_5/HCC70_SHELF3_high_low_scRNA.rds")
ELF3_high_low_df            <- data.frame(cors=c(HCC70_SHELF3_high_low_scRNA$Basal,
                                                 HCC70_SHELF3_high_low_scRNA$Claudin_low),
                                          group=c(HCC70_SHELF3_high_low_scRNA$group,
                                                  HCC70_SHELF3_high_low_scRNA$group),
                                          centroid=c(rep("Basal-like",200),rep("Claudin-low",200)))

ELF3_high_low_df        <- subset(ELF3_high_low_df,centroid=="Basal-like")
ELF3_high_low_df$group  <- factor(ELF3_high_low_df$group,levels = c("ELF3_high","ELF3_low"))
library(ggpubr)
fig_5f     <- ggboxplot(ELF3_high_low_df,x = "group",y = "cors",add = "point",size= 1,add.params=list(size=3))+ 
  ggplot.style+stat_compare_means(method = "wilcox.test")+
  theme(panel.grid.major=element_line(colour=NA),
        panel.grid.minor = element_blank())+xlab("")+ylab("Correlation coefficient")
ggsave(fig_5f,filename = "Pro_TNBC/paper/plot/section_5/fig5d_boxplot.of.correlation.in.two.group.pdf",width = 20,height = 15)

####Supplementary Figure 8b: GSEA for cl-low genes in de results of ELF3-low ####
library(clusterProfiler)
library(enrichplot)
alldiff               <- subset(HCC70_SHELF3_low_VS_high_deres_symbol_highCount,sig=="up"|sig=="down")
alldiff               <- alldiff[order(alldiff$avg_log2FC,decreasing = T),]

genelist              <- alldiff$avg_log2FC
names(genelist)       <- alldiff$ENSEMBL
cl_low_genes_df       <- data.frame(term=rep("Claudin-low low genes",length(cl_low_genes)),cl_low_genes=cl_low_genes)
cl.low.gsea.re        <- GSEA(genelist,  
                              TERM2GENE = cl_low_genes_df,  
                              pvalueCutoff = 1 , 
                              pAdjustMethod='BH')  
pdf("Pro_TNBC/paper/plot/section_5/gsea.cl.low.pdf",width =12,height =12)
gseaplot2(cl.low.gsea.re,
          cl.low.gsea.re$ID,
          color = "red",
          base_size = 20,
          rel_heights = c(1.5, 0.5, 1),
          subplots = 1:3, 
          ES_geom = "line",
          pvalue_table = T) 
dev.off()
save(cl.low.gsea.re,file = "Pro_TNBC/paper/data/results/section_4/cl.low.gsea.re.RData")





de_results_HDQP1_GSE173634  <- within(de_results_HDQP1_GSE173634,{
  sig <- NA
  sig[p_val_adj < 0.05& avg_log2FC >0.5] <- "up"
  sig[p_val_adj < 0.05& avg_log2FC < -0.5] <- "down"
})
de_GSE173634 <- de_results_HDQP1_GSE173634[!is.na(de_results_HDQP1_GSE173634$sig),]

de_results_HDQP1_GSE202771  <- within(de_results_HDQP1_GSE202771,{
  sig <- NA
  sig[p_val_adj < 0.05& avg_log2FC >0.5] <- "up"
  sig[p_val_adj < 0.05& avg_log2FC < -0.5] <- "down"
})
de_GSE202771 <- de_results_HDQP1_GSE202771[!is.na(de_results_HDQP1_GSE202771$sig),]
de_results_TNBC0554  <- within(de_results_TNBC0554,{
  sig <- NA
  sig[p_val_adj < 0.05& avg_log2FC >0.5] <- "up"
  sig[p_val_adj < 0.05& avg_log2FC < -0.5] <- "down"
})
de_TNBC0554  <- de_results_TNBC0554[!is.na(de_results_TNBC0554$sig),]

de_intersect <- merge(de_GSE173634,de_GSE202771,by="SYMBOL")
de_intersect <- merge(de_intersect,de_TNBC0554,by="SYMBOL")
identical(de_intersect$sig,de_intersect$sig.y)
identical(de_intersect$sig.x,de_intersect$sig.z)
de_intersect_genes  <- de_intersect[,c(1,2,4,7,8,10,13,14,18,21,22)]
colnames(de_intersect_genes)[1:11] <- c("SYMBOL","ENSEMBL","log2fc_HDQP1_GSE173634","padj_HDQP1_GSE173634","sig.x","log2fc_HDQP1_GSE202771","padj_HDQP1_GSE202771","sig.y","log2fc_0554","padj_0554","sig")
de_intersect_genes <- de_intersect_genes %>%
  filter(sig.x == sig.y) 
de_intersect_genes <- de_intersect_genes %>%
  filter(sig.x == sig) 
de_intersect_genes  <- de_intersect_genes[,-c(5,8)]
de_intersect_genes <- de_intersect_genes[-6,]
write.csv(de_intersect_genes,file = "Pro_TNBC/paper/data/results/section_5/de_intersect_genes.csv")
