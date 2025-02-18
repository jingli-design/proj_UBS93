####################################################################################
####Identification of Claudin-low cancer cells in Basal-like cell line and patient
####################################################################################

library(ggplot2)
library(ggrepel)
source("Pro_TNBC/paper/code/ggplot.style.R")
##################################################################################################################
####Figure 4: Discovery of “Basal-like/Claudin-low mixed” subtype in HDQP1 cell line (from GSE173634 dataset).####
##################################################################################################################

####*Figure 4a: the t-SNE plot for HDQP1 data####
load("Pro_TNBC/paper/data/results/section_4/GSE173634_HDQP1_scRNA.RData")
df1                              <- as.data.frame(GSE173634.HDQP1.scRNA@reductions$tsne@cell.embeddings)
df1$cell.id                      <- rownames(df1)
subtype                          <- as.data.frame(GSE173634.HDQP1.scRNA@meta.data)
subtype$cell.id                  <- rownames(subtype)
df1                              <- merge(df1,subtype,by="cell.id")
fig_4a                           <- ggplot(data = df1,aes(x=tSNE_1,y=tSNE_2,col=subtype))+
  geom_point(size=5)+
  ggplot.style+theme(plot.title = element_text(hjust=0.5,size = 44))+scale_color_manual(values = c(Basal="Red",Claudin_low="Blue"))
ggsave(fig_4a,filename = "Pro_TNBC/paper/plot/section_4/tsne.plot.of.HDQP1.by.using.all.genes.pdf",width=23,height=15)


####*Figure 4b: the valcano plot of different expression genes in HDQP1 cell line (GSE173634)####
load("Pro_TNBC/paper/data/results/section_5/de_results_HDQP1_GSE173634_deseq2.RData")
marker_genes                               <- c("CLDN3","CLDN4","CLDN7","EPCAM","VIM","CD24","CD44")
de_results_HDQP1_GSE173634                 <- within(de_results_HDQP1_GSE173634,{
  sig                                      <- NA
  sig[p_val_adj < 0.05 &avg_log2FC >0.5]   <- "up"
  sig[p_val_adj < 0.05 &avg_log2FC < -0.5] <- "down"
})

fig5b                                            <- ggplot(de_results_HDQP1_GSE173634,aes(avg_log2FC, -log10(p_val_adj),color=sig))+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black")+
  geom_vline(xintercept = c(-0.5,0.5), linetype = "dashed", color = "black")+
  geom_point(size=5)+
  scale_color_manual(values=c("#2f5688","#BBBBBB","#CC0000"))+
  scale_size_continuous(range = c(0,1))+
  theme(legend.title = element_blank() )+
  geom_label_repel(data = subset(de_results_HDQP1_GSE173634, SYMBOL %in% marker_genes),
                   max.overlaps = getOption("ggrepel.max.overlaps", default = 20),
                   aes(label = SYMBOL),
                   size = 8, 
                   color = 'black',
                   fontface="italic",
                   segment.color = "black",   
                   segment.linetype = "solid", 
                   segment.size = 0.5 ,
                   box.padding = 3,        
                   point.padding = 0.2 ) +
  xlab("log2fc")+
  ylab("-log10p-adjust")+ggplot.style

ggsave(fig5b,filename = "Pro_TNBC/paper/plot/section_4/HDQP1_GSE173634/Volcano.map.of.GSE173634.HDQP1.pdf",width = 20,height = 15)

####*Figure 4c: GSEA plot for EMT signature in DE genes in GSE173634 HDQP1####
library(clusterProfiler)
library(enrichplot)
library(msigdbr)
alldiff               <- subset(de_results_HDQP1_GSE173634,sig=="up"|sig=="down")
alldiff               <- alldiff[order(alldiff$avg_log2FC,decreasing = T),]
genelist              <- alldiff$avg_log2FC
names(genelist)       <- alldiff$SYMBOL
hallmark              <- msigdbr(species = "Homo sapiens", category = "H")
hallmark              <- hallmark[,c(3,4)]
colnames(hallmark)    <- c('term','gene')
hallmark              <- subset(hallmark,term=="HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")
hallmark$term         <- "EMT"
GSE173634.HDQP1.gsea.re              <- GSEA(genelist,  
                TERM2GENE = hallmark,  
                pvalueCutoff = 1 , 
                pAdjustMethod='BH')  
pdf("Pro_TNBC/paper/plot/section_4/HDQP1_GSE173634/gsea.EMT.pdf",width =12,height =12)
p1 <- gseaplot2(GSE173634.HDQP1.gsea.re,
          GSE173634.HDQP1.gsea.re$ID,
                    color = "red",
                    base_size = 20,
                    rel_heights = c(1.5, 0.5, 1),
                    subplots =1,
                    ES_geom = "line",
                    pvalue_table = F)
p2 <- gseaplot2(GSE173634.HDQP1.gsea.re,
                GSE173634.HDQP1.gsea.re$ID,
                color = "red",
                base_size = 20,
                rel_heights = c(1.5, 0.5, 1),
                subplots =2,
                ES_geom = "line",
                pvalue_table = F)
p3 <- plot_grid(p1,p2,ncol = 1)
ggsave(p3,filename ="Pro_TNBC/paper/plot/section_4/HDQP1_GSE173634/gsea.EMT.pdf",width =12,height =12)

save(GSE173634.HDQP1.gsea.re,file = "Pro_TNBC/paper/data/results/section_4/GSE173634.HDQP1.gsea.re.RDa")


####*Figure 4d: CD44/CD24 ratio BL and CL cells in GSE173634 HDQP1####
library(ggpubr)
load("Pro_TNBC/paper/data/results/section_4/GSE173634_HDQP1_scRNA.RData")
GSE173634.HDQP1.cpm                     <- as.matrix(GSE173634.HDQP1.scRNA@assays$RNA@data) %>% t() %>% as.data.frame()
GSE173634.HDQP1.CD44.CD24.ratio         <- GSE173634.HDQP1.cpm[,c("ENSG00000026508","ENSG00000272398")] 
colnames(GSE173634.HDQP1.CD44.CD24.ratio)[1:2]   <- c("CD44","CD24")
GSE173634.HDQP1.CD44.CD24.ratio$subtype <- GSE173634.HDQP1.scRNA$subtype
GSE173634.HDQP1.CD44.CD24.ratio$ratio   <- (GSE173634.HDQP1.CD44.CD24.ratio$CD44+0.1)/(GSE173634.HDQP1.CD44.CD24.ratio$CD24+0.1)
save(GSE173634.HDQP1.CD44.CD24.ratio,file = "Pro_TNBC/paper/data/results/section_4/GSE173634.HDQP1.CD44.CD24.ratio.RData")

fig_4d     <- ggboxplot(GSE173634.HDQP1.CD44.CD24.ratio,x = "subtype",y = "ratio",add = "point",size= 1,add.params=list(size=3))+ 
  ggplot.style+
  theme(panel.grid.major=element_line(colour=NA),
        panel.grid.minor = element_blank())+
  xlab("")+ylab("CD44/CD24 ratio")+stat_compare_means(method = "wilcox.test")+
  scale_x_discrete("", labels = c("Basal" = paste("Basal-like",sep = "\n","(n=58)"),"Claudin_low" = paste("Claudin-low",sep = "\n","(n=48)")))
ggsave(fig_4d,filename = "Pro_TNBC/paper/plot/section_4/HDQP1_GSE173634/fig4d_boxplot.of.ratio.in.GSE173634.HDQP1.pdf",width = 15,height = 15)

####*Figure 4e: correlation coefficents between claudin-low cells  and CCLE ####
load("Pro_TNBC/paper/data/results/section_4/cor.with.CCLE.RData")
CL.df                 <- CL_rs$correlation.matrix %>% t()%>% as.data.frame()
CL.df$CCLE_Name       <- rownames(CL.df)
brca_info             <- read_csv("Pro_TNBC/output/data/CCLE/brca_info.csv")
brca_info             <- brca_info[!is.na(brca_info$lineage_molecular_subtype),]
CL.df                 <- merge(CL.df,brca_info[,c(2,5,22)],by="CCLE_Name")
CL.df$subtype         <- ifelse(CL.df$lineage_molecular_subtype=="basal_B","Claudin-low",
                                ifelse(CL.df$lineage_molecular_subtype=="basal_A","Basal-like","Others"))
CL.df[is.na(CL.df$subtype),7]  <- "Others"
CL.df$HDQP1_GSE173634_CL_rank               <- rank(CL.df$HDQP1_GSE173634_CL) 
CL.df                 <- CL.df[order(CL.df$HDQP1_GSE173634_CL_rank,decreasing = F),]
fig1e <- ggplot(CL.df, aes(x = HDQP1_GSE173634_CL_rank, y = HDQP1_GSE173634_CL,color=subtype)) +
  geom_point(size = 6) +  
  ylim(-0.1,0.8)+
  labs(x = "Rank", y = "Correlation coefficients") +
  scale_color_manual(values = c("Basal-like" = "red","Claudin-low"="blue" ,"others" = "black"))  +  
  theme_bw()+
  geom_point(data = subset(CL.df, lineage_molecular_subtype == "basal_B"), color = "blue", size = 6) +
  geom_point(data = subset(CL.df, lineage_molecular_subtype == "basal_A"), color = "red", size = 6)+ggplot.style +
  theme(panel.grid = element_blank())+
  geom_text_repel(data = CL.df[46:50,], 
                  aes(label = stripped_cell_line_name), 
                  nudge_y = 0.08,  # 标签在纵坐标上的偏移量
                  color = "black",
                  size = 6,
                  arrow = arrow(length = unit(0.02, "npc")))
ggsave(fig1e,filename = "Pro_TNBC/paper/plot/section_4/HDQP1_GSE173634/cor.between.CCLE.and.GSE173634.HDQP1.CL.pdf",width = 20,height = 15)

####*Figure 4f: correlation coefficents between basal-like cells  and CCLE ####
load("Pro_TNBC/paper/data/results/section_4/cor.with.CCLE.RData")
BL.df                 <- BL_rs$correlation.matrix %>% t()%>% as.data.frame()
BL.df$CCLE_Name       <- rownames(BL.df)
brca_info             <- read_csv("Pro_TNBC/output/data/CCLE/brca_info.csv")
brca_info             <- brca_info[!is.na(brca_info$lineage_molecular_subtype),]
BL.df                 <- merge(BL.df,brca_info[,c(2,5,22)],by="CCLE_Name")
BL.df$subtype         <- ifelse(BL.df$lineage_molecular_subtype=="basal_B","Claudin-low",
                                ifelse(BL.df$lineage_molecular_subtype=="basal_A","Basal-like","Others"))
BL.df[is.na(BL.df$subtype),7]  <- "Others"
BL.df$HDQP1_GSE173634_BL_rank               <- rank(BL.df$HDQP1_GSE173634_BL) 
BL.df                 <- BL.df[order(BL.df$HDQP1_GSE173634_BL_rank,decreasing = F),]
fig4f                 <- ggplot(BL.df, aes(x = HDQP1_GSE173634_BL_rank, y = HDQP1_GSE173634_BL,color=subtype)) +
  geom_point(size = 6) + 
  ylim(-0.1,0.8)+
  labs(x = "Rank", y = "Correlation coefficients") +
  scale_color_manual(values = c("Basal-like" = "red","Claudin-low"="blue" ,"others" = "black"))  +  
  theme_bw()+
  geom_point(data = subset(BL.df, lineage_molecular_subtype == "basal_B"), color = "blue", size = 6) +
  geom_point(data = subset(BL.df, lineage_molecular_subtype == "basal_A"), color = "red", size = 6)+ggplot.style +
  theme(panel.grid = element_blank())+
  geom_text_repel(data = BL.df[46:50,], 
                  aes(label = stripped_cell_line_name), 
                  nudge_y = 0.08,  # 标签在纵坐标上的偏移量
                  color = "black",
                  size = 6,
                  arrow = arrow(length = unit(0.02, "npc")))
ggsave(fig4f,filename = "Pro_TNBC/paper/plot/section_4/HDQP1_GSE173634/cor.between.CCLE.and.GSE173634.HDQP1.BL.pdf",width = 20,height = 15)
####*Figure 4g: correlation coefficients between Claudin-low centroid and Claudin-low cells in TNBC ####
load("Pro_TNBC/paper/data/results/section_4/TNBC/TNBC_CL_cor.RData")
TNBC_cells                       <- table(TNBC_CL_cor$patient.id)
TNBC_cells                       <- TNBC_cells[TNBC_cells>20]
TNBC_CL_cor                      <- TNBC_CL_cor[TNBC_CL_cor$patient.id %in% names(TNBC_cells),]
fig_1b                           <- ggplot(TNBC_CL_cor,aes(x=reorder(patient.id,Claudin_low,FUN = median),y=Claudin_low))+
  geom_boxplot() + ggplot.style+
  theme(axis.text = element_text(size=12, face="bold",angle = 90),
        axis.title = element_text(size=45, face="bold"))+
  geom_point()+
  theme(panel.grid.major=element_line(colour=NA),
        panel.grid.minor = element_blank())+xlab("Patient id")+ylab("Correlation coefficients")+
  scale_x_discrete("Patient ID", labels = c(
    "BIOKEY_1" = paste("BIOKEY_1", sep = "\n", "(n=55)"),
    "BIOKEY_10" = paste("BIOKEY_10", sep = "\n", "(n=456)"),
    "BIOKEY_11" = paste("BIOKEY_11", sep = "\n", "(n=29)"),
    "BIOKEY_14" = paste("BIOKEY_14", sep = "\n", "(n=25)"),
    "BIOKEY_15" = paste("BIOKEY_15", sep = "\n", "(n=29)"),
    "BIOKEY_16" = paste("BIOKEY_16", sep = "\n", "(n=49)"),
    "BIOKEY_19" = paste("BIOKEY_19", sep = "\n", "(n=38)"),
    "BIOKEY_2" = paste("BIOKEY_2", sep = "\n", "(n=44)"),
    "BIOKEY_25" = paste("BIOKEY_25", sep = "\n", "(n=2)"),
    "BIOKEY_26" = paste("BIOKEY_26", sep = "\n", "(n=465)"),
    "BIOKEY_31" = paste("BIOKEY_31", sep = "\n", "(n=35)"),
    "BIOKEY_33" = paste("BIOKEY_33", sep = "\n", "(n=13)"),
    "BIOKEY_35" = paste("BIOKEY_35", sep = "\n", "(n=159)"),
    "BIOKEY_36" = paste("BIOKEY_36", sep = "\n", "(n=4)"),
    "BIOKEY_39" = paste("BIOKEY_39", sep = "\n", "(n=5)"),
    "BIOKEY_41" = paste("BIOKEY_41", sep = "\n", "(n=127)"),
    "BIOKEY_8" = paste("BIOKEY_8", sep = "\n", "(n=2)"),
    "BIOKEY_9" = paste("BIOKEY_9", sep = "\n", "(n=305)"),
    "CID4495" = paste("CID4495", sep = "\n", "(n=93)"),
    "CID44971" = paste("CID44971", sep = "\n", "(n=11)"),
    "CID44991" = paste("CID44991", sep = "\n", "(n=7)"),
    "CID4513" = paste("CID4513", sep = "\n", "(n=979)"),
    "CID4523" = paste("CID4523", sep = "\n", "(n=8)"),
    "PT039" = paste("PT039", sep = "\n", "(n=2)"),
    "PT084" = paste("PT084", sep = "\n", "(n=3)"),
    "PT089" = paste("PT089", sep = "\n", "(n=10)"),
    "PT126" = paste("PT126", sep = "\n", "(n=15)"),
    "sc5rJUQ033" = paste("sc5rJUQ033", sep = "\n", "(n=40)"),
    "sc5rJUQ039" = paste("sc5rJUQ039", sep = "\n", "(n=1)"),
    "sc5rJUQ042" = paste("sc5rJUQ042", sep = "\n", "(n=6)"),
    "sc5rJUQ045" = paste("sc5rJUQ045", sep = "\n", "(n=347)"),
    "sc5rJUQ053" = paste("sc5rJUQ053", sep = "\n", "(n=36)"),
    "SU4" = paste("SU4", sep = "\n", "(n=21)"),
    "SU58" = paste("SU58", sep = "\n", "(n=1)"),
    "TNBC_0126" = paste("TNBC_0126", sep = "\n", "(n=12)"),
    "TNBC_0131" = paste("TNBC_0131", sep = "\n", "(n=25)"),
    "TNBC_0135" = paste("TNBC_0135", sep = "\n", "(n=3)"),
    "TNBC_0554" = paste("TNBC_0554", sep = "\n", "(n=215)")
  ))
ggsave(fig_1b,filename = "Pro_TNBC/paper/plot/section_4/Box.plot.of.cor.cl.cells.and.centroid.pdf",height = 15,width = 20)





############################################################################################
####Figure 5:  Discovery of “Basal-like/Claudin-low mixed” subtype in TNBC patient 0554.####
############################################################################################

####*Figure 5a: the t-SNE plot for TNBC 0554 data####
load("Pro_TNBC/paper/data/results/section_4/TNBC_BRCA1_0554_tumor_scRNA(tsne_by_cor).RData")
df1             <- as.data.frame(TNBC_BRCA1_0554_tumor_scRNA@reductions$tsne@cell.embeddings)
df1$cell.id     <- rownames(df1)
subtype         <- as.data.frame(TNBC_BRCA1_0554_tumor_scRNA@meta.data)
subtype$cell.id <- rownames(subtype)
df1             <- merge(df1,subtype,by="cell.id")
df1             <- subset(df1, subtype=="Basal"|subtype=="Claudin_low")
table(df1$subtype)
fig_4c              <- ggplot(data = df1,aes(x=tSNE_1,y=tSNE_2,col=subtype))+
  geom_point(size=5)+xlab("")+ylab("")+
  ggplot.style+theme(plot.title = element_text(hjust=0.5,size = 44))+scale_color_manual(values = c("Red","Blue"))
ggsave(fig_4c,filename = "Pro_TNBC/paper/plot/section_4/TNBC_BRCA1_0554/tsne.plot.of.TNBC_BRCA1_0554_tumor.by.using.all.genes.pdf",width=23,height=15)





####*Figure 5b: the valcano plot of different expression genes in TNBC 0554####
load("Pro_TNBC/paper/data/results/section_5/de_results_TNBC0554_deseq2.RData")
library(ggrepel)
fig5a                                     <- ggplot(de_results_TNBC0554,aes(avg_log2FC, -log10(p_val_adj),color=sig))+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black")+
  geom_vline(xintercept = c(-0.5,0.5), linetype = "dashed", color = "black")+
  geom_point(size=5)+
  scale_color_manual(values=c("#2f5688","#BBBBBB","#CC0000"))+
  scale_size_continuous(range = c(0,1))+
  theme(legend.title = element_blank() )+
  geom_label_repel(data = subset(de_results_TNBC0554, SYMBOL %in% marker_genes),
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

ggsave(fig5a,filename = "Pro_TNBC/paper/plot/section_4/TNBC_BRCA1_0554/Volcano.map.of.0554.pdf",width = 20,height = 15)

####*Figure 5c: GSEA plot for EMT signature in DE genes in TNBC 0554####
library(clusterProfiler)
library(enrichplot)
library(msigdbr)
alldiff               <- subset(de_results_TNBC0554,sig=="up"|sig=="down")
alldiff               <- alldiff[order(alldiff$avg_log2FC,decreasing = T),]

genelist              <- alldiff$avg_log2FC
names(genelist)       <- alldiff$SYMBOL
hallmark              <- msigdbr(species = "Homo sapiens", category = "H")
hallmark              <- hallmark[,c(3,4)]
colnames(hallmark)    <- c('term','gene')
hallmark              <- subset(hallmark,term=="HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")
hallmark$term         <- "EMT"
TNBC0554.gsea.re              <- GSEA(genelist,  
                                             TERM2GENE = hallmark,  
                                             pvalueCutoff = 1 , 
                                             pAdjustMethod='BH')  
p1 <- gseaplot2(TNBC0554.gsea.re,
                TNBC0554.gsea.re$ID,
                color = "red",
                base_size = 20,
                rel_heights = c(1.5, 0.5, 1),
                subplots =1,
                ES_geom = "line",
                pvalue_table = F)
p2 <- gseaplot2(TNBC0554.gsea.re,
                TNBC0554.gsea.re$ID,
                color = "red",
                base_size = 20,
                rel_heights = c(1.5, 0.5, 1),
                subplots =2,
                ES_geom = "line",
                pvalue_table = F)
p3 <- plot_grid(p1,p2,ncol = 1)
ggsave(p3,filename ="Pro_TNBC/paper/plot/section_4/TNBC_BRCA1_0554/gsea.EMT.pdf",width =12,height =12)
save(TNBC0554.gsea.re,file = "Pro_TNBC/paper/data/results/section_4/TNBC0554.gsea.re.RData")


####*Figure 5d: CD44/CD24 ratio BL and CL cells in GSE173634 HDQP1####
library(ggpubr)
load("Pro_TNBC/paper/data/results/section_4/TNBC_BRCA1_0554_tumor_scRNA.RData")
TNBC0554.cpm                     <- as.matrix(TNBC_BRCA1_0554_tumor_scRNA@assays$RNA@data) %>% t() %>% as.data.frame()
TNBC0554.CD44.CD24.ratio         <- TNBC0554.cpm[,c("ENSG00000026508","ENSG00000272398")] 
colnames(TNBC0554.CD44.CD24.ratio)[1:2]   <- c("CD44","CD24")
TNBC0554.CD44.CD24.ratio$subtype <- TNBC_BRCA1_0554_tumor_scRNA$subtype
TNBC0554.CD44.CD24.ratio$ratio   <- (TNBC0554.CD44.CD24.ratio$CD44+0.1)/(TNBC0554.CD44.CD24.ratio$CD24+0.1)
table(TNBC0554.CD44.CD24.ratio$subtype)
TNBC0554.CD44.CD24.ratio  <- subset(TNBC0554.CD44.CD24.ratio,subtype=="Basal"|subtype=="Claudin_low")
save(TNBC0554.CD44.CD24.ratio,file = "Pro_TNBC/paper/data/results/section_4/TNBC0554.CD44.CD24.ratio.RData")

fig_4d     <- ggboxplot(TNBC0554.CD44.CD24.ratio,x = "subtype",y = "ratio",add = "point",size= 1,add.params=list(size=3))+ 
  ggplot.style+
  theme(panel.grid.major=element_line(colour=NA),
        panel.grid.minor = element_blank())+
  xlab("")+ylab("CD44/CD24 ratio")+stat_compare_means(method = "wilcox.test")+
  scale_x_discrete("", labels = c("Basal" = paste("Basal-like",sep = "\n","(n=1961)"),"Claudin_low" = paste("Claudin-low",sep = "\n","(n=170)")))
ggsave(fig_4d,filename = "Pro_TNBC/paper/plot/section_4/TNBC_BRCA1_0554/fig4d_boxplot.of.ratio.in.TNBC0554.pdf",width = 15,height = 15)

####*Figure 5e: correlation coefficents between claudin-low cells  and CCLE ####
load("Pro_TNBC/paper/data/results/section_4/cor.with.CCLE.RData")
CL.df                 <- CL_rs$correlation.matrix %>% t()%>% as.data.frame()
CL.df$CCLE_Name       <- rownames(CL.df)
brca_info             <- read_csv("Pro_TNBC/output/data/CCLE/brca_info.csv")
brca_info             <- brca_info[!is.na(brca_info$lineage_molecular_subtype),]
CL.df                 <- merge(CL.df,brca_info[,c(2,5,22)],by="CCLE_Name")
CL.df$subtype         <- ifelse(CL.df$lineage_molecular_subtype=="basal_B","Claudin-low",
                                ifelse(CL.df$lineage_molecular_subtype=="basal_A","Basal-like","Others"))
CL.df[is.na(CL.df$subtype),7]  <- "Others"
CL.df$TNBC0554_CL_rank              <- rank(CL.df$TNBC0554_CL) 
CL.df                 <- CL.df[order(CL.df$TNBC0554_CL_rank,decreasing = F),]
range(CL.df$TNBC0554_CL)
fig5e <- ggplot(CL.df, aes(x =TNBC0554_CL_rank , y = TNBC0554_CL,color=subtype)) +
  geom_point(size = 6) +
  ylim(0,0.55)+
  labs(x = "Rank", y = "Correlation coefficients") +
  scale_color_manual(values = c("Basal-like" = "red","Claudin-low"="blue" ,"others" = "black"))  +  
  theme_bw()+
  geom_point(data = subset(CL.df, lineage_molecular_subtype == "basal_B"), color = "blue", size = 6) +
  geom_point(data = subset(CL.df, lineage_molecular_subtype == "basal_A"), color = "red", size = 6)+ggplot.style+
  theme(panel.grid = element_blank())+
  geom_text_repel(data = CL.df[46:50,], 
                  aes(label = stripped_cell_line_name), 
                  nudge_y = 0.04,  # 标签在纵坐标上的偏移量
                  color = "black",
                  size = 6,
                  arrow = arrow(length = unit(0.02, "npc")))
ggsave(fig5e,filename = "Pro_TNBC/paper/plot/section_4/TNBC_BRCA1_0554/cor.between.CCLE.and.TNBC0554.CL.pdf",width = 20,height = 15)

####*Figure 5f: correlation coefficents between basal-like cells  and CCLE ####
load("Pro_TNBC/paper/data/results/section_4/cor.with.CCLE.RData")
BL.df                 <- BL_rs$correlation.matrix %>% t()%>% as.data.frame()
BL.df$CCLE_Name       <- rownames(BL.df)
brca_info             <- read_csv("Pro_TNBC/output/data/CCLE/brca_info.csv")
brca_info             <- brca_info[!is.na(brca_info$lineage_molecular_subtype),]
BL.df                 <- merge(BL.df,brca_info[,c(2,5,22)],by="CCLE_Name")
BL.df$subtype         <- ifelse(BL.df$lineage_molecular_subtype=="basal_B","Claudin-low",
                                ifelse(BL.df$lineage_molecular_subtype=="basal_A","Basal-like","Others"))
BL.df[is.na(BL.df$subtype),7]  <- "Others"
BL.df$TNBC0554_BL_rank              <- rank(BL.df$TNBC0554_BL) 
BL.df                 <- BL.df[order(BL.df$TNBC0554_BL_rank,decreasing = F),]
range(BL.df$TNBC0554_BL)
fig5f                 <- ggplot(BL.df, aes(x = TNBC0554_BL_rank, y = TNBC0554_BL,color=subtype)) +
  geom_point(size = 6) +  
  ylim(0,0.55)+
  labs(x = "Rank", y = "Correlation coefficients") +
  scale_color_manual(values = c("Basal-like" = "red","Claudin-low"="blue" ,"others" = "black"))  +  
  theme_bw()+
  geom_point(data = subset(BL.df, lineage_molecular_subtype == "basal_B"), color = "blue", size = 6) +
  geom_point(data = subset(BL.df, lineage_molecular_subtype == "basal_A"), color = "red", size = 6)+ggplot.style +
  theme(panel.grid = element_blank())+
  geom_text_repel(data = BL.df[46:50,], 
                  aes(label = stripped_cell_line_name), 
                  nudge_y = 0.04,  # The offset of the label on the vertical axis
                  color = "black",
                  size = 6,
                  arrow = arrow(length = unit(0.02, "npc")))
ggsave(fig5f,filename = "Pro_TNBC/paper/plot/section_4/TNBC_BRCA1_0554/cor.between.CCLE.and.TNBC0554.BL.pdf",width = 20,height = 15)

####supplementary Figure 3: feature plot of CLDN3,CLDN4,CLDN7, EPCAM, KRT19, VIM in GSE173634 HDQP1####
load("Pro_TNBC/paper/data/results/section_4/GSE173634_HDQP1_scRNA.RData")
library(ggplot2)
library(Seurat)
p1 <- FeaturePlot(GSE173634.HDQP1.scRNA,features = c("ENSG00000165215"),
                  reduction = "tsne",pt.size = 5)+
  theme( 
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
    axis.ticks = element_blank(),axis.text.x = element_blank(),axis.title.y = element_blank(),
    legend.position = "none", 
    axis.title = element_text( size=48, face="bold"),
    axis.text  = element_text( size=48, face="bold"),
    plot.title = element_text(size=48,hjust=0.5))+ggtitle("CLDN3")+theme_bw(base_size = 48) + theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))
p2 <- FeaturePlot(GSE173634.HDQP1.scRNA,features = c("ENSG00000189143"),
                  reduction = "tsne",pt.size = 5)+
  theme( 
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
    axis.ticks = element_blank(),axis.text.x = element_blank(),axis.title.y = element_blank(),
    legend.position = "none", 
    plot.title = element_text(size=48,hjust=0.5),
    axis.title = element_text( size=48, face="bold"),
    axis.text  = element_text( size=48, face="bold"),
    plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))+ggtitle("CLDN4")+theme_bw(base_size = 48) 
p3 <- FeaturePlot(GSE173634.HDQP1.scRNA,features = c("ENSG00000181885"),
                  reduction = "tsne",pt.size = 5)+
  theme( 
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
    axis.ticks = element_blank(),axis.text.x = element_blank(),axis.title.y = element_blank(),
    legend.position = "none", 
    plot.title = element_text(size=48,hjust=0.5),
    axis.title = element_text( size=48, face="bold"),
    axis.text  = element_text( size=48, face="bold"),
    plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))+ggtitle("CLDN7")+theme_bw(base_size = 48)
p4 <- FeaturePlot(GSE173634.HDQP1.scRNA,features = c("ENSG00000119888"),
                  reduction = "tsne",pt.size = 5)+
  theme( 
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
    axis.ticks = element_blank(),axis.text.x = element_blank(),axis.title.y = element_blank(),
    legend.position = "none", 
    axis.title = element_text( size=48, face="bold"),
    axis.text  = element_text( size=48, face="bold"),
    plot.title = element_text(size=48,hjust=0.5),
    plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))+ggtitle("EPCAM")+theme_bw(base_size = 48)
p5 <- FeaturePlot(GSE173634.HDQP1.scRNA,features = c("ENSG00000171345"),
                  reduction = "tsne",pt.size = 5)+
  theme( 
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
    axis.ticks = element_blank(),axis.text.x = element_blank(),axis.title.y = element_blank(),
    legend.position = "none", 
    axis.title = element_text( size=48, face="bold"),
    axis.text  = element_text( size=48, face="bold"),
    plot.title = element_text(size=48,hjust=0.5),
    plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))+ggtitle("KRT19")+theme_bw(base_size = 48)
p6 <- FeaturePlot(GSE173634.HDQP1.scRNA,features = c("ENSG00000026025"),
                  reduction = "tsne",pt.size = 5)+
  theme( 
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
    axis.ticks = element_blank(),axis.text.x = element_blank(),axis.title.y = element_blank(),
    legend.position = "none", 
    axis.title = element_text( size=48, face="bold"),
    axis.text  = element_text( size=48, face="bold"),
    plot.title = element_text(size=48,hjust=0.5),
    plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))+ggtitle("VIM")+theme_bw(base_size = 48)

figs3 <- plot_grid(p1,p2,p3,p4,p5,p6,ncol = 3)
ggsave(figs3,filename = "Pro_TNBC/paper/plot/section_4/HDQP1_GSE173634/feature.plot.in.GSE173634.HDQP1.pdf",width=40,height=25)




library(ggplot2)
library(ggrepel)
source("Pro_TNBC/paper/code/ggplot.style.R")
####Supplementary Figure 4: Discovery of “Basal-like/Claudin-low mixed” subtype in HDQP1 cell line (from GSE202771 dataset).#### 
####*Figure S4a: the t-SNE plot for HDQP1 data####
df1                              <- as.data.frame(GSE202771_HDQP1_scRNA@reductions$tsne@cell.embeddings)
df1$cell.id                      <- rownames(df1)
subtype                          <- as.data.frame(GSE202771_HDQP1_scRNA@meta.data)
subtype$cell.id                  <- rownames(subtype)
df1                              <- merge(df1,subtype,by="cell.id")
fig_4a                           <- ggplot(data = df1,aes(x=tSNE_1,y=tSNE_2,col=subtype))+
  geom_point(size=5)+
  ggplot.style+theme(plot.title = element_text(hjust=0.5,size = 44))+scale_color_manual(values = c(Basal="Red",Claudin_low="Blue"))
ggsave(fig_4a,filename = "Pro_TNBC/paper/plot/section_4/HDQP1_GSE202771/tsne.plot.of.HDQP1.pdf",width=23,height=15)


####*Figure S4b: the valcano plot of different expression genes in HDQP1 cell line (GSE173634)####
library(ggrepel)
load("Pro_TNBC/paper/data/results/section_5/de_results_HDQP1_GSE202771_deseq2.RData")
marker_genes                               <- c("CLDN3","CLDN4","CLDN7","EPCAM","VIM","CD24","CD44")

fig5b                                            <- ggplot(de_results_HDQP1_GSE202771,aes(avg_log2FC, -log10(p_val_adj),color=sig))+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black")+
  geom_vline(xintercept = c(-0.5,0.5), linetype = "dashed", color = "black")+
  geom_point(size=5)+
  scale_color_manual(values=c("#2f5688","#BBBBBB","#CC0000"))+
  scale_size_continuous(range = c(0,1))+
  theme(legend.title = element_blank() )+
  geom_label_repel(data = subset(de_results_HDQP1_GSE202771, SYMBOL %in% marker_genes),
                   max.overlaps = getOption("ggrepel.max.overlaps", default = 20),
                   aes(label = SYMBOL),
                   size = 8, 
                   color = 'black',
                   fontface="bold.italic",
                   segment.color = "black",   
                   segment.linetype = "solid", 
                   segment.size = 0.5 ,
                   box.padding = 3,        
                   point.padding = 0.2 ) +
  xlab("log2fc")+
  ylab("-log10p-adjust")+ggplot.style

ggsave(fig5b,filename = "Pro_TNBC/paper/plot/section_4/HDQP1_GSE202771/Volcano.map.of.GSE202771.HDQP1.pdf",width = 20,height = 15)

####*Figure S4c: GSEA plot for EMT signature in DE genes in GSE173634 HDQP1####
library(clusterProfiler)
library(enrichplot)
library(msigdbr)
alldiff               <- subset(de_results_HDQP1_GSE202771,sig=="up"|sig=="down")
alldiff               <- alldiff[!duplicated(alldiff$SYMBOL),]
alldiff               <- alldiff[order(alldiff$avg_log2FC,decreasing = T),]
genelist              <- alldiff$avg_log2FC
names(genelist)       <- alldiff$SYMBOL
hallmark              <- msigdbr(species = "Homo sapiens", category = "H")
hallmark              <- hallmark[,c(3,4)]
colnames(hallmark)    <- c('term','gene')
hallmark              <- subset(hallmark,term=="HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")
hallmark$term         <- "EMT"
GSE202771.HDQP1.gsea.re              <- GSEA(genelist,  
                                             TERM2GENE = hallmark,  
                                             pvalueCutoff = 1 , 
                                             pAdjustMethod='BH')  
pdf("Pro_TNBC/paper/plot/section_4/HDQP1_GSE202771/gsea.EMT.pdf",width =20,height =15)
gseaplot2(GSE202771.HDQP1.gsea.re,
          GSE202771.HDQP1.gsea.re$ID,
          color = "red",
          base_size = 20,
          rel_heights = c(1.5, 0.5, 1),
          subplots =1:2,
          ES_geom = "line",
          pvalue_table = F)
dev.off()
#p-value = 1.44e-07
save(GSE202771.HDQP1.gsea.re,file = "Pro_TNBC/paper/data/results/section_4/GSE202771.HDQP1.gsea.re.RDa")


####Figure S4d: CD44/CD24 ratio BL and CL cells in GSE173634 HDQP1####
library(ggpubr)
load("Pro_TNBC/paper/data/results/section_4/GSE202771_HDQP1_scRNA.RData")
GSE202771.HDQP1.cpm                     <- as.matrix(GSE202771_HDQP1_scRNA@assays$RNA@data) %>% t() %>% as.data.frame()
GSE202771.HDQP1.CD44.CD24.ratio         <- GSE202771.HDQP1.cpm[,c("CD44"),drop=F] #no CD24
colnames(GSE202771.HDQP1.CD44.CD24.ratio)[1:2]   <- c("CD44","CD24")
GSE202771.HDQP1.CD44.CD24.ratio$subtype <- GSE202771_HDQP1_scRNA$subtype
GSE202771.HDQP1.CD44.CD24.ratio$ratio   <- (GSE202771.HDQP1.CD44.CD24.ratio$CD44+0.1)/(GSE202771.HDQP1.CD44.CD24.ratio$CD24+0.1)
save(GSE202771.HDQP1.CD44.CD24.ratio,file = "Pro_TNBC/paper/data/results/section_4/GSE202771.HDQP1.CD44.CD24.ratio.RData")

fig_4d     <- ggboxplot(GSE202771.HDQP1.CD44.CD24.ratio,x = "subtype",y = "CD44",add = "point",size= 1,add.params=list(size=3))+ 
  ggplot.style+
  theme(panel.grid.major=element_line(colour=NA),
        panel.grid.minor = element_blank())+
  xlab("")+ylab("CD44")+stat_compare_means(method = "wilcox.test")+
  scale_x_discrete("", labels = c("Basal" = paste("Basal-like",sep = "\n","(n=6729)"),"Claudin_low" = paste("Claudin-low",sep = "\n","(n=239)")))
ggsave(fig_4d,filename = "Pro_TNBC/paper/plot/section_4/HDQP1_GSE202771/fig4d_boxplot.of.ratio.in.GSE202771.HDQP1.pdf",width = 15,height = 15)

####Figure S4e: correlation coefficents between claudin-low cells  and CCLE ####
load("Pro_TNBC/paper/data/results/section_4/cor.with.CCLE.RData")
CL.df                 <- CL_rs$correlation.matrix %>% t()%>% as.data.frame()
CL.df$CCLE_Name       <- rownames(CL.df)
brca_info             <- read_csv("Pro_TNBC/output/data/CCLE/brca_info.csv")
brca_info             <- brca_info[!is.na(brca_info$lineage_molecular_subtype),]
CL.df                 <- merge(CL.df,brca_info[,c(2,5,22)],by="CCLE_Name")
CL.df$subtype         <- ifelse(CL.df$lineage_molecular_subtype=="basal_B","Claudin-low",
                                ifelse(CL.df$lineage_molecular_subtype=="basal_A","Basal-like","Others"))
CL.df[is.na(CL.df$subtype),7]  <- "Others"
CL.df$HDQP1_GSE202771_CL_rank               <- rank(CL.df$HDQP1_GSE202771_CL) 
CL.df                 <- CL.df[order(CL.df$HDQP1_GSE202771_CL_rank,decreasing = F),]
range(CL.df$HDQP1_GSE202771_CL)
fig1e <- ggplot(CL.df, aes(x = HDQP1_GSE202771_CL_rank, y = HDQP1_GSE202771_CL,color=subtype)) +
  geom_point(size = 6) + ylim(0,0.9) +
  labs(x = "Rank", y = "Correlation coefficient") +
  scale_color_manual(values = c("Basal-like" = "red","Claudin-low"="blue" ,"others" = "black"))  +  
  theme_bw()+
  geom_point(data = subset(CL.df, lineage_molecular_subtype == "basal_B"), color = "blue", size = 6) +
  geom_point(data = subset(CL.df, lineage_molecular_subtype == "basal_A"), color = "red", size = 6)+ggplot.style +
  theme(panel.grid = element_blank())+
  geom_text_repel(data = CL.df[46:50,], 
                  aes(label = stripped_cell_line_name), 
                  nudge_y = 0.06,  
                  color = "black",
                  fontface="bold",
                  size = 6,
                  arrow = arrow(length = unit(0.02, "npc")))
ggsave(fig1e,filename = "Pro_TNBC/paper/plot/section_4/HDQP1_GSE202771/cor.between.CCLE.and.GSE202771.HDQP1.CL.pdf",width = 20,height = 15)

####Figure S4f: correlation coefficents between basal-like cells  and CCLE ####
load("Pro_TNBC/paper/data/results/section_4/cor.with.CCLE.RData")
BL.df                 <- BL_rs$correlation.matrix %>% t()%>% as.data.frame()
BL.df$CCLE_Name       <- rownames(BL.df)
brca_info             <- read_csv("Pro_TNBC/output/data/CCLE/brca_info.csv")
brca_info             <- brca_info[!is.na(brca_info$lineage_molecular_subtype),]
BL.df                 <- merge(BL.df,brca_info[,c(2,5,22)],by="CCLE_Name")
BL.df$subtype         <- ifelse(BL.df$lineage_molecular_subtype=="basal_B","Claudin-low",
                                ifelse(BL.df$lineage_molecular_subtype=="basal_A","Basal-like","Others"))
BL.df[is.na(BL.df$subtype),7]  <- "Others"
BL.df$HDQP1_GSE202771_BL_rank               <- rank(BL.df$HDQP1_GSE202771_BL) 
BL.df                 <- BL.df[order(BL.df$HDQP1_GSE202771_BL_rank,decreasing = F),]
range(BL.df$HDQP1_GSE202771_BL)
fig4f                 <- ggplot(BL.df, aes(x = HDQP1_GSE202771_BL_rank, y = HDQP1_GSE202771_BL,color=subtype)) +
  geom_point(size = 6) + ylim(0,0.9)+ 
  labs(x = "Rank", y = "Correlation coefficient") +
  scale_color_manual(values = c("Basal-like" = "red","Claudin-low"="blue" ,"others" = "black"))  +  
  theme_bw()+
  geom_point(data = subset(BL.df, lineage_molecular_subtype == "basal_B"), color = "blue", size = 6) +
  geom_point(data = subset(BL.df, lineage_molecular_subtype == "basal_A"), color = "red", size = 6)+ggplot.style +
  theme(panel.grid = element_blank())+
  geom_text_repel(data = BL.df[46:50,], 
                  aes(label = stripped_cell_line_name), 
                  nudge_y = 0.06,  
                  color = "black",
                  fontface="bold",
                  size = 6,
                  arrow = arrow(length = unit(0.02, "npc")))
ggsave(fig4f,filename = "Pro_TNBC/paper/plot/section_4/HDQP1_GSE202771/cor.between.CCLE.and.GSE202771.HDQP1.BL.pdf",width = 20,height = 15)


####supplementary Figure 5: feature plot of CLDN3,CLDN4,CLDN7, EPCAM, KRT19, VIM in GSE202771 HDQP1####
load("Pro_TNBC/paper/data/results/section_4/GSE202771_HDQP1_scRNA.RData")
table(GSE202771_HDQP1_scRNA$subtype)
library(ggplot2)
library(Seurat)
library(cowplot)
p1 <- FeaturePlot(GSE202771_HDQP1_scRNA,features = c("CLDN3"),
                  reduction = "tsne",pt.size = 5)+
  theme( 
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
    axis.ticks = element_blank(),axis.text.x = element_blank(),axis.title.y = element_blank(),
    legend.position = "none", 
    axis.title = element_text( size=48, face="bold"),
    axis.text  = element_text( size=48, face="bold"),
    plot.title = element_text(size=48,hjust=0.5))+ggtitle("CLDN3")+theme_bw(base_size = 48) + theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))
p2 <- FeaturePlot(GSE202771_HDQP1_scRNA,features = c("CLDN4"),
                  reduction = "tsne",pt.size = 5)+
  theme( 
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
    axis.ticks = element_blank(),axis.text.x = element_blank(),axis.title.y = element_blank(),
    legend.position = "none", 
    plot.title = element_text(size=48,hjust=0.5),
    axis.title = element_text( size=48, face="bold"),
    axis.text  = element_text( size=48, face="bold"),
    plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))+ggtitle("CLDN4")+theme_bw(base_size = 48) 
p3 <- FeaturePlot(GSE202771_HDQP1_scRNA,features = c("CLDN7"),
                  reduction = "tsne",pt.size = 5)+
  theme( 
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
    axis.ticks = element_blank(),axis.text.x = element_blank(),axis.title.y = element_blank(),
    legend.position = "none", 
    plot.title = element_text(size=48,hjust=0.5),
    axis.title = element_text( size=48, face="bold"),
    axis.text  = element_text( size=48, face="bold"),
    plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))+ggtitle("CLDN7")+theme_bw(base_size = 48)
p4 <- FeaturePlot(GSE202771_HDQP1_scRNA,features = c("EPCAM"),
                  reduction = "tsne",pt.size = 5)+
  theme( 
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
    axis.ticks = element_blank(),axis.text.x = element_blank(),axis.title.y = element_blank(),
    legend.position = "none", 
    axis.title = element_text( size=48, face="bold"),
    axis.text  = element_text( size=48, face="bold"),
    plot.title = element_text(size=48,hjust=0.5),
    plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))+ggtitle("EPCAM")+theme_bw(base_size = 48)
p5 <- FeaturePlot(GSE202771_HDQP1_scRNA,features = c("KRT19"),
                  reduction = "tsne",pt.size = 5)+
  theme( 
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
    axis.ticks = element_blank(),axis.text.x = element_blank(),axis.title.y = element_blank(),
    legend.position = "none", 
    axis.title = element_text( size=48, face="bold"),
    axis.text  = element_text( size=48, face="bold"),
    plot.title = element_text(size=48,hjust=0.5),
    plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))+ggtitle("KRT19")+theme_bw(base_size = 48)
p6 <- FeaturePlot(GSE202771_HDQP1_scRNA,features = c("VIM"),
                  reduction = "tsne",pt.size = 5)+
  theme( 
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
    axis.ticks = element_blank(),axis.text.x = element_blank(),axis.title.y = element_blank(),
    legend.position = "none", 
    axis.title = element_text( size=48, face="bold"),
    axis.text  = element_text( size=48, face="bold"),
    plot.title = element_text(size=48,hjust=0.5),
    plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))+ggtitle("VIM")+theme_bw(base_size = 48)

figs5 <- plot_grid(p1,p2,p3,p4,p5,p6,ncol = 3)
ggsave(figs5,filename = "Pro_TNBC/paper/plot/section_4/HDQP1_GSE202771/feature.plot.in.GSE202771.HDQP1.pdf",width=40,height=25)


####supplementary Figure 6: feature plot of CLDN3,CLDN4,CLDN7, EPCAM, KRT19, VIM in TNBC 0554####
load("Pro_TNBC/paper/data/results/section_4/TNBC_BRCA1_0554_tumor_scRNA(tsne_by_cor).RData")
table(TNBC_BRCA1_0554_tumor_scRNA$subtype)
TNBC_BRCA1_0554_tumor_scRNA  <- subset(TNBC_BRCA1_0554_tumor_scRNA,subtype=="Basal"|subtype=="Claudin_low")
library(ggplot2)
library(Seurat)
library(cowplot)
p1 <- FeaturePlot(TNBC_BRCA1_0554_tumor_scRNA,features = c("ENSG00000165215"),
                  reduction = "tsne",pt.size = 5)+
  theme( 
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
    axis.ticks = element_blank(),axis.text.x = element_blank(),axis.title.y = element_blank(),
    legend.position = "none", 
    plot.title = element_text(size=55,hjust=0.5))+ggtitle("CLDN3")+theme_bw(base_size = 55) + theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))
p2 <- FeaturePlot(TNBC_BRCA1_0554_tumor_scRNA,features = c("ENSG00000189143"),
                  reduction = "tsne",pt.size = 5)+
  theme( 
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
    axis.ticks = element_blank(),axis.text.x = element_blank(),axis.title.y = element_blank(),
    legend.position = "none", 
    plot.title = element_text(size=55,hjust=0.5),
    plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))+ggtitle("CLDN4")+theme_bw(base_size = 55) 
p3 <- FeaturePlot(TNBC_BRCA1_0554_tumor_scRNA,features = c("ENSG00000181885"),
                  reduction = "tsne",pt.size = 5)+
  theme( 
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
    axis.ticks = element_blank(),axis.text.x = element_blank(),axis.title.y = element_blank(),
    legend.position = "none", 
    plot.title = element_text(size=55,hjust=0.5),
    plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))+ggtitle("CLDN7")+theme_bw(base_size = 55)
p4 <- FeaturePlot(TNBC_BRCA1_0554_tumor_scRNA,features = c("ENSG00000119888"),
                  reduction = "tsne",pt.size = 5)+
  theme( 
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
    axis.ticks = element_blank(),axis.text.x = element_blank(),axis.title.y = element_blank(),
    legend.position = "none", 
    plot.title = element_text(size=55,hjust=0.5),
    plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))+ggtitle("EPCAM")+theme_bw(base_size = 55)
p5 <- FeaturePlot(TNBC_BRCA1_0554_tumor_scRNA,features = c("ENSG00000171345"),
                  reduction = "tsne",pt.size = 5)+
  theme( 
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
    axis.ticks = element_blank(),axis.text.x = element_blank(),axis.title.y = element_blank(),
    legend.position = "none", 
    plot.title = element_text(size=55,hjust=0.5),
    plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))+ggtitle("KRT19")+theme_bw(base_size = 55)
p6 <- FeaturePlot(TNBC_BRCA1_0554_tumor_scRNA,features = c("ENSG00000026025"),
                  reduction = "tsne",pt.size = 5)+
  theme( 
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
    axis.ticks = element_blank(),axis.text.x = element_blank(),axis.title.y = element_blank(),
    legend.position = "none", 
    plot.title = element_text(size=55,hjust=0.5),
    plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))+ggtitle("VIM")+theme_bw(base_size = 55)

figs3 <- plot_grid(p1,p2,p3,p4,p5,p6,ncol = 2)
ggsave(figs3,filename = "Pro_TNBC/paper/plot/section_4/TNBC_BRCA1_0554/feature.plot.in.TNBC0554.pdf",width=25,height=30)


                                  
