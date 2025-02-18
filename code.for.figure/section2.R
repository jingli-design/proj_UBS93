####################################################################
####Figure 2. Evaluation of UBS93 classifier in human and mouse.####
####################################################################

####Figure 2a: Prediction outcomes of the UBS93 classifier on an assembled bulk RNA-seq dataset of CCLE breast cancer cell lines.####
subtype                                    <- c("Basal-like","Claudin-low","HER2-amp","Luminal")
CCLE.UBS93.confusion.matrix                <- table(CCLE_brca.subtype$lineage_molecular_subtype,CCLE_brca.subtype$UBS93.subtype) %>% as.matrix()
rownames(CCLE.UBS93.confusion.matrix)[1:4] <- subtype
colnames(CCLE.UBS93.confusion.matrix)[1:4] <- subtype
CCLE.UBS93.confusion.matrix.norm           <- apply(CCLE.UBS93.confusion.matrix,1,function(x){x/sum(x)})
CCLE.UBS93.confusion.matrix.norm           <- t(CCLE.UBS93.confusion.matrix.norm)
library(pheatmap)
p2          <-pheatmap(CCLE.UBS93.confusion.matrix.norm,cluster_rows = F,
                       cluster_cols = F,
                       color = colorRampPalette(c("navy","white","firebrick3"))(100),
                       show_colnames = T,border_color = NA,
                       scale = "none",
                       show_rownames = T,
                       fontsize=30,cellwidth = 55, cellheight = 55,display_numbers = T,
                       number_format = "%.2f",fontsize_number = 20, number_color = "white",
                       angle_col = 0,breaks = seq(0,1,length.out=101)
)
ggsave(p2,filename = "Pro_TNBC/paper/plot/section_2/CCLE.USB93.confusion.matrix.pdf",width = 20,height = 15)

####Figure 2b: Prediction outcomes of the UBS93 classifier on bulk RNA-seq data of TCGA breast cancer samples.####
subtype                                    <- c("Basal-like","Claudin-low","HER2-amp","Luminal")
TCGA.UBS93.confusion.matrix                <- table(tcga.brca.subtype.compare$Imputed.Subtype,tcga.brca.subtype.compare$UBS93.subtype) %>% as.matrix()
rownames(TCGA.UBS93.confusion.matrix)[1:4] <- subtype
colnames(TCGA.UBS93.confusion.matrix)[1:4] <- subtype
TCGA.UBS93.confusion.matrix.norm           <- apply(TCGA.UBS93.confusion.matrix,1,function(x){x/sum(x)})
TCGA.UBS93.confusion.matrix.norm           <- t(TCGA.UBS93.confusion.matrix.norm)
library(pheatmap)
p1          <-pheatmap(TCGA.UBS93.confusion.matrix.norm,cluster_rows = F,
                       cluster_cols = F,
                       color = colorRampPalette(c("navy","white","firebrick3"))(100),
                       show_colnames = T,border_color = NA,
                       scale = "none",
                       show_rownames = T,
                       fontsize=30,cellwidth = 55, cellheight = 55,display_numbers = T,
                       number_format = "%.2f",fontsize_number = 20, number_color = "white",
                       angle_col = 0,breaks = seq(0,1,length.out=101)
)
ggsave(p1,filename = "Pro_TNBC/paper/plot/section_2/TCGA.USB93.confusion.matrix.pdf",width = 20,height = 15)


####Figure 2c: Performance of the UBS93 classifier on an assembled scRNA-seq dataset.####
scRNAseq_UBS93_subtype$lineage_molecular_subtype  <- as.factor(scRNAseq_UBS93_subtype$lineage_molecular_subtype)
fig2a <- ggboxplot(scRNAseq_UBS93_subtype,x = "lineage_molecular_subtype",y = "ratio",add = "point",size= 1,add.params=list(size=5))+ 
  ggplot.style+
  theme(panel.grid.major=element_line(colour=NA),
        panel.grid.minor = element_blank())+ 
  xlab('') + ylab('Accuracy')+
  scale_x_discrete("subtype", labels = c("Basal" = paste("Basal-like",sep = "\n","(n=12)"),
                                         "Claudin_low" = paste("Claudin-low",sep = "\n","(n=5)"),
                                         "Luminal"=paste("Luminal",sep = "\n","(n=8)"),
                                         "HER2_amp" =paste("HER2-amp",sep = "\n","(n=5)")))
median_ratio <- scRNAseq_UBS93_subtype %>% group_by(lineage_molecular_subtype) %>% summarise(median_ratio=median(ratio))
ggsave(fig2c,filename = "Pro_TNBC/paper/plot/section_2/fig2c.scRNAseq.subtype.pdf",width = 20,height = 15)

####Figure 2d: Heat map used to distinguish between Basal and Cl in mice(Figure 2b)####
load("Pro_TNBC/data/mouse/mouse_model/gene.expression.RData")
GSE157333_log2tpm             <- log2.tpm.matrix
load("~/Pro_TNBC/data/mouse/mouse_claudin_low_ref/mouse_claudin_low_ref.RData")
CL_ref_log2tpm                <- log2.tpm.matrix
identical(rownames(GSE157333_log2tpm),rownames(CL_ref_log2tpm))
mouse_log2tpm                 <- cbind(GSE157333_log2tpm,CL_ref_log2tpm)
gene_mart                     <- read.delim("Pro_TNBC/data/mouse/mart_export.txt",sep = ",")
gene.id                       <- c("Cldn3","Cldn13","Cldn7","Epcam","Vim")
gene.5                        <- gene_mart[gene_mart$Mouse.gene.name %in% gene.id,4:5]
mouse_log2tpm_df              <- as.data.frame(mouse_log2tpm)
mouse_log2tpm_df$Mouse.gene.stable.ID               <- rownames(mouse_log2tpm_df)
mouse_log2tpm_df_5            <- merge(gene.5,mouse_log2tpm_df,by="Mouse.gene.stable.ID")
mouse_log2tpm_df_5            <- mouse_log2tpm_df_5[match(gene.id,mouse_log2tpm_df_5$Mouse.gene.name),]
identical(gene.id,mouse_log2tpm_df_5$Mouse.gene.name)
rown                          <- mouse_log2tpm_df_5$Mouse.gene.name
mouse_log2tpm_df_5            <- mouse_log2tpm_df_5[,-(1:2)]
rownames(mouse_log2tpm_df_5)  <- rown
mouse_log2tpm_df_5            <- as.matrix(mouse_log2tpm_df_5)
library(pheatmap)
fig_2c <- pheatmap(mouse_log2tpm_df_5,cluster_rows = F,
                   cluster_cols = T,
                   color = colorRampPalette(c("navy","white","firebrick3"))(100),
                   show_colnames = F,border_color = NA,
                   scale = "row",
                   show_rownames = T,
                   fontsize=30,cellwidth = 30, cellheight = 55
                   
)
col_order                   <- fig_2c$tree_col$order
mouse_log2tpm_info          <- data.frame(run_accession=colnames(mouse_log2tpm_df_5))
mouse_log2tpm_info          <- mouse_log2tpm_info[col_order,,drop=F]
mouse_log2tpm_info$subtype  <- c(rep("Basal",5),rep("Claudin_low",22))
#fontsize_row = 8,fontsize_col=12 ,
ggsave(fig_2c,filename = "Pro_TNBC/paper/plot/section_2/fig.2c.second.heat.map.used.to.distinguish.between.Basal.and.Cl.in.mice.pdf",width = 20,height = 15)
save(mouse_log2tpm_df_5,file = "Pro_TNBC/paper/data/results/section_2/fig.2c.mouse_log2tpm_df_5gene.RData")
save(mouse_log2tpm_info,file = "Pro_TNBC/paper/data/results/section_2/mouse_log2tpm_info.RData")








