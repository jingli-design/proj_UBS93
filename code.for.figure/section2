#######################################################
####2.1: Validating UBS93 with bulk RNAseq data
#######################################################

#### 15 cell lines samples encoded as GSE212143####
library(readr)
library(GSVA)
library(dplyr)
load("Pro_TNBC/paper/data/method/UBS93.data.RData")
load("Pro_TNBC/data/CCLE/bulk_RNAseq/cell_lines/processed_by_myself/GSE212143.gene.expression.RData")
GSE212143.log2tpm                 <- log2.tpm.matrix
GSE212143_info                    <- read_csv("Pro_TNBC/data/CCLE/bulk_RNAseq/cell_lines/processed_by_myself/GSE212143_info.csv")
GSE212143_subtype                 <- breast.cancer.predictor(expr.of.sample = GSE212143.log2tpm,
                                                             expr.of.centroid = UBS93.data$UBS93.centroid,
                                                             marker.gene = UBS93.data$UBS93.gene.df$ENSEMBL,
                                                             HER2.amp.signature.genes = UBS93.data$HER2.amp.signature.genes,
                                                             ssgsea.cutoff = UBS93.data$ssgsea.cutoff)
GSE212143.subtype                 <- GSE212143_subtype$subtype
GSE212143.subtype$run_accession   <- rownames(GSE212143.subtype)
GSE212143.subtype                 <- merge(GSE212143.subtype,GSE212143_info,by="run_accession")
GSE212143.cor                     <- GSE212143_subtype$cor.matrix %>% as.data.frame()
GSE212143.cor$run_accession       <- rownames(GSE212143.cor)
GSE212143.subtype                 <- merge(GSE212143.subtype,GSE212143.cor,by="run_accession")
GSE212143.subtype                 <- GSE212143.subtype[,-3]

#### 62 cell lines samples encoded as GSE48213#### 
library(readr)
load("Pro_TNBC/data/CCLE/bulk_RNAseq/cell_lines/processed_by_myself/GSE48213.gene.expression.RData")
GSE48213_info                       <- read_csv("Pro_TNBC/data/CCLE/bulk_RNAseq/cell_lines/processed_by_myself/GSE48213_info.csv")
GSE48213_log2tpm                    <- log2.tpm.matrix
GSE48213_subtype                    <- breast.cancer.predictor(expr.of.sample = GSE48213_log2tpm,
                                                               expr.of.centroid = UBS93.data$UBS93.centroid,
                                                               marker.gene = UBS93.data$UBS93.gene.df$ENSEMBL,
                                                               HER2.amp.signature.genes = UBS93.data$HER2.amp.signature.genes,
                                                               ssgsea.cutoff = UBS93.data$ssgsea.cutoff)
GSE48213.subtype                 <- GSE48213_subtype$subtype
GSE48213.subtype$run_accession   <- rownames(GSE48213.subtype)
GSE48213.subtype                 <- merge(GSE48213.subtype,GSE48213_info,by="run_accession")
GSE48213.cor                     <- GSE48213_subtype$cor.matrix %>% as.data.frame()
GSE48213.cor$run_accession       <- rownames(GSE48213.cor)
GSE48213.subtype                 <- merge(GSE48213.subtype,GSE48213.cor,by="run_accession")
GSE48213.subtype                 <- GSE48213.subtype[,-3]
GSE48213.subtype                 <- subset(GSE48213.subtype,lineage_molecular_subtype=="luminal"|
                                             lineage_molecular_subtype=="HER2_amp"|
                                             lineage_molecular_subtype=="basal_A"|
                                             lineage_molecular_subtype=="basal_B")

#### 82 cell lines samples encoded as GSE73526#####
GSE73526_all_fpkm        <- read.delim("Pro_TNBC/data/CCLE/bulk_RNAseq/cell_lines/GSE73526_all_FPKM.txt",sep = "\t")
rown                     <- GSE73526_all_fpkm$ensembl_id
GSE73526_fpkm            <- GSE73526_all_fpkm[,-c(1:3)]  %>% as.matrix()
rownames(GSE73526_fpkm)  <- rown
GSE73526_log2fpkm        <- log2(GSE73526_fpkm +1 )
GSE73526_info            <- data.frame(sample.id=colnames(GSE73526_fpkm))
GSE73526_info$CCLE_Name  <- paste(GSE73526_info$sample.id,"_BREAST",sep = "")
GSE73526_info            <- left_join(GSE73526_info,ccle_info[,c(4,21)],by="CCLE_Name")
GSE73526_info[65:74,3]   <- c("basal","basal","basal_B","basal_B","luminal","basal_A","basal","basal","luminal","luminal")
save(GSE73526_info,file = "Pro_TNBC/data/CCLE/bulk_RNAseq/cell_lines/GSE73526_info.RData")
save(GSE73526_log2fpkm,file = "Pro_TNBC/data/CCLE/bulk_RNAseq/cell_lines/processed_by_myself/GSE73526_log2fpkm.RData")
load("Pro_TNBC/data/CCLE/bulk_RNAseq/cell_lines/GSE73526_info.RData")
load("Pro_TNBC/data/CCLE/bulk_RNAseq/cell_lines/processed_by_myself/GSE73526_log2fpkm.RData")
GSE73526_subtype                 <- breast.cancer.predictor(expr.of.sample = GSE73526_log2fpkm,
                                                               expr.of.centroid = UBS93.data$UBS93.centroid,
                                                               marker.gene = UBS93.data$UBS93.gene.df$ENSEMBL,
                                                               HER2.amp.signature.genes = UBS93.data$HER2.amp.signature.genes,
                                                               ssgsea.cutoff = UBS93.data$ssgsea.cutoff)
GSE73526.subtype                 <- GSE73526_subtype$subtype
GSE73526.subtype$sample.id       <- rownames(GSE73526.subtype)
GSE73526.subtype                 <- merge(GSE73526.subtype,GSE73526_info,by="sample.id")
GSE73526.cor                     <- GSE73526_subtype$cor.matrix %>% as.data.frame()
GSE73526.cor$sample.id           <- rownames(GSE73526.cor)
GSE73526.subtype                 <- merge(GSE73526.subtype,GSE73526.cor,by="sample.id")
GSE73526.subtype                 <- subset(GSE73526.subtype,lineage_molecular_subtype=="luminal"|
                                             lineage_molecular_subtype=="HER2_amp"|
                                             lineage_molecular_subtype=="basal_A"|
                                             lineage_molecular_subtype=="basal_B")

#### TCGA TNBC####
library(readxl)
org_clinical_patient_brca           <- read_excel("Pro_TNBC/data/TCGA/TCGA/org_clinical_patient_brca.xlsx")
load("~/Pro_TNBC/data/TCGA/TCGA/Breast Invasive Carcinoma.RData")
org_clinical_patient_brca$sample.id <- paste(org_clinical_patient_brca$bcr_patient_barcode,"-01",sep = "")
tcga_brca_meta                      <- org_clinical_patient_brca[,c(2,44,50,56,113)]          
tcga_brca_meta                      <- tcga_brca_meta[-c(1,2),]
tcga_tnbc_meta                      <- subset(tcga_brca_meta,er_status_by_ihc=="Negative"&pr_status_by_ihc=="Negative"&her2_status_by_ihc=="Negative")
tcga_tnbc_log2tpm                   <- log2.tpm.matrix[,colnames(log2.tpm.matrix)%in% tcga_tnbc_meta$sample.id]
tcga_tnbc_subtype                   <- breast.cancer.predictor(expr.of.sample = tcga_tnbc_log2tpm,
                                                               expr.of.centroid = UBS93.data$UBS93.centroid,
                                                               marker.gene = UBS93.data$UBS93.gene.df$ENSEMBL,
                                                               HER2.amp.signature.genes = UBS93.data$HER2.amp.signature.genes,
                                                               ssgsea.cutoff = UBS93.data$ssgsea.cutoff)
tcga_tnbc.subtype                   <- tcga_tnbc_subtype$subtype
tcga_tnbc.subtype$sample.id         <- rownames(tcga_tnbc.subtype)
tcga_tnbc.cor                       <- tcga_tnbc_subtype$cor.matrix %>% as.data.frame()
tcga_tnbc.cor$sample.id             <- rownames(tcga_tnbc.cor)
tcga_tnbc.subtype                   <- merge(tcga_tnbc.subtype,tcga_tnbc.cor,by="sample.id")


#################################################
####2.2 Validating UBS93 with scRNAseq data
#################################################


#### single cell RNAseq of cell lines(GSE173634)####
library(Seurat)
load("Pro_TNBC/output/data/CCLE/GSE173634_scRNA.RData")
CCLE.Name                       <- names(table(GSE173634_scRNA$orig.ident))
load("~/Pro_TNBC/output/data/CCLE/GSE173634_info.RData")
GSE173634_subtype_result        <- NULL
for (i in 1:length(CCLE.Name)) {
  A                  <- GSE173634_info[i,3] %>% as.character()
  molecular_subtype  <- GSE173634_info[i,2] %>% as.character()
  #normalized
  scRNA              <- subset(GSE173634_scRNA, orig.ident == A)
  scRNA              <- NormalizeData(scRNA)
  exprmat            <- as.matrix(scRNA@assays$RNA@data) # data after normalizing
  predictor.res      <- breast.cancer.predictor(expr.of.sample = exprmat,
                                                expr.of.centroid = UBS93.data$UBS93.centroid,
                                                marker.gene = UBS93.data$UBS93.gene.df$ENSEMBL,
                                                HER2.amp.signature.genes = UBS93.data$HER2.amp.signature.genes,
                                                ssgsea.cutoff = UBS93.data$ssgsea.cutoff)
  
  Subtype                    <- table(predictor.res$subtype) %>% as.data.frame() 
  colnames(Subtype)[2]       <- GSE173634_info[i,3] 
  if (is.null(GSE173634_subtype_result)) {
    GSE173634_subtype_result <- Subtype
  } else {
    GSE173634_subtype_result <- merge(GSE173634_subtype_result, Subtype,by="subtype",all=T)
  }
}

rown                                    <- as.character(GSE173634_subtype_result$subtype)
GSE173634_subtype_result                <- GSE173634_subtype_result[,-1]    
rownames(GSE173634_subtype_result)      <- rown
GSE173634_subtype_result                <- as.matrix(GSE173634_subtype_result)  %>% t()  %>% as.data.frame()
GSE173634_subtype_result                <- GSE173634_subtype_result[rownames(GSE173634_subtype_result) %in% GSE173634_info$ccle.name,]
GSE173634_info                          <- GSE173634_info[match(GSE173634_info$ccle.name,rownames(GSE173634_subtype_result)),]
GSE173634_subtype_result                <- cbind(GSE173634_info[,2],GSE173634_subtype_result)
GSE173634_subtype_result$total_numbers  <- apply(GSE173634_subtype_result[,2:5],1,function(x){sum(x,na.rm =T)})
GSE173634_subtype_basal                 <- subset(GSE173634_subtype_result,lineage_molecular_subtype=="basal_A")
GSE173634_subtype_basal$ratio           <- GSE173634_subtype_basal$Basal/GSE173634_subtype_basal$total_numbers
GSE173634_subtype_cl                    <- subset(GSE173634_subtype_result,lineage_molecular_subtype=="basal_B")
GSE173634_subtype_cl$ratio              <- GSE173634_subtype_cl$Claudin_low/GSE173634_subtype_cl$total_numbers
GSE173634_subtype_her2                  <- subset(GSE173634_subtype_result,lineage_molecular_subtype=="HER2_amp")
GSE173634_subtype_her2$ratio            <- GSE173634_subtype_her2$HER2_amp/GSE173634_subtype_her2$total_numbers
GSE173634_subtype_luminal               <- subset(GSE173634_subtype_result,lineage_molecular_subtype=="luminal")
GSE173634_subtype_luminal$ratio         <- GSE173634_subtype_luminal$Luminal/GSE173634_subtype_luminal$total_numbers
GSE173634_subtype                       <- rbind(GSE173634_subtype_basal,GSE173634_subtype_cl)
GSE173634_subtype                       <- rbind(GSE173634_subtype,GSE173634_subtype_her2)
GSE173634_subtype                       <- rbind(GSE173634_subtype,GSE173634_subtype_luminal)
GSE173634_subtype$celllines_name        <- rownames(GSE173634_subtype)
GSE173634_subtype                       <- arrange(GSE173634_subtype,lineage_molecular_subtype)
GSE173634_subtype[GSE173634_subtype$lineage_molecular_subtype=="basal_A",1]  <- c(rep("Basal",11))
GSE173634_subtype[GSE173634_subtype$lineage_molecular_subtype=="basal_B",1]  <- c(rep("Claudin_low",4))
GSE173634_subtype[GSE173634_subtype$lineage_molecular_subtype=="luminal",1]  <- c(rep("Luminal",8))
save(GSE173634_subtype,file = "Pro_TNBC/paper/data/results/section_2/fig2b.GSE17634.SUBTYPE.RData")
ggplot(GSE173634_subtype,aes(x=lineage_molecular_subtype,y=ratio))+geom_boxplot()
fig_2a             <- ggplot(GSE173634_subtype,aes(x=lineage_molecular_subtype,y=ratio))+
  geom_boxplot(outlier.shape=NA,lwd=1.2)+
  geom_point(size=4,position="jitter") + ggplot.style + 
  xlab('Subtype') + ylab('Ratio')+
  scale_x_discrete("subtype", labels = c("Basal" = paste("Basal",sep = "\n","(n=13)"),
                                         "Claudin_low" = paste("Claudin_low",sep = "\n","(n=4)"),
                                         "HER2_amp"=paste("HER2_amp",sep = "\n","(n=5)"),
                                         "Luminal"=paste("Luminal",sep = "\n","(n=9)")))
ggsave(fig_2a,filename = "Pro_TNBC/paper/plot/section_2/fig_2b.GSE173634.subtype.pdf",width = 20,height = 15)

#### scRNAseq data of cell lines (GSE202771)####
library(GSVA)
table(GSE202771_scRNA$orig.ident)
CCLE.Name                       <- names(table(GSE202771_scRNA$orig.ident))
GSE202771_info                  <- data.frame(sample.id=CCLE.Name)
GSE202771_info$CCLE.Name        <- apply(GSE202771_info[,1,drop=F],1,function(x){strsplit(x, "_")[[1]][1]})
GSE202771_info$CCLE_Name        <- paste(GSE202771_info$CCLE.Name,"_BREAST",sep ="")
GSE202771_infor                 <- left_join(GSE202771_info,ccle_info[,c(4,21)],by="CCLE_Name")
GSE202771_infor[c(13,16,17,21:23),4]           <- c("basal_A","basal_B","basal_B","basal_A","luminal","basal_A")
GSE202771_infor[c(13,21:23),4]                 <- c(NA,"basal","luminal","basal")
save(GSE202771_infor,file = "Pro_TNBC/data/CCLE/scRNAseq/GSE202771_infor.RData")

load("Pro_TNBC/output/data/scRNASeq/GSE202771/GSE202771_scRNA.RData")
load("Pro_TNBC/data/CCLE/scRNAseq/GSE202771_infor.RData")
GSE202771_subtype_result        <- NULL
for (i in 1:length(CCLE.Name)) {
  A                  <- GSE202771_infor[i,1] %>% as.character()
  molecular_subtype  <- GSE202771_infor[i,4] %>% as.character()
  #normalized
  scRNA              <- subset(GSE202771_scRNA, orig.ident == A)
  scRNA              <- NormalizeData(scRNA)
  exprmat            <- as.matrix(scRNA@assays$RNA@data) # data after normalizing
  exprmat_UBS93                   <- as.data.frame(exprmat)
  exprmat_UBS93$SYMBOL            <- rownames(exprmat_UBS93)
  exprmat_UBS93                   <- merge(UBS93.data$UBS93.gene.df,exprmat_UBS93,by="SYMBOL")
  rown                            <- exprmat_UBS93$ENSEMBL
  exprmat_UBS93                   <- exprmat_UBS93[,-c(1:3)]
  rownames(exprmat_UBS93)         <- rown
  exprmat_UBS93                   <- as.matrix(exprmat_UBS93)
  predictor.res      <- breast.cancer.predictor(expr.of.sample = exprmat_UBS93,
                                                expr.of.centroid = UBS93.data$UBS93.centroid,
                                                marker.gene = UBS93.data$UBS93.gene.df$ENSEMBL,
                                                HER2.amp.signature.genes = UBS93.data$HER2.amp.signature.genes,
                                                ssgsea.cutoff = UBS93.data$ssgsea.cutoff)
  
  Subtype                    <- table(predictor.res$subtype) %>% as.data.frame() 
  colnames(Subtype)[2]       <- GSE202771_infor[i,1] 
  if (is.null(GSE202771_subtype_result)) {
    GSE202771_subtype_result <- Subtype
  } else {
    GSE202771_subtype_result <- merge(GSE202771_subtype_result, Subtype,by="subtype",all=T)
  }
}

rown                                    <- as.character(GSE202771_subtype_result$subtype)
GSE202771_subtype_result                <- GSE202771_subtype_result[,-1]    
rownames(GSE202771_subtype_result)      <- rown
GSE202771_subtype_result                <- as.matrix(GSE202771_subtype_result)  %>% t()  %>% as.data.frame()
GSE202771_infor                         <- GSE202771_infor[match(GSE202771_infor$sample.id,rownames(GSE202771_subtype_result)),]
GSE202771_subtype_result                <- cbind(GSE202771_infor[,4],GSE202771_subtype_result)
GSE202771_subtype_result$total_numbers  <- apply(GSE202771_subtype_result[,2:5],1,function(x){sum(x,na.rm =T)})
colnames(GSE202771_subtype_result)[1]   <- "lineage_molecular_subtype"
GSE202771_subtype_basal                 <- subset(GSE202771_subtype_result,lineage_molecular_subtype=="basal_A")
GSE202771_subtype_basal$ratio           <- GSE202771_subtype_basal$Basal/GSE202771_subtype_basal$total_numbers
GSE202771_subtype_cl                    <- subset(GSE202771_subtype_result,lineage_molecular_subtype=="basal_B")
GSE202771_subtype_cl$ratio              <- GSE202771_subtype_cl$Claudin_low/GSE202771_subtype_cl$total_numbers
GSE202771_subtype_luminal               <- subset(GSE202771_subtype_result,lineage_molecular_subtype=="luminal")
GSE202771_subtype_luminal$ratio         <- GSE202771_subtype_luminal$Luminal/GSE202771_subtype_luminal$total_numbers
GSE202771_subtype                       <- rbind(GSE202771_subtype_basal,GSE202771_subtype_cl)
GSE202771_subtype                       <- rbind(GSE202771_subtype,GSE202771_subtype_luminal)
GSE202771_subtype$celllines_name        <- rownames(GSE202771_subtype)
GSE202771_subtype[GSE202771_subtype$lineage_molecular_subtype=="basal_A",1]  <- c(rep("Basal",7))
GSE202771_subtype[GSE202771_subtype$lineage_molecular_subtype=="basal_B",1]  <- c(rep("Claudin_low",6))
GSE202771_subtype[GSE202771_subtype$lineage_molecular_subtype=="luminal",1]  <- c(rep("Luminal",2))
save(GSE202771_subtype,file = "Pro_TNBC/paper/data/results/section_2/supfig_4.GSE202771.SUBTYPE.RData")
ggplot(GSE202771_subtype,aes(x=lineage_molecular_subtype,y=ratio))+geom_boxplot()
supfig_4             <- ggplot(GSE202771_subtype,aes(x=lineage_molecular_subtype,y=ratio))+
  geom_boxplot(outlier.shape=NA,lwd=1.2)+
  geom_point(size=4,position="jitter") + ggplot.style + 
  xlab('Subtype') + ylab('Ratio')+
  scale_x_discrete("subtype", labels = c("Basal" = paste("Basal",sep = "\n","(n=10)"),
                                         "Claudin_low" = paste("Claudin_low",sep = "\n","(n=6)"),
                                         "Luminal"=paste("Luminal",sep = "\n","(n=2)")))
ggsave(supfig_4,filename = "Pro_TNBC/paper/plot/section_2/supfig_4.GSE202771.subtype.pdf",width = 20,height = 15)



#####*Figure S6: The expression values of ERBB2 in breast cancer cell lines of CCLE.####
library(readr)
ERBB2_Expression_Public_23Q4          <- read_csv("Pro_TNBC/data/CCLE/ERBB2 Expression Public 23Q4.csv")
ERBB2_Expression_BRCA                 <- subset(ERBB2_Expression_Public_23Q4,Lineage=="Breast")
ERBB2_Expression_BRCA                 <- subset(ERBB2_Expression_BRCA,`Lineage Subtype`=="Invasive Breast Carcinoma")
colnames(ERBB2_Expression_BRCA)[3]    <- "stripped_cell_line_name"
ERBB2_Expression_BRCA                 <- merge(ERBB2_Expression_BRCA,ccle_info[,c(3,4,21)],by="stripped_cell_line_name")
ERBB2_Expression_BRCA                 <- ERBB2_Expression_BRCA[!is.na(ERBB2_Expression_BRCA$lineage_molecular_subtype),]
colnames(ERBB2_Expression_BRCA)[8]    <- "subtype"
ERBB2_Expression_BRCA                 <- subset(ERBB2_Expression_BRCA,subtype!= "basal")
ERBB2_Expression_BRCA[ERBB2_Expression_BRCA$subtype=="basal_A",8]  <- c(rep("Basal",14))
ERBB2_Expression_BRCA[ERBB2_Expression_BRCA$subtype=="basal_B",8]  <- c(rep("Claudin_low",11))
ERBB2_Expression_BRCA[ERBB2_Expression_BRCA$subtype=="luminal",8]  <- c(rep("Luminal",14))
supfig_4              <- ggplot(ERBB2_Expression_BRCA,aes(x=Lineage,y=`Expression Public 23Q4`))+
  geom_boxplot(outlier.shape=NA,lwd=1.2)+
  geom_point(size=4,position="jitter",aes(colour = factor(subtype))) + ggplot.style + 
  xlab('') + ylab('ERBB2[log2(TPM+1)]')+scale_color_manual(values = c("Red","#FFCC00","Blue","Purple","#996633"))
ggsave(supfig_4,filename = "Pro_TNBC/paper/plot/section_2/supfig_5.ERBB2.in.CCLE.pdf",width = 20,height = 15)
save(ERBB2_Expression_BRCA,file = "Pro_TNBC/paper/data/results/section_2/ERBB2_Expression_BRCA.RData")


###########################################################
####2.3 Validating UBS93 with mouse bulk RNAseq data
###########################################################

load("Pro_TNBC/data/mouse/mouse_model/gene.expression.RData")
GSE157333_log2tpm             <- log2.tpm.matrix
load("~/Pro_TNBC/data/mouse/mouse_claudin_low_ref/mouse_claudin_low_ref.RData")
CL_ref_log2tpm                <- log2.tpm.matrix

#### Heat map used to distinguish between Basal and Cl in mice####
identical(rownames(GSE157333_log2tpm),rownames(CL_ref_log2tpm))
mouse_log2tpm                 <- cbind(GSE157333_log2tpm,CL_ref_log2tpm)
gene_mart                     <- read.delim("Pro_TNBC/data/mouse/mart_export.txt",sep = ",")
gene.id                       <- c("Cldn3","Cldn13","Cldn7","Vim","Epcam")
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

#### run our predictor for RNAseq data of mouse data ####
UBS93.mouse.gene                                      <- gene_mart[gene_mart$Gene.stable.ID %in% UBS93.data$UBS93.gene.df$ENSEMBL,]
UBS93.mouse.gene                                      <- UBS93.mouse.gene[!is.na(UBS93.mouse.gene$Mouse.gene.stable.ID),]
mouse_log2tpm_UBS93                                   <- mouse_log2tpm[rownames(mouse_log2tpm) %in% UBS93.mouse.gene$Mouse.gene.stable.ID,]
mouse_log2tpm_UBS93                                   <- as.data.frame(mouse_log2tpm_UBS93)
mouse_log2tpm_UBS93$Mouse.gene.stable.ID              <- rownames(mouse_log2tpm_UBS93)
mouse_log2tpm_UBS93                                   <- merge(UBS93.mouse.gene,mouse_log2tpm_UBS93,by="Mouse.gene.stable.ID")
rown                                                  <- mouse_log2tpm_UBS93$Gene.stable.ID
mouse_log2tpm_UBS93                                   <- mouse_log2tpm_UBS93[,-(1:5)] %>% as.matrix()
rownames(mouse_log2tpm_UBS93)                         <- rown
mouse_data_subtype                                    <-  breast.cancer.predictor(expr.of.sample = mouse_log2tpm_UBS93,
                                                                                  expr.of.centroid = UBS93.data$UBS93.centroid,
                                                                                  marker.gene = UBS93.data$UBS93.gene.df$ENSEMBL,
                                                                                  HER2.amp.signature.genes = UBS93.data$HER2.amp.signature.genes,
                                                                                  ssgsea.cutoff = UBS93.data$ssgsea.cutoff)
mouse.subtype                                         <- mouse_data_subtype$subtype
mouse.subtype$run_accession                           <- rownames(mouse.subtype)
colnames(mouse.subtype)[1]                            <- "UBS93_subtype"
mouse.subtype                                         <- merge(mouse.subtype,mouse_log2tpm_info,by="run_accession")
save(mouse.subtype,file="Pro_TNBC/paper/data/results/section_2/mouse.subtype.RData")

