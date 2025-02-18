####################################################
####Subtyping for celllines, TCGA BRCA and mouse####
####################################################

####the function used to calculate F1 score####
#pre：classified results predicted
#y：true classified results
f1_fun = function(pre,tru){
  class = sort(unique(tru))
  tp=NA
  fp=NA
  fn=NA
  for(i in 1:length(class)){
    tp[i] = sum(pre==class[i] & tru==class[i])
    fp[i] = sum(pre==class[i] & tru!=class[i])
    fn[i] = sum(pre!=class[i] & tru==class[i])
  }
  Precision <- tp / (tp + fp)  
  Recall    <- tp / (tp + fn)  
  F1        <- 2 * (Precision * Recall) / (Precision + Recall)  
  names(F1) = class
  return(F1)
}




####1.Cell lines from other papers(bulk)#####
library(readr)
GSE212143_info                    <- read_csv("Pro_TNBC/data/CCLE/bulk_RNAseq/cell_lines/processed_by_myself/GSE212143_info.csv")
GSE48213_info                     <- read_csv("Pro_TNBC/data/CCLE/bulk_RNAseq/cell_lines/processed_by_myself/GSE48213_info.csv")
load("Pro_TNBC/data/CCLE/bulk_RNAseq/cell_lines/GSE73526_info.RData")
load("Pro_TNBC/data/CCLE/bulk_RNAseq/cell_lines/processed_by_myself/GSE212143.gene.expression.RData")
GSE212143.log2tpm                 <- log2.tpm.matrix
load("Pro_TNBC/data/CCLE/bulk_RNAseq/cell_lines/processed_by_myself/GSE48213.gene.expression.RData")
GSE48213_log2tpm                  <- log2.tpm.matrix
load("Pro_TNBC/data/CCLE/bulk_RNAseq/cell_lines/processed_by_myself/GSE73526.gene.expression.RData")
GSE73526_log2tpm                  <- log2.tpm.matrix
load("Pro_TNBC/data/CCLE/CCLE.RData")
ccle_info                         <- read.csv("Pro_TNBC/data/CCLE/ccle_info.csv")
brca_meta                         <- CCLE.sample.meta[CCLE.sample.meta$cell.line.tumor.site=="BREAST",]
brca_info                         <- ccle_info[ccle_info$CCLE_Name %in% brca_meta$cell.line.name,]
brca_info                         <- subset(brca_info,lineage_molecular_subtype !="luminal_HER2_amp")
brca_info                         <- subset(brca_info,CCLE_Name!="HMEL_BREAST")
save(brca_info,file = "Pro_TNBC/paper/data/results/section_2/brca_info.RData")
GSE212143_info$celllines          <- gsub("-", "", GSE212143_info$sample_title)
GSE212143_info[8,5]               <- "HS578T"
CCLE_GSE212143                    <- intersect(GSE212143_info$celllines,ccle_info$stripped_cell_line_name)
CCLE_lest                         <- setdiff(brca_info$stripped_cell_line_name,CCLE_GSE212143)
CCLE_lest_GSE48213                <- intersect(GSE48213_info$cell_lines,CCLE_lest)
CCLE_lest_2                       <- setdiff(CCLE_lest,CCLE_lest_GSE48213)
CCLE_lest_GSE73526                <- intersect(CCLE_lest_2,GSE73526_info$sample.id)
CCLE_lest_3                       <- setdiff(CCLE_lest_2,CCLE_lest_GSE73526)
GSE212143_infor                   <- GSE212143_info[,c(2,5)]
colnames(GSE212143_infor)[2]      <- "stripped_cell_line_name"
GSE48213_infor                    <- GSE48213_info[GSE48213_info$cell_lines %in% CCLE_lest_GSE48213,c(5,2)]
GSE48213_infor                    <- GSE48213_infor[-c(14,15,23,25),]
colnames(GSE48213_infor)[2]       <- "stripped_cell_line_name"
GSE73526_infor                    <- read.delim("/home/BioklabData/lijing_Data/Cell_line/GSE73526_run_file.txt",sep = " ")
GSE73526_infor                    <- GSE73526_infor[,c(1,11)]
colnames(GSE73526_infor)[2]       <- "stripped_cell_line_name"

CCLE_test_infor                   <- bind_rows(GSE212143_infor,GSE48213_infor,GSE73526_infor)
celllines                         <- data.frame(run_accession = c("SRR19170149"),stripped_cell_line_name = c("HCC2157"))
CCLE_test_infor                   <- rbind(CCLE_test_infor,celllines)
CCLE_test_infor                   <- merge(CCLE_test_infor,brca_info[,c(3,21)],by="stripped_cell_line_name")
CCLE_test_infor                   <- arrange(CCLE_test_infor,lineage_molecular_subtype)
GSE48213_log2tpm_test             <- GSE48213_log2tpm[,colnames(GSE48213_log2tpm) %in% CCLE_test_infor$run_accession]
GSE73526_log2tpm_test             <- GSE73526_log2tpm[,colnames(GSE73526_log2tpm) %in% CCLE_test_infor$run_accession]
identical(rownames(GSE212143.log2tpm),rownames(GSE48213_log2tpm_test))
identical(rownames(GSE48213_log2tpm_test),rownames(GSE73526_log2tpm_test))
identical(rownames(GSE212143.log2tpm),rownames(GSE48213_log2tpm_test))
GSE212143.log2tpm                 <- GSE212143.log2tpm[match(rownames(GSE48213_log2tpm_test),rownames(GSE212143.log2tpm)),]
CCLE_test_log2tpm                 <- bind_cols(GSE212143.log2tpm,GSE48213_log2tpm_test,GSE73526_log2tpm_test) %>% as.matrix() 
rownames(CCLE_test_log2tpm)       <- rownames(GSE48213_log2tpm_test)
CCLE_test_log2tpm                 <- CCLE_test_log2tpm[,CCLE_test_infor$run_accession]
save(CCLE_test_log2tpm,file = "Pro_TNBC/paper/data/results/section_2/CCLE_test_log2tpm.RData")
save(CCLE_test_infor,file = "Pro_TNBC/paper/data/results/section_2/CCLE_test_infor.RData")



####*UBS93 predicting####
load("Pro_TNBC/paper/data/results/section_2/CCLE_test_log2tpm.RData")
load("Pro_TNBC/paper/data/results/section_2/CCLE_test_infor.RData")
load("Pro_TNBC/paper/data/method/UBS93.data.RData")
CCLE_brca_subtype                   <- breast.cancer.predictor(expr.of.sample = CCLE_test_log2tpm,
                                                               expr.of.centroid = UBS93.data$UBS93.centroid,
                                                               marker.gene = UBS93.data$UBS93.gene.df$ENSEMBL,
                                                               HER2.amp.signature.genes = UBS93.data$HER2.amp.signature.genes,
                                                               ssgsea.cutoff = UBS93.data$ssgsea.cutoff)
CCLE_brca.subtype                   <- CCLE_brca_subtype$subtype
CCLE_brca.subtype$run_accession     <- rownames(CCLE_brca.subtype)
CCLE_brca.subtype                   <- merge(CCLE_brca.subtype,CCLE_test_infor,by="run_accession")
colnames(CCLE_brca.subtype)[2]      <- "UBS93.subtype"
CCLE_brca.subtype                   <- CCLE_brca.subtype[,c(1:2,4,3)]
save(CCLE_brca.subtype,file = "Pro_TNBC/paper/data/results/section_2/CCLE_brca.subtype.RData")

# F1 score 
CCLE_brca.subtype[CCLE_brca.subtype$lineage_molecular_subtype=="basal_A",3] <- "Basal"
CCLE_brca.subtype[CCLE_brca.subtype$lineage_molecular_subtype=="basal_B",3] <- "Claudin_low"
CCLE_brca.subtype[CCLE_brca.subtype$lineage_molecular_subtype=="luminal",3] <- "Luminal"
CCLE_brca.subtype.F1score <- f1_fun(pre = CCLE_brca.subtype$UBS93.subtype,tru=CCLE_brca.subtype$lineage_molecular_subtype)
#Basal      Claudin_low    HER2_amp     Luminal 
#0.9655172   0.8000000    1.0000000    0.9166667 
save(CCLE_brca.subtype,file = "Pro_TNBC/paper/data/results/section_2/CCLE_brca.subtype.RData")

####*Claudin low typing for CCLE####
library(genefu)
data(claudinLowData)
# Testing Set
test                <- medianCtr(CCLE_test_log2tpm)
#Training Set
train               <- claudinLowData
train$xd            <- medianCtr(train$xd)
train_exp           <- train$xd  %>% as.data.frame()
train_exp$ENTREZID  <- rownames(train_exp)
library(tidyverse)
library(clusterProfiler)
name                  <- bitr(rownames(train_exp),fromType = 'ENTREZID',toType = 'ENSEMBL',OrgDb = 'org.Hs.eg.db')
table(rownames(test) %in% name$ENSEMBL)
gene.name             <- name[name$ENSEMBL %in% rownames(test),]
train_exp             <- merge(gene.name,train_exp,by="ENTREZID")
rown                  <- train_exp$ENSEMBL
train_exp             <- train_exp[,-c(1:2)]
rownames(train_exp)   <- rown
train_exp             <- as.matrix(train_exp)

common.gene           <- intersect(rownames(train_exp),rownames(test))
train_exp             <- train_exp[common.gene,]
test                  <- test[common.gene,]
identical(rownames(train_exp),rownames(test))

# Generate Predictions
predout                 <- claudinLow(x=train_exp, classes=as.matrix(train$classes$Group,ncol=1), y=test,distm = "spearman")
# Obtain results
results                 <- predout$predictions %>% as.data.frame()
colnames(results)[1]    <- "run_accession"
CCLE_brca.subtype       <- merge(CCLE_brca.subtype,results,by="run_accession")

CCLE.brca.subtype.compare  <- merge(CCLE.brca.subtype.compare,results,by="run_accession")
colnames(CCLE.brca.subtype.compare)[7] <- "CL.subtype.NOcenter"###centroid and sample :no center
save(CCLE.brca.subtype.compare,file = "Pro_TNBC/paper/data/results/section_3/CCLE.brca.subtype.compare.RData.RData")


####*PAM50 subtyping####
library(genefu)
library(dplyr)
data("pam50")
pam50.gene.df                  <- read.delim("Pro_TNBC/data/TCGA/TCGA/pam50_entrez_gene_id_and_ensemble_gene_id.txt", stringsAsFactors=FALSE,comment.char = '#') 
pam50.gene                     <- pam50.gene.df$ensemble.gene.id %>% as.character
colnames(pam50.gene.df)[1]     <- 'EntrezGene.ID'
colnames(pam50.gene.df)[2]     <- 'probe' 
pam50.gene.df$EntrezGene.ID    <- as.character(pam50.gene.df$EntrezGene.ID)

sampleID                       <- colnames(CCLE_test_log2tpm)
pam50.gene.expr                <- CCLE_test_log2tpm[pam50.gene.df$probe %>% as.character,sampleID] %>% t 
annot.matrix                   <- pam50.gene.df[,1:2] %>% as.matrix()
rownames(annot.matrix)         <- annot.matrix[,'probe']
pam50.subtype                  <- intrinsic.cluster.predict(sbt.model = pam50,data = pam50.gene.expr,annot = annot.matrix,mapping=annot.matrix,do.mapping = TRUE )
CCLE.pam50.subtype             <- pam50.subtype["subtype"] %>% as.data.frame()
CCLE.pam50.subtype$run_accession   <- rownames(CCLE.pam50.subtype)
CCLE_brca.subtype              <- merge(CCLE_brca.subtype,CCLE.pam50.subtype,by="run_accession")
CCLE_brca.subtype              <- CCLE_brca.subtype[,c(1:2,5:6,3:4)]
colnames(CCLE_brca.subtype)[3:4]  <- c("CL.subtype","PAM50.subtype")
CCLE.brca.subtype.compare         <- CCLE_brca.subtype
save(CCLE.brca.subtype.compare,file = "Pro_TNBC/paper/data/results/section_3/CCLE.brca.subtype.compare.RData")


CCLE.brca.subtype.compare  <- within(CCLE.brca.subtype.compare,{
  genefu.pre.subtype       <- NA
  genefu.pre.subtype[CL.subtype.NOcenter=="Claudin"] <- "Claudin_low"
  genefu.pre.subtype[CL.subtype.NOcenter=="Others"&PAM50.subtype=="Basal"]  <-"Basal"
  genefu.pre.subtype[CL.subtype.NOcenter=="Others"&PAM50.subtype=="Her2"]  <- "HER2_amp"
  genefu.pre.subtype[is.na(genefu.pre.subtype)]  <- "Luminal"
})
CCLE.brca.subtype.compare  <- within(CCLE.brca.subtype.compare,{
  genefu.pre.subtype.center       <- NA
  genefu.pre.subtype.center[CL.subtype=="Claudin"] <- "Claudin_low"
  genefu.pre.subtype.center[CL.subtype=="Others"&PAM50.subtype=="Basal"]  <-"Basal"
  genefu.pre.subtype.center[CL.subtype=="Others"&PAM50.subtype=="Her2"]  <- "HER2_amp"
  genefu.pre.subtype.center[is.na(genefu.pre.subtype.center)]  <- "Luminal"
})
CCLE.brca.subtype.compare[CCLE.brca.subtype.compare$lineage_molecular_subtype=="basal_A",5] <- "Basal"
CCLE.brca.subtype.compare[CCLE.brca.subtype.compare$lineage_molecular_subtype=="basal_B",5] <- "Claudin_low"
CCLE.brca.subtype.compare[CCLE.brca.subtype.compare$lineage_molecular_subtype=="luminal",5] <- "Luminal"
CCLE.brca.subtype.compare.F1score <- f1_fun(pre = CCLE.brca.subtype.compare$genefu.pre.subtype,tru=CCLE.brca.subtype.compare$lineage_molecular_subtype)
#Basal Claudin_low    HER2_amp     Luminal 
#0.6956522   0.8000000   0.7407407   0.8333333
CCLE.brca.subtype.compare.center.F1score <- f1_fun(pre = CCLE.brca.subtype.compare$genefu.pre.subtype.center,tru=CCLE.brca.subtype.compare$lineage_molecular_subtype)
#   Basal Claudin_low    HER2_amp     Luminal 
#0.3333333   0.5714286   0.8333333   0.8333333
save(CCLE.brca.subtype.compare,file = "Pro_TNBC/paper/data/results/section_3/CCLE.brca.subtype.compare.RData")

####2.TCGA####
####imputed subtypes for breast cancer in TCGA####
library(readxl)
load("~/Pro_TNBC/data/TCGA/TCGA/Breast Invasive Carcinoma.RData")
tcga_log2tpm               <- log2.tpm.matrix
brca_tcga_clinical_data    <- read.delim("Pro_TNBC/data/TCGA/TCGA/brca_tcga_clinical_data.tsv",sep = "\t")
brca_tcga_meta             <- subset(brca_tcga_clinical_data,Sample.ID %in% colnames(tcga_log2tpm))
MBC_meta                   <- subset(brca_tcga_meta,Cancer.Type.Detailed=="Metaplastic Breast Cancer")
TNBC_meta                  <- subset(brca_tcga_meta,ER.Status.By.IHC=="Negative"&PR.status.by.ihc=="Negative"&IHC.HER2=="Negative")
HER2_meta                  <- read.delim("Pro_TNBC/data/TCGA/TCGA/brca_HER2amp_tcga_clinical_data.tsv",sep = "\t")
cols                       <- c("Study.ID","Patient.ID","Sample.ID","Cancer.Type.Detailed","ER.Status.By.IHC","PR.status.by.ihc","IHC.HER2" ) 
brca_tcga_metadata         <- brca_tcga_meta[,cols]
brca_tcga_metadata         <- within(brca_tcga_metadata,{
  HER2_AMP <- NA
  HER2_AMP[ Sample.ID%in% HER2_meta$Sample.ID ]          <- "YES"
})
brca_tcga_metadata[is.na(brca_tcga_metadata$HER2_AMP),8] <- "NO"
brca_tcga_metadata[is.na(brca_tcga_metadata)]            <- "unknown"
brca_tcga_metadata$subtype  <- ifelse(brca_tcga_metadata$Cancer.Type.Detailed == "Metaplastic Breast Cancer", "MBC",
                                      ifelse(brca_tcga_metadata$ER.Status.By.IHC == "Negative" & brca_tcga_metadata$PR.status.by.ihc == "Negative" & brca_tcga_metadata$IHC.HER2 == "Negative", "TNBC",
                                             ifelse(brca_tcga_metadata$HER2_AMP == "YES", "HER2_AMP",
                                                    ifelse(brca_tcga_metadata$ER.Status.By.IHC == "Positive" | brca_tcga_metadata$PR.status.by.ihc == "Positive" | brca_tcga_metadata$IHC.HER2 == "Positive",
                                               "Luminal","unknown")))) 

table(brca_tcga_metadata$subtype)
save(brca_tcga_metadata,file = "Pro_TNBC/output/data/TCGA/brca_tcga_metadata.RData")
####imputed subtypes for MBC####
TCGA_MBC_log2tpm                      <- log2.tpm.matrix[,colnames(log2.tpm.matrix) %in% TCGA_MBC_subtype$Sample.ID]
TCGA_MBC_log2tpm_df                   <- as.data.frame(TCGA_MBC_log2tpm)
TCGA_MBC_log2tpm_df$ENSEMBL           <- rownames(TCGA_MBC_log2tpm_df) 
symbol                                <- bitr(TCGA_MBC_log2tpm_df$ENSEMBL,fromType = "ENSEMBL",toType = "SYMBOL",OrgDb="org.Hs.eg.db")
TCGA_MBC_log2tpm_df                   <- merge(symbol,TCGA_MBC_log2tpm_df,by="ENSEMBL")
marker.gene                           <- c("CLDN3","CLDN4","CLDN7","EPCAM","VIM")
TCGA_MBC_log2tpm_df                   <- TCGA_MBC_log2tpm_df[TCGA_MBC_log2tpm_df$SYMBOL %in% marker.gene,]
rown                                  <- TCGA_MBC_log2tpm_df$SYMBOL
TCGA_MBC_log2tpm_df                   <- TCGA_MBC_log2tpm_df[,-c(1,2)]
rownames(TCGA_MBC_log2tpm_df)         <- rown
TCGA_MBC_log2tpm_matrix               <- as.matrix(TCGA_MBC_log2tpm_df)
TCGA_MBC_log2tpm_matrix               <- TCGA_MBC_log2tpm_matrix[marker.gene,]
library(pheatmap)
fig4a1 <- pheatmap(TCGA_MBC_log2tpm_matrix,cluster_rows = F,
                   cluster_cols = T,
                   color = colorRampPalette(c("navy","white","firebrick3"))(100),
                   show_colnames = T,border_color = NA,
                   scale = "row",
                   show_rownames = T,
                   fontsize=30,cellwidth = 45, cellheight = 55
)
ggsave(fig4a1,filename = "Pro_TNBC/paper/plot/section_2/heat.map.of.marker.gene.in.TCGA.MBC.pdf",width = 20,height = 15)

CL.id                               <- c("TCGA-AC-A7VC-01","TCGA-A2-A4S1-01","TCGA-BH-A6R9-01")
colnames(brca_tcga_metadata)[9]     <- "Imputed.Subtype"
MBC.id                     <- brca_tcga_metadata[brca_tcga_metadata$Imputed.Subtype=="MBC",]$Sample.ID
LEST.id                    <- setdiff(MBC.id,CL.id)
unknown.id                 <- subset(brca_tcga_metadata,Imputed.Subtype=="unknown")$Sample.ID
ID                         <- c(LEST.id,unknown.id)
brca_tcga_metadata         <- subset(brca_tcga_metadata,!(Sample.ID  %in% ID))
brca_tcga_metadata[brca_tcga_metadata$Imputed.Subtype=="MBC",9]  <- "Claudin_low"
save(brca_tcga_metadata,file = "Pro_TNBC/paper/data/results/section_2/brca_tcga_metadata.RData")

load("Pro_TNBC/paper/data/results/section_2/brca_tcga_metadata.RData")
load("~/Pro_TNBC/data/TCGA/TCGA/Breast Invasive Carcinoma.RData")
tcga_log2tpm                        <- log2.tpm.matrix
tcga_log2tpm                        <- tcga_log2tpm[,colnames(tcga_log2tpm) %in% brca_tcga_metadata$Sample.ID]
####*UBS93 predicting####
tcga_brca_subtype                   <- breast.cancer.predictor(expr.of.sample = tcga_log2tpm,
                                                               expr.of.centroid = UBS93.data$UBS93.centroid,
                                                               marker.gene = UBS93.data$UBS93.gene.df$ENSEMBL,
                                                               HER2.amp.signature.genes = UBS93.data$HER2.amp.signature.genes,
                                                               ssgsea.cutoff = UBS93.data$ssgsea.cutoff)
tcga_brca.subtype                   <- tcga_brca_subtype$subtype
tcga_brca.subtype$Sample.ID         <- rownames(tcga_brca.subtype)
tcga_brca.subtype                   <- merge(tcga_brca.subtype,brca_tcga_metadata,by="Sample.ID")
colnames(tcga_brca.subtype)[2]      <- "UBS93.subtype"
tcga_brca.subtype                   <- tcga_brca.subtype[,c(1:2,10,3:9)]
save(tcga_brca.subtype,file = "Pro_TNBC/paper/data/results/section_2/tcga_brca.subtype.RData")

####*Claudin low typing for TCGA(center)####
library(genefu)
data(claudinLowData)
# Testing Set
test                <- medianCtr(tcga_log2tpm)
#Training Set
train               <- claudinLowData
train$xd            <- medianCtr(train$xd)
train_exp           <- train$xd  %>% as.data.frame()
train_exp$ENTREZID  <- rownames(train_exp)
library(tidyverse)
library(clusterProfiler)
name                  <- bitr(rownames(train_exp),fromType = 'ENTREZID',toType = 'ENSEMBL',OrgDb = 'org.Hs.eg.db')
table(rownames(test) %in% name$ENSEMBL)
gene.name             <- name[name$ENSEMBL %in% rownames(test),]
train_exp             <- merge(gene.name,train_exp,by="ENTREZID")
rown                  <- train_exp$ENSEMBL
train_exp             <- train_exp[,-c(1:2)]
rownames(train_exp)   <- rown
train_exp             <- as.matrix(train_exp)

common.gene           <- intersect(rownames(train_exp),rownames(test))
train_exp             <- train_exp[common.gene,]
test                  <- test[common.gene,]
identical(rownames(train_exp),rownames(test))

# Generate Predictions
predout                 <- claudinLow(x=train_exp, classes=as.matrix(train$classes$Group,ncol=1), y=test,distm = "spearman")
# Obtain results
results                 <- predout$predictions %>% as.data.frame()
colnames(results)[1]    <- "Sample.ID"
tcga_brca.subtype       <- merge(tcga_brca.subtype,results,by="Sample.ID")

save(tcga.brca.subtype.compare,file = "Pro_TNBC/paper/data/results/section_3/tcga.brca.subtype.compare.RData")


####*PAM50 subtyping####
library(genefu)
library(dplyr)
data("pam50")
pam50.centroid                 <- pam50$centroids
pam50.centroid                 <- pam50.centroid[,-5]
pam50$centroids                <- pam50.centroid
pam50.gene.df                  <- read.delim("Pro_TNBC/data/TCGA/TCGA/pam50_entrez_gene_id_and_ensemble_gene_id.txt", stringsAsFactors=FALSE,comment.char = '#') 
pam50.gene                     <- pam50.gene.df$ensemble.gene.id %>% as.character
colnames(pam50.gene.df)[1]     <- 'EntrezGene.ID'
colnames(pam50.gene.df)[2]     <- 'probe' 
pam50.gene.df$EntrezGene.ID    <- as.character(pam50.gene.df$EntrezGene.ID)

sampleID                       <- colnames(tcga_log2tpm)
pam50.gene.expr                <- tcga_log2tpm[pam50.gene.df$probe %>% as.character,sampleID] %>% t 
annot.matrix                   <- pam50.gene.df[,1:2] %>% as.matrix()
rownames(annot.matrix)         <- annot.matrix[,'probe']
pam50.subtype                  <- intrinsic.cluster.predict(sbt.model = pam50,data = pam50.gene.expr,annot = annot.matrix,mapping=annot.matrix,do.mapping = TRUE )
tcga.pam50.subtype             <- pam50.subtype["subtype"] %>% as.data.frame()
tcga.pam50.subtype$Sample.ID   <- rownames(tcga.pam50.subtype)
tcga_brca.subtype              <- merge(tcga_brca.subtype,tcga.pam50.subtype,by="Sample.ID")
tcga_brca.subtype              <- tcga_brca.subtype[,c(1:3,11:12,4:10)]
colnames(tcga_brca.subtype)[4:5]  <- c("CL.subtype","PAM50.subtype")
tcga.brca.subtype.compare         <- tcga_brca.subtype
save(tcga.brca.subtype.compare,file = "Pro_TNBC/paper/data/results/section_3/tcga.brca.subtype.compare.RData")

####*Claudin low typing for TCGA(no center)####
library(genefu)
data(claudinLowData)
# Testing Set
test                <- tcga_log2tpm
#Training Set
train               <- claudinLowData
train_exp           <- train$xd  %>% as.data.frame()
train_exp$ENTREZID  <- rownames(train_exp)
library(tidyverse)
library(clusterProfiler)
name                  <- bitr(rownames(train_exp),fromType = 'ENTREZID',toType = 'ENSEMBL',OrgDb = 'org.Hs.eg.db')
table(rownames(test) %in% name$ENSEMBL)
gene.name             <- name[name$ENSEMBL %in% rownames(test),]
train_exp             <- merge(gene.name,train_exp,by="ENTREZID")
rown                  <- train_exp$ENSEMBL
train_exp             <- train_exp[,-c(1:2)]
rownames(train_exp)   <- rown
train_exp             <- as.matrix(train_exp)

common.gene           <- intersect(rownames(train_exp),rownames(test))
train_exp             <- train_exp[common.gene,]
test                  <- test[common.gene,]
identical(rownames(train_exp),rownames(test))

# Generate Predictions
predout                 <- claudinLow(x=train_exp, classes=as.matrix(train$classes$Group,ncol=1), y=test,distm = "spearman")
# Obtain results
results                 <- predout$predictions %>% as.data.frame()
colnames(results)[1]    <- "Sample.ID"
tcga.brca.subtype.compare  <- merge(tcga.brca.subtype.compare,results,by="Sample.ID")
colnames(tcga.brca.subtype.compare)[13] <- "CL.subtype.NOcenter" ###centroid and sample :no center
colnames(tcga.brca.subtype.compare)[4]  <- "CL.subtype.center"
save(tcga.brca.subtype.compare,file = "Pro_TNBC/paper/data/results/section_3/tcga.brca.subtype.compare.RData")

#F1 score
tcga.brca.subtype.compare  <- within(tcga.brca.subtype.compare,{
  genefu.pre.subtype       <- NA
  genefu.pre.subtype[CL.subtype.NOcenter=="Claudin"] <- "Claudin_low"
  genefu.pre.subtype[CL.subtype.NOcenter=="Others"&PAM50.subtype=="Basal"]  <-"Basal"
  genefu.pre.subtype[CL.subtype.NOcenter=="Others"&PAM50.subtype=="Her2"]  <- "HER2_amp"
  genefu.pre.subtype[is.na(genefu.pre.subtype)]  <- "Luminal"
})
tcga.brca.subtype.compare  <- within(tcga.brca.subtype.compare,{
  genefu.pre.subtype.center       <- NA
  genefu.pre.subtype.center[CL.subtype.center=="Claudin"] <- "Claudin_low"
  genefu.pre.subtype.center[CL.subtype.center=="Others"&PAM50.subtype=="Basal"]  <-"Basal"
  genefu.pre.subtype.center[CL.subtype.center=="Others"&PAM50.subtype=="Her2"]  <- "HER2_amp"
  genefu.pre.subtype.center[CL.subtype.center=="Others"&PAM50.subtype=="Normal"]  <- "Normal"
  genefu.pre.subtype.center[is.na(genefu.pre.subtype.center)]  <- "Luminal"
})
tcga.brca.subtype.compare[tcga.brca.subtype.compare$Imputed.Subtype=="TNBC",3] <- "Basal"
tcga.brca.subtype.compare[tcga.brca.subtype.compare$Imputed.Subtype=="HER2_AMP",3] <- "HER2_amp"
tcga.brca.subtype.compare.F1score <- f1_fun(pre = tcga.brca.subtype.compare$genefu.pre.subtype,tru=tcga.brca.subtype.compare$Imputed.Subtype)
#Basal       Claudin_low    HER2_amp     Luminal 
# 0.6990291    0.1621622   0.4812834   0.9156785 
tcga.brca.subtype.compare.center.F1score <- f1_fun(pre = tcga.brca.subtype.compare$genefu.pre.subtype.center,tru=tcga.brca.subtype.compare$Imputed.Subtype)
#   Basal    Claudin_low    HER2_amp     Luminal 
#0.10169492    0.01239669  0.32051282    0.67394958  

tcga.brca.subtype.compare.UBS93.F1score  <- f1_fun(pre = tcga.brca.subtype.compare$UBS93.subtype,tru = tcga.brca.subtype.compare$Imputed.Subtype)
#Basal Claudin_low    HER2_amp     Luminal 
#0.7755102   0.6666667   0.5126354   0.8873239 
save(tcga.brca.subtype.compare,file = "Pro_TNBC/paper/data/results/section_3/tcga.brca.subtype.compare.RData")
genefu_confusion <- table(tcga.brca.subtype.compare$Imputed.Subtype,tcga.brca.subtype.compare$genefu.pre.subtype)
UBS93_confusion  <- table(tcga.brca.subtype.compare$Imputed.Subtype,tcga.brca.subtype.compare$UBS93.subtype)
F1.score.compare <- data.frame(F1.score=c(tcga.brca.subtype.compare.F1score,tcga.brca.subtype.compare.UBS93.F1score),gene_panel=c(rep("genefu",4),rep("UBS93",4)),subtype=c(rep(c("Basal-like","Claudin-low","HER2-amp","Luminal"),2)))
save(F1.score.compare,file = "Pro_TNBC/paper/data/results/section_3/F1.score.compare.RData")


####3. mouse model####
load("Pro_TNBC/data/mouse/mouse_model/gene.expression.RData")
GSE157333_log2tpm             <- log2.tpm.matrix
load("~/Pro_TNBC/data/mouse/mouse_claudin_low_ref/mouse_claudin_low_ref.RData")
CL_ref_log2tpm                    <- log2.tpm.matrix
identical(rownames(GSE157333_log2tpm),rownames(CL_ref_log2tpm))
mouse_log2tpm                     <- cbind(GSE157333_log2tpm,CL_ref_log2tpm)
load("Pro_TNBC/paper/data/results/section_2/mouse_log2tpm_info.RData")
####*UBS93####
gene_mart                                             <- read.delim("Pro_TNBC/data/mouse/mart_export.txt",sep = ",")
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

####*Claudin-low typing(center)####
library(genefu)
data(claudinLowData)

#Training Set
train               <- claudinLowData
train$xd            <- medianCtr(train$xd)
train_exp           <- train$xd  %>% as.data.frame()
train_exp$ENTREZID  <- rownames(train_exp)
library(tidyverse)
library(clusterProfiler)
name                  <- bitr(rownames(train_exp),fromType = 'ENTREZID',toType = 'ENSEMBL',OrgDb = 'org.Hs.eg.db')
table(rownames(test) %in% name$ENSEMBL)
train_exp             <- merge(name,train_exp,by="ENTREZID")
rown                  <- train_exp$ENSEMBL
train_exp             <- train_exp[,-c(1:2)]
rownames(train_exp)   <- rown
train_exp             <- as.matrix(train_exp)
##test
mouse_gene                 <- gene_mart[gene_mart$Mouse.gene.stable.ID %in% rownames(mouse_log2tpm),]
test                       <- mouse_log2tpm[rownames(mouse_log2tpm) %in% mouse_gene$Mouse.gene.stable.ID,] %>% as.data.frame()
test$Mouse.gene.stable.ID  <- rownames(test)
test                       <- merge(mouse_gene,test,by="Mouse.gene.stable.ID")
rown                       <- test$Gene.stable.ID
test                       <- test[,-(1:5)]
# Find duplicate row names
duplicated_rows             <- duplicated(rown)
duplicated_indices          <- which(duplicated_rows)
#Output index of duplicate rows Output index of duplicate rows
print(duplicated_indices)
#Delete duplicate row names 
test                        <- test[-duplicated_indices, ]
rown                        <- rown[-duplicated_indices]
rownames(test)              <- rown
test                        <- as.matrix(test)
test                        <- medianCtr(test)

common.gene           <- intersect(rownames(train_exp),rownames(test))
train_exp             <- train_exp[common.gene,]
test                  <- test[common.gene,]
identical(rownames(train_exp),rownames(test))

# Generate Predictions
predout                     <- claudinLow(x=train_exp, classes=as.matrix(train$classes$Group,ncol=1), y=test,distm = "spearman")
# Obtain results
results                     <- cbind(predout$predictions, predout$distances)
colnames(results)[1]        <- "run_accession"
mouse.subtype               <- merge(mouse.subtype,results[,1:2],by="run_accession")
colnames(mouse.subtype)[4]  <- "CL.center.subtype"       

####*PAM50####
library(genefu)
library(dplyr)
data("pam50")
pam50.gene.df                                         <- read.delim("Pro_TNBC/data/TCGA/TCGA/pam50_entrez_gene_id_and_ensemble_gene_id.txt", stringsAsFactors=FALSE,comment.char = '#') 
PAM50.mouse.gene                                      <- gene_mart[gene_mart$Gene.stable.ID %in% pam50.gene.df$ensemble.gene.id,]
PAM50.mouse.gene                                      <- PAM50.mouse.gene[!is.na(PAM50.mouse.gene$Mouse.gene.stable.ID),]
Mouse_log2tpm_PAM50                                   <- mouse_log2tpm[rownames(mouse_log2tpm) %in% PAM50.mouse.gene$Mouse.gene.stable.ID,]
Mouse_log2tpm_PAM50                                   <- as.data.frame(Mouse_log2tpm_PAM50)
Mouse_log2tpm_PAM50$Mouse.gene.stable.ID              <- rownames(Mouse_log2tpm_PAM50)
Mouse_log2tpm_PAM50                                   <- merge(PAM50.mouse.gene,Mouse_log2tpm_PAM50,by="Mouse.gene.stable.ID")
rown                                                  <- Mouse_log2tpm_PAM50$Gene.stable.ID
Mouse_log2tpm_PAM50                                   <- Mouse_log2tpm_PAM50[,-(1:5)] %>% as.matrix()
rownames(Mouse_log2tpm_PAM50)                         <- rown
pam50.gene                                            <- pam50.gene.df$ensemble.gene.id %>% as.character#将数字转化成字符,Ensemble ID形式：ENSG00000223972
colnames(pam50.gene.df)[1]                            <- 'EntrezGene.ID'
colnames(pam50.gene.df)[2]                            <- 'probe' 
pam50.gene.df$EntrezGene.ID                           <- as.character(pam50.gene.df$EntrezGene.ID)

sampleID                                              <- colnames(Mouse_log2tpm_PAM50)
Mouse_log2tpm_PAM50                                   <- Mouse_log2tpm_PAM50[pam50.gene.df$probe %>% as.character,sampleID] %>% t #t函数为转置函数，可以将行转成列，将列转成行#
annot.matrix                                          <- pam50.gene.df[,1:2] %>% as.matrix#将pam50.gene.df中的第一和二列转化成矩阵
rownames(annot.matrix)                                <- annot.matrix[,'probe']
Mouse.pam50.subtype                                   <- intrinsic.cluster.predict(sbt.model = pam50,data = Mouse_log2tpm_PAM50,annot = annot.matrix,mapping=annot.matrix,do.mapping = T)
table(Mouse.pam50.subtype$subtype)
Mouse_pam50_subtype                                   <- Mouse.pam50.subtype$subtype %>% as.data.frame()
Mouse_pam50_subtype$run_accession                     <- rownames(Mouse_pam50_subtype)
colnames(Mouse_pam50_subtype)[1]                      <- "PAM50.subtype"
mouse.subtype                                         <- merge(Mouse_pam50_subtype,mouse.subtype,by="run_accession")
colnames(mouse.subtype)[4]                            <- "Imputed.subtype"

####*Claudin-low typing(no center)####
library(genefu)
data(claudinLowData)

#Training Set
train               <- claudinLowData
train_exp           <- train$xd  %>% as.data.frame()
train_exp$ENTREZID  <- rownames(train_exp)
library(tidyverse)
library(clusterProfiler)
name                  <- bitr(rownames(train_exp),fromType = 'ENTREZID',toType = 'ENSEMBL',OrgDb = 'org.Hs.eg.db')
table(rownames(test) %in% name$ENSEMBL)
gene.name             <- name[name$ENSEMBL %in% rownames(test),]
train_exp             <- merge(gene.name,train_exp,by="ENTREZID")
rown                  <- train_exp$ENSEMBL
train_exp             <- train_exp[,-c(1:2)]
rownames(train_exp)   <- rown
train_exp             <- as.matrix(train_exp)
##test
mouse_gene                 <- gene_mart[gene_mart$Mouse.gene.stable.ID %in% rownames(mouse_log2tpm),]
test                       <- mouse_log2tpm[rownames(mouse_log2tpm) %in% mouse_gene$Mouse.gene.stable.ID,] %>% as.data.frame()
test$Mouse.gene.stable.ID  <- rownames(test)
test                       <- merge(mouse_gene,test,by="Mouse.gene.stable.ID")
rown                       <- test$Gene.stable.ID
test                       <- test[,-(1:5)]
# Find duplicate row names
duplicated_rows             <- duplicated(rown)
duplicated_indices          <- which(duplicated_rows)
#Output index of duplicate rows Output index of duplicate rows
print(duplicated_indices)
#Delete duplicate row names 
test                        <- test[-duplicated_indices, ]
rown                        <- rown[-duplicated_indices]
rownames(test)              <- rown
test                        <- as.matrix(test)


common.gene           <- intersect(rownames(train_exp),rownames(test))
train_exp             <- train_exp[common.gene,]
test                  <- test[common.gene,]
identical(rownames(train_exp),rownames(test))

# Generate Predictions
predout                     <- claudinLow(x=train_exp, classes=as.matrix(train$classes$Group,ncol=1), y=test,distm = "spearman")
# Obtain results
results                     <- cbind(predout$predictions, predout$distances)
colnames(results)[1]        <- "run_accession"
mouse.subtype               <- merge(mouse.subtype,results[,1:2],by="run_accession")
colnames(mouse.subtype)[6]  <- "CL.nocenter.subtype"   
mouse.subtype.compare       <- mouse.subtype  
save(mouse.subtype.compare,file = "Pro_TNBC/paper/data/results/section_3/mouse.subtype.compare.RData")


####4.Cell lines from other papers(single-cell)####
######*predicting the single cell RNAseq of cell lines(GSE173634) by UBS93 ####
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

#########*predicting the scRNAseq data of cell lines (GSE202771) by UBS93####
library(GSVA)
table(GSE202771_scRNA$orig.ident)
CCLE.Name                       <- names(table(GSE202771_scRNA$orig.ident))
GSE202771_info                  <- data.frame(sample.id=CCLE.Name)
GSE202771_info$CCLE.Name        <- apply(GSE202771_info[,1,drop=F],1,function(x){strsplit(x, "_")[[1]][1]})
GSE202771_info$CCLE_Name        <- paste(GSE202771_info$CCLE.Name,"_BREAST",sep ="")
GSE202771_infor                 <- left_join(GSE202771_info,ccle_info[,c(4,21)],by="CCLE_Name")
GSE202771_infor[c(13,16,17,21:23),4]           <- c("basal_A","basal_B","basal_B","basal_A","luminal","basal_A")
GSE202771_infor[c(13,21:23),4]           <- c(NA,"basal","luminal","basal")
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
colnames(GSE202771_subtype_result)[1] <- "lineage_molecular_subtype"
GSE202771_subtype_basal                 <- subset(GSE202771_subtype_result,lineage_molecular_subtype=="basal_A")
GSE202771_subtype_basal$ratio           <- GSE202771_subtype_basal$Basal/GSE202771_subtype_basal$total_numbers
GSE202771_subtype_cl                    <- subset(GSE202771_subtype_result,lineage_molecular_subtype=="basal_B")
GSE202771_subtype_cl$ratio              <- GSE202771_subtype_cl$Claudin_low/GSE202771_subtype_cl$total_numbers
GSE202771_subtype_luminal                  <- subset(GSE202771_subtype_result,lineage_molecular_subtype=="luminal")
GSE202771_subtype_luminal$ratio            <- GSE202771_subtype_luminal$Luminal/GSE202771_subtype_luminal$total_numbers
GSE202771_subtype                       <- rbind(GSE202771_subtype_basal,GSE202771_subtype_cl)
GSE202771_subtype                       <- rbind(GSE202771_subtype,GSE202771_subtype_luminal)
GSE202771_subtype$celllines_name        <- rownames(GSE202771_subtype)
GSE202771_subtype[GSE202771_subtype$lineage_molecular_subtype=="basal_A",1]  <- c(rep("Basal",7))
GSE202771_subtype[GSE202771_subtype$lineage_molecular_subtype=="basal_B",1]  <- c(rep("Claudin_low",6))
GSE202771_subtype[GSE202771_subtype$lineage_molecular_subtype=="luminal",1]  <- c(rep("Luminal",2))
save(GSE202771_subtype,file = "Pro_TNBC/paper/data/results/section_2/supfig_4.GSE202771.SUBTYPE.RData")



####*assemble####
library(Seurat)
load("Pro_TNBC/output/data/CCLE/GSE173634_scRNA.RData")
CCLE.Name                       <- names(table(GSE173634_scRNA$orig.ident))
load("Pro_TNBC/paper/data/results/section_2/brca_info.RData")
brca.id            <- intersect(CCLE.Name,brca_info$stripped_cell_line_name)
brca.lest.id       <- setdiff(brca_info$stripped_cell_line_name,brca.id)

load("Pro_TNBC/output/data/scRNASeq/GSE202771/GSE202771_scRNA.RData")
GSE202771.id       <- names(table(GSE202771_scRNA$orig.ident))
GSE202771.id       <- strsplit(GSE202771.id,"[_]")
GSE202771.id       <- sapply(GSE202771.id,function(x)x[1])
GSE202771.id.CCLE  <- intersect(brca.lest.id,GSE202771.id)
GSE202771.id.LEST  <- setdiff(brca.lest.id,GSE202771.id.CCLE)

load("~/Pro_TNBC/paper/data/results/section_2/fig2b.GSE17634.SUBTYPE.RData")
load("~/Pro_TNBC/paper/data/results/section_2/supfig_4.GSE202771.SUBTYPE.RData")
GSE173634_subtype           <- GSE173634_subtype[GSE173634_subtype$celllines_name %in% brca.id,]
GSE173634_subtype[is.na(GSE173634_subtype)]  <- 0
GSE202771_subtype[is.na(GSE202771_subtype)]  <- 0
GSE202771_subtype$CCLE_name <- apply(GSE202771_subtype,1,function(x){strsplit(x[9],"[_]")[[1]][1]})
GSE202771_subtype           <- GSE202771_subtype[GSE202771_subtype$CCLE_name %in% GSE202771.id.CCLE,]
GSE202771_subtype           <- GSE202771_subtype[c(1,3),]
GSE202771_subtype           <- GSE202771_subtype[,c(1,8,9)]
GSE173634_subtype           <- GSE173634_subtype[,c(1,7,8)]
scRNAseq_UBS93_subtype      <- rbind(GSE173634_subtype,GSE202771_subtype)
save(scRNAseq_UBS93_subtype,file = "Pro_TNBC/paper/data/results/section_2/scRNAseq_UBS93_subtype.RData")



load("~/Pro_TNBC/paper/data/results/section_3/GSE173634_subtype_compare.RData")
load("~/Pro_TNBC/paper/data/results/section_3/GSE202771_subtype_compare.RData")
GSE173634_subtype_compare  <- GSE173634_subtype_compare[GSE173634_subtype_compare$celllines_name %in% brca.id,]
GSE202771_subtype_compare  <- GSE202771_subtype_compare[GSE202771_subtype_compare$celllines_name %in% c("HCC1806_B2_scRNA"),]
scRNAseq_UBS93_subtype_compare      <- rbind(GSE173634_subtype_compare,GSE202771_subtype_compare)
save(scRNAseq_UBS93_subtype_compare,file = "Pro_TNBC/paper/data/results/section_3/scRNAseq_UBS93_subtype_compare.RData")
