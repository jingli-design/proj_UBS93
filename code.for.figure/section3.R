#################################################
####Figure 3: Comparison to existing methods.####
#################################################

####Figure 3a: Performance comparison between genefu and UBS93 classifier on TCGA breast cancer samples.####
load("Pro_TNBC/paper/data/results/section_3/F1.score.compare.RData")
fig_3                   <- ggplot(F1.score.compare, aes(x = gene_panel, y = F1.score))+
  geom_point(size = 6,color='#6600FF') + 
  geom_point(size = 6,shape=21) + 
  geom_line(aes(group = subtype), color = 'black', lwd = 0.5) +
  labs(x = '', y = 'F1 score')+ggplot.style

ggsave(fig_3,filename = "Pro_TNBC/paper/plot/section_3/F1.score.compear.pdf",height = 15,width = 20)


####Figure 3b: Performance comparison between SCSubtype and UBS93 classifier on an assembled scRNA-seq dataset. ####
library(ggpubr)
fig_3e     <- ggboxplot(scRNAseq_UBS93_subtype_compare,x = "Subtype",y = "ratio",color = "classify",add = "point",size= 1,add.params=list(size=3))+ 
  ggplot.style+
  scale_color_manual(values =c(SCSubtype="Blue",UBS93="Red"))+
  theme(panel.grid.major=element_line(colour=NA),
        panel.grid.minor = element_blank())+xlab("Subtype")+ylab("Accuracy")+
  scale_x_discrete("Subtype", labels = c("Basal" = paste("Basal-like",sep = "\n","(n=12)"),"Luminal" = paste("Luminal",sep = "\n","(n=8)"),"HER2_amp"=paste0("HER2-amp",sep="\n","(n=5)")))
ggsave(fig_3e,filename = "Pro_TNBC/paper/plot/section_3/fig_3e.scRNA.compare.subtype.pdf",width = 20,height = 15)

scRNAseq_UBS93_subtype_compare           <- arrange(scRNAseq_UBS93_subtype_compare,Subtype)
scRNAseq_UBS93_subtype_compare_basal     <- subset(scRNAseq_UBS93_subtype_compare,Subtype=="Basal")
wilcox.test(ratio~classify,scRNAseq_UBS93_subtype_compare_basal,paired=T)#p-value = 0.004883
scRNAseq_UBS93_subtype_compare_luminal   <- subset(scRNAseq_UBS93_subtype_compare,Subtype=="Luminal")
wilcox.test(ratio~classify,scRNAseq_UBS93_subtype_compare_luminal,paired=T)#p-value = 0.3828
scRNAseq_UBS93_subtype_compare_HER2      <- subset(scRNAseq_UBS93_subtype_compare,Subtype=="HER2_amp")
wilcox.test(ratio~classify,scRNAseq_UBS93_subtype_compare_HER2,paired=T)#p-value = 0.125
scRNAseq_UBS93_subtype_compare$Subtype   <- factor(scRNAseq_UBS93_subtype_compare$Subtype,levels =c("Luminal","HER2_amp","Basal"))


#### fig S2b: Box plot of the mean score of UBS93 gene panel and PAM50 gene panel.####
library(readr)
library(ggplot2)
load("Pro_TNBC/output/data/scRNASeq/26_sample/pam50.gene/PAM50.gene.name.RData")
load("~/Pro_TNBC/output/data/scRNASeq/26_sample/pam50.gene/PAM50.gene.name.RData")
mean_score         <- read_csv("Pro_TNBC/paper/data/results/section_1/mean_score.csv")
UBS93_score        <- mean_score[mean_score$gene %in% UBS93.gene.df$SYMBOL,]
PAM50_score        <- mean_score[mean_score$gene %in% pam50.gene,]
UBS93_pam50_score  <- data.frame(s_mean=c(UBS93_score$S_mean,PAM50_score$S_mean),group=c(rep("UBS93_score",93),rep("pam50_gene",50)))
UBS93_pam50_score$group  <- factor(UBS93_pam50_score$group,levels = c("pam50_gene","UBS93_score"))
fig_3a <- ggboxplot(UBS93_pam50_score ,x = "group",y = "s_mean",add = "point",size= 1,add.params=list(size=5))+ 
  ggplot.style+
  theme(panel.grid.major=element_line(colour=NA),
        panel.grid.minor = element_blank())+ 
  xlab('') + ylab('C-M-ratio')+
  scale_x_discrete("", labels = c("UBS93_score" = paste("UBS93",sep = "\n","(n=93)"),"pam50_gene" = paste("PAM50",sep = "\n","(n=50)")))+
  stat_compare_means(method = "wilcox.test")
ggsave(fig_3a,filename = "Pro_TNBC/paper/plot/section_3/the.score.of.UBS93.gene.panel.and.pam50.gene.panel.pdf",width=20,height=15)
wilcox.test(s_mean~group,data = UBS93_pam50_score)#p-value < 2.2e-16
save(UBS93_pam50_score,file="Pro_TNBC/paper/data/section_3/fig3a.UBS93.pam50.score.RData")











