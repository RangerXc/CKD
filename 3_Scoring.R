
library(Seurat)
library(dplyr)
library(monocle)
library(ggsci)

###
gene_list <- read.csv("Core ECM score.csv", header = F)
gene_list <- gene_list$V1
gene_list[1] <- "COL1A1"
geneList <- list()
geneList <- c(geneList, list(gene_list))

###AddModuleScore
Subgroup <- AddModuleScore(Subgroup, features = geneList) # 结果在matadata的cluster1属性里面
colnames(Subgroup@meta.data)[1] <- "ECM-hi Fibs"
Subgroup@meta.data[1:5,]
FeaturePlot(Subgroup, features = c("Cluster1"), cols = c("#FAEAD3", "#FF0000")) #color

###Scoring
Subgroup$cluster<-Idents(Subgroup)###将细胞命加入对象
Subgroup@meta.data[1:5,]
cell_score <- Subgroup@meta.data[which(Subgroup@meta.data$cluster %in% c("ecm-hi Fibs", "ecm-med Fibs", "CCL1921-hi Fibs", "CCL2-hi Fibs","Actin Fibs","Metabol-Fibs","Fibs")),c(8,9)]
colnames(cell_score)[1] <- "score"

library(ggpubr)
ggviolin(cell_score, x = "cluster", y = "score", fill = "cluster",
         add = "boxplot",
         add.params = list(width = 0.1,
                           fill = "white"),width = 1.1) +
  labs(x = "", y = "ECM-hi Fibs", fill = "Cell Types") +
  stat_compare_means(method = "anova")  + scale_fill_manual(values = mycolor) 