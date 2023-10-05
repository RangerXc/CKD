
library(Seurat)
library(dplyr)
library(dplyr)
library(msigdbr)
library(clusterProfiler)
library(devtools)
library(harmony)

###Fibroblasts extraction
Fibroblasts<-CKD_val[,Idents(CKD_val)%in%c("Fibroblasts")]

###Saverds
saveRDS(Fibroblasts, file = "Fibroblasts.rds")

###Reclustering
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
Fibroblasts <- CellCycleScoring(Fibroblasts, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
Fibroblasts@meta.data[1:5,]

Fibroblasts<-NormalizeData(Fibroblasts,verbose = T) 
Fibroblasts<-FindVariableFeatures(Fibroblasts,selection.method = "vst", nfeatures = 2000)
Fibroblasts<-ScaleData(Fibroblasts,vars.to.regress = c("percent.mt","S.Score","G2M.Score"),verbose = FALSE)
Fibroblasts<-RunPCA(Fibroblasts,verbose = T,npcs = 50)
ElbowPlot(Fibroblasts,ndims = 100)
p1 <- DimPlot(object = Fibroblasts, reduction = "pca", pt.size = .1, group.by = "orig.ident")
p2 <- VlnPlot(object = Fibroblasts, features = "PC_1", group.by = "orig.ident", pt.size = .1)
CombinePlots(plots=list(p1,p2))
Fibroblasts<-RunHarmony(Fibroblasts,"orig.ident", plot_convergence = TRUE)

harmony_embeddings <- Embeddings(Fibroblasts, 'harmony')
dim(harmony_embeddings)
p3 <- DimPlot(object = Fibroblasts, reduction = "harmony", pt.size = .1, group.by = "orig.ident")
p4 <- VlnPlot(object = Fibroblasts, features = "harmony_1", group.by = "orig.ident", pt.size = .1)
CombinePlots(plots=list(p3,p4))

p1+p3

Fibroblasts <- Fibroblasts %>% 
  RunUMAP(reduction = "harmony", dims = 1:50) %>% 
  RunTSNE(reduction = "harmony", dims = 1:50) %>%
  FindNeighbors(reduction = "harmony", dims = 1:50)


Fibroblasts<-FindClusters(Fibroblasts,resolution = 1)

####vlnplot and dotplot
features = c("POSTN","VCAN","FN1","COL1A1","COMP","COL3A1","CCL19","CCL21","CXCL12","CCL2","CD74","HLA-DRB1","AIF","CLIC5","PODXL","IGF1","IGF2","IGFBP3")
length(features)
VlnPlot(Fibroblasts,features = features[1:3],cols = cors,pt.size = 0)
VlnPlot(Fibroblasts,features = features[4:6],cols = cors,pt.size = 0)
VlnPlot(Fibroblasts,features = features[7:9],cols = cors,pt.size = 0)
VlnPlot(Fibroblasts,features = features[10:12],cols = cors,pt.size = 0)
VlnPlot(Fibroblasts,features = features[13:15],cols = cors,pt.size = 0)
VlnPlot(Fibroblasts,features = features[16:18],cols = cors,pt.size = 0)

###Bubble plot
features = c("FN1","VCAN","COL1A1","COL3A1","POSTN","BGN","COMP","TIMP1","C7","CCL19","CCL21","CCL2","CXCL12")
DotPlot(Subgroup,features = features,dot.min=0,col.max = 2,col.min = -1,cols = c("grey","#D53E4F"))
###Color。
brewer.pal(11,'Spectral')
display.brewer.all()
mycolor<-colorRampPalette(brewer.pal(11,'Spectral'))(2)

###Doheatmap
features<-unique(markers$gene[markers$p_val_adj<0.001&markers$avg_log2FC>1.2]) 
features
Fibroblasts$cluster<-Idents(Fibroblasts)
Fibroblast<- ScaleData(Fibroblasts, vars.to.regress = c("percent.mt","S.Score","G2M.Score"),features = features)
top10 <- markers%>%group_by(cluster) %>%top_n(n=10,wt=avg_log2FC)
Fibroblasts$cluster <- factor(x = Fibroblasts$cluster, levels = c("ecm-hi Fibs", "ecm-med Fibs", "CCL1921-hi Fibs", "CCL2-hi Fibs","Actin Fibs","Metabol-Fibs","Fibs"))
DoHeatmap(Fibroblasts,
          features = as.character(unique(top10$gene)),
          group.by = "cluster",
          assay = "RNA",
          group.colors = cors)+ 
  scale_fill_gradientn(colors = c("navy","white","firebrick3"))

###Volcano plot
library(ggplot2)
library(tidyverse)
library(ggrepel)

###Input
DEGs <- read.csv("Fibroblasts_markers.csv",header = T)
head(DEGs)
#Subset_1；
top10sig0 <- filter(DEGs,cluster=="0") %>% distinct(geneID,.keep_all = T) %>% top_n(20,abs(log2FC))
head(top10sig0)
#Subset_2；
top10sig1 <- filter(DEGs,cluster=="1") %>% distinct(geneID,.keep_all = T) %>% top_n(20,abs(log2FC))
head(top10sig1)
#Subset_3；
top10sig2 <- filter(DEGs,cluster=="2") %>% distinct(geneID,.keep_all = T) %>% top_n(20,abs(log2FC))
head(top10sig2)
#Subset_4；
top10sig3 <- filter(DEGs,cluster=="3") %>% distinct(geneID,.keep_all = T) %>% top_n(20,abs(log2FC))
head(top10sig3)
#Subset_5；
top10sig4 <- filter(DEGs,cluster=="4") %>% distinct(geneID,.keep_all = T) %>% top_n(20,abs(log2FC))
head(top10sig4)
#Subset_6；
top10sig5 <- filter(DEGs,cluster=="5") %>% distinct(geneID,.keep_all = T) %>% top_n(20,abs(log2FC))
head(top10sig5)
#Subset_7；
top10sig6 <- filter(DEGs,cluster=="6") %>% distinct(geneID,.keep_all = T) %>% top_n(20,abs(log2FC))
head(top10sig6)
###Combination：
top10sig <- rbind(top10sig0,top10sig1,top10sig2,top10sig3,top10sig4,top10sig5,top10sig6)
###
DEGs$size <- case_when(!(DEGs$geneID %in% top10sig$geneID)~ 1,
                     DEGs$geneID %in% top10sig$geneID ~ 2)
###
dt <- filter(DEGs,size==1)
head(dt)

#绘制每个Cluster Top10以外基因的散点火山图：
p <- ggplot()+
  geom_jitter(data = dt,
              aes(x = cluster, y = log2FC, color = "red"),
              size = 0.85,
              width =0.4)
p
###
p <- ggplot()+
  geom_jitter(data = dt,
              aes(x = cluster, y = log2FC,color="dt"),
              size = 0.85,###散点大小
              width =0.4)+###每个cluster宽度
  geom_jitter(data = top10sig,
              aes(x = cluster, y = log2FC,color="top10sig"),
              size = 1,
              width =0.4)
p
###：
DEGsbar<-data.frame(x=c(0,1,2,3,4,5,6),
                  y=c(2,2.6,2.8,2.5,5.4,3.2,2.1))
DEGsbar1<-data.frame(x=c(0,1,2,3,4,5,6),
                   y=c(-3.3,-3.5,-2.2,-3.7,-5.8,-3.1,-2))
###
p1 <- ggplot()+
  geom_col(data = DEGsbar,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)+
  geom_col(data = DEGsbar1,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)
p1
###
p2 <- ggplot()+
  geom_col(data = DEGsbar,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)+
  geom_col(data = DEGsbar1,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)+
  geom_jitter(data = dt,
              aes(x = cluster, y = log2FC,color = "dt",alpha=0.4),
              size = 1.5,
              width =0.4)+
  geom_jitter(data = top10sig,
              aes(x = cluster, y = log2FC,color = "top10sig",alpha=0.4),
              size = 2.5,
              width =0.4)
p2
###
library(RColorBrewer)
mycolor<-colorRampPalette(brewer.pal(9,'Set3'))(8)

DEGscol<-data.frame(x=c(0:6),
                  y=0,
                  label=c(0:6))
mycol <- cors
p3 <- p2 + geom_tile(data = DEGscol,
                     aes(x=x,y=y),
                     height=0.4,
                     color = "black",
                     fill = mycol,
                     alpha = 0.6,
                     show.legend = F)
p3
###
p4 <- p3+
  geom_text_repel(
    data=top10sig,
    aes(x=cluster,y=log2FC,label=geneID,
    ),size = 3.5,###字体大小
    force = 1.2,color = "black",segment.color = "white",
    show.legend = FALSE ,point.size = 5,
    arrow = NULL
  )
p4
###
library(RColorBrewer)
display.brewer.all()
brewer.pal(9,"Set1")
p5 <- p4 +
  scale_color_manual(values=c("#377EB8","#4DAF4A","#FC4E2A"))###散点颜色
p5
###：
p6 <- p5+
  theme_minimal()+
  theme(
    axis.title = element_text(size = 1,
                              color = "black",
                              face = "bold"),
    axis.line.y = element_line(color = "black",
                               size = 1.2),
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "top",
    legend.direction = "vertical",
    legend.justification = c(1,0),
    legend.text = element_text(size = 0.1)
  )
p6



