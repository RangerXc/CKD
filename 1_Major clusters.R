
library(Seurat)
library(dplyr)
library(monocle)
library(ggsci)

####CKD
CKD<-Read10X(data.dir = "./ScRNA_CKD/")
CKD[1:5,1:5]

###cell quality control
CKD_object<- CreateSeuratObject(counts = CKD, project = "CKD", min.cells = 0, min.features = 200)
CKD_object[["percent.mt"]] <- PercentageFeatureSet(CKD_object, pattern = "^MT-")
hist(CKD_object[["percent.mt"]]$percent.mt)
VlnPlot(CKD_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(CKD_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(CKD_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1,plot2))
CKD_object

###select the cells with 300 genes at least and 4000 at most, the percent of mitochondrion genes is less than 10%
CKD_val<- subset(CKD_object, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 20)
CKD_val
CKD_val@meta.data[1:5,]

###Normalization
CKD_val <- NormalizeData(CKD_val, normalization.method = "LogNormalize", scale.factor = 10000)
CKD_val <- FindVariableFeatures(CKD_val, selection.method = "vst", nfeatures = 2000)
CKD_val@assays$RNA@var.features
length(CKD_val@assays$RNA@var.features)

###scaling the data###
all.genes <- rownames(CKD_val)
CKD_val <- ScaleData(CKD_val,vars.to.regress = c("percent.mt"))###处理线粒体基因
CKD_val

###perform linear dimensional reduction###
CKD_val <- RunPCA(CKD_val, features = VariableFeatures(object = CKD_val))
print(CKD_val[["pca"]], dims = 1:5, nfeatures = 5)

###PCA
VizDimLoadings(CKD_val, dims = 1:2, reduction = "pca")
ElbowPlot(CKD_val,ndims = 50)

####cluster the cells###
CKD_val <- FindNeighbors(CKD_val, dims = 1:30)

####Run non-linear dimensional reduction (UMAP/tSNE
CKD_val <- RunUMAP(CKD_val, dims = 1:30)
CKD_val<-RunTSNE(CKD_val,dims=1:30)

###resolution seting
CKD_val <- FindClusters(CKD_val, resolution = 0.8)

DimPlot(CKD_val, reduction = "umap",label = F)
DimPlot(CKD_val, reduction = "tsne",label = F)

###Find marker for subclusters
CKD.markers <- FindAllMarkers(CKD_val, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(CKD.markers,"./Subgroup.markers.csv")

###Rename
CKD_val<-RenameIdents(CKD_val,"5"="Epithelial","11"="Epithelial","17"="Epithelial","18"="Epithelial","19"="Epithelial","16"="Epithelial","29"="Epithelial","23"="Epithelial",
                         "0"="Endothelial","1"="Endothelial","2"="Endothelial","3"="Endothelial","7"="Endothelial","9"="Endothelial","10"="Endothelial","22"="Endothelial","30"="Endothelial","6"="Endothelial","12"="Endothelial","28"="Endothelial","34"="Endothelial",
                         "8"="Myeloid","25"="Myeloid","27"="Myeloid","21"="Myeloid" ,"15"="Macrophages","13"="Pericytes","4"="B cells","31"="Plasma","26"="T cells","24"="NK", "32"="Mast")

###
table(Idents(CKD_val))







