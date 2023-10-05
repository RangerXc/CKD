

library(Seurat)
library(ggplot2)
library(forcats)
library(ggstance)
library(clusterProfiler)
library(msigdbr)
library(enrichplot)
library(ggstatsplot)

Fibroblasts <- readRDS("Fibroblasts.rds")

###
Fibroblasts@meta.data$CD248 <- Fibroblasts@assays$RNA@counts["CD248",]
Fibroblasts@meta.data$type <- ifelse(Fibroblasts@meta.data$CD248 > 0, "CD248+", "CD248-")
Fibroblasts@meta.data[1:5,]

###
Fibroblasts <- ScaleData(Fibroblasts, features = rownames(Fibroblasts))
###
Idents(Fibroblasts) <- as.factor(Fibroblasts@meta.data$type)
Fibroblasts@meta.data[1:5,]
###DEGs for CD248+ vs CD248-
diff <- FindMarkers(Fibroblasts, ident.1 = "CD248+", ident.2 = "CD248-", logfc.threshold = 0)
diff <- diff[order(diff$avg_log2FC, decreasing = T),]
diff_list <- diff$avg_log2FC
names(diff_list) <- rownames(diff)

###Vocano plot
library(EnhancedVolcano)
#设置工作目录并载入本地数据：
data <- read.csv("aaa.csv",header=T)###提取csv格式文件

#突出所关注标签变量：
vals<-c('CD74','S100A6','HLA-DRB1','POSTN','FN1','COL3A1','COMP','COL8A1','COL1A1','VCAN','COL6A3')

#自定义颜色，现在将log2FC>1.5（up）指定为橙色，将log2FC<-1.5（down）指定为蓝色，其余（nodiff）为灰色:
#指定标签向量group，通过ifelse条件判断函数指定标签和对应颜色
group<-ifelse(
  data$log2FC<(-0.5)&data$Pvalue<0.05,'#4D4398',
  ifelse(data$log2FC>(0.5)&data$Pvalue<0.05,'#F18D00',
         '#b5b5b5'))
group[is.na(group)]<-'#b5b5b5'
names(group)[group=='#F18D00']<-'Up'
names(group)[group=='#b5b5b5']<-'Nodiff'
names(group)[group=='#4D4398']<-'Down'

downvals<-c('CD74','S100A6','HLA-DRB1')
upvals<-c('POSTN','FN1','COL3A1','COMP','COL8A1','COL1A1','VCAN','COL6A3')

#默认原始参数绘制基本火山图：
EnhancedVolcano(data,
                x="log2FC",
                y="Pvalue",
                lab=data$id,
                pCutoff=10e-1/20,#y轴阈值线(水平)
                FCcutoff=0.5,#x轴阈值线（垂直）
                labSize=3,#标签大小
                xlim=c(-2, 2),#限制X轴范围
                ylim=c(0,30),#限制Y轴范围
                selectLab=c(vals),#使用selectLab参数选定所关注的标签
                xlab=bquote(~Log[2]~'fold change'),#将内容传递给xlab
                labCol='black',#标签颜色
                labFace='bold',#标签字体
                boxedLabels=TRUE,#是否在框中绘制标签
                drawConnectors=TRUE,#是否通过连线将标签连接到对应的点上
                widthConnectors=0.8,#连线的宽度
                endsConnectors= "first",#连线绘制箭头的方向，可选first、both、last
                colConnectors='grey',#连线的颜色
                colCustom=group,#用group覆盖默认配色方案
                colAlpha=0.6,#调整透明度
                cutoffLineType='twodash',#阈值线类型，可选“blank”、“solid”、“dashed”、“dotted”、“dotdash”、“longdash”和“twodash”
                cutoffLineCol='pink',#阈值线颜色
                cutoffLineWidth=0.88,#阈值线粗细
                title="Volcano Plot Exp",#主标题
                subtitle="Differential expression",#副标题
                caption=bquote(~Log[2]~"fold change cutoff,1;p-value cutoff,0.05"),#注释说明
                legendPosition='right',#图例位置
                legendLabSize=12,#图例文字大小
                legendIconSize=6,#图例符号大小
                pointSize=c(ifelse(data$id %in% vals,5,3)))+#通过ifelse条件判断函数指定散点大小
  coord_flip()###翻转



