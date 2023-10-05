
library(Seurat)
library(dplyr)
library(monocle)
library(ggsci)
library(RColorBrewer)

##get the monoand fibroblasts
Fibroblasts$cell_type_val<-Idents(Fibroblasts)
Fibroblasts<-FindVariableFeatures(Fibroblasts)

matrix<-as.matrix(Fibroblasts@assays$RNA@counts)
dim(matrix)
matrix[1:500,1:10]
gene_ann <- data.frame(
  gene_short_name = row.names(matrix), 
  row.names = row.names(matrix)
)
#oct[["cell_group"]]<-Idents(oct)
sample_ann <- Fibroblasts@meta.data

fd <- new("AnnotatedDataFrame",data=gene_ann)
pd<-new("AnnotatedDataFrame",data=sample_ann)
#?newCellDataSet
sc_cds_2 <- newCellDataSet(matrix,  phenoData = pd,featureData =fd,expressionFamily = negbinomial.size(),lowerDetectionLimit=0.01)
sc_cds_2 <- estimateSizeFactors(sc_cds_2)
sc_cds_2 <- estimateDispersions(sc_cds_2)
####QC#####
sc_cds_2 <- detectGenes(sc_cds_2, min_expr = 1)
expressed_genes <- row.names(subset(fData(sc_cds_2), num_cells_expressed >= 5))

fData(sc_cds_2)[1:5,]

#####cluster analysis######
diff_test_res <- differentialGeneTest(sc_cds_2[expressed_genes,],
                                      fullModelFormulaStr = "~cell_type_val",cores = 7,
                                      verbose = T) #+num_genes_expressed+orig.ident
ordering_genes <- row.names(subset(diff_test_res, qval < 1e-2)) #1e-200

#ordering_genes<-markers %>% group_by(cluster) %>% top_n(n=200,wt=avg_logFC)
sc_cds2 <- setOrderingFilter(sc_cds_2, ordering_genes)
plot_ordering_genes(sc_cds2)
sc_cds2<- reduceDimension(sc_cds2, max_components = 2, num_dim=6,reduction_method  = "DDRTree")#,residualModelFormulaStr = "~orig.ident")
sc_cds2 <- orderCells(sc_cds2)
plot_cell_trajectory(sc_cds2, color_by = "cell_type_val",show_branch_points = T)

#p1<-plot_cell_trajectory(sc_cds2, color_by = "RNA_snn_res.0.8",show_branch_points = F)#+facet_wrap(~RNA_snn_res.0.8)
plot_cell_trajectory(sc_cds2, color_by = "State",show_branch_points = T)#+facet_wrap(~Site)
sc_cds2=orderCells(sc_cds2,root_state = 3) 
plot_cell_trajectory(sc_cds2,color_by = "Pseudotime",show_branch_points = F)
plot_cell_trajectory(sc_cds2, markers = c("CCL21"),use_color_gradient = T,show_branch_points = F)
to_be_tested <- row.names(subset(fData(sc_cds2), gene_short_name %in% c("CCL2")))
cds_subset <- sc_cds2[to_be_tested,]
plot_genes_in_pseudotime(cds_subset, color_by = "cell_type_val")

####get the genes change along with the Pseudotime###
diff_test_res <- differentialGeneTest(sc_cds2[expressed_genes,],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)",cores = 7)
sig_gene_names <- row.names(subset(diff_test_res, qval < 0.05))
#
pseudoplot<-plot_pseudotime_heatmap(sc_cds2[sig_gene_names,],
                                    num_clusters = 2,
                                    cores = 6,
                                    show_rownames = T,return_heatmap = T,
                                    hmcols = colorRampPalette(rev(brewer.pal(9, "PRGn")))(62))

####get the gene names in each cluster##
clusters <- cutree(pseudoplot$tree_row, k = 2)
clustering <- data.frame(clusters)

clustering[,1] <- as.character(clustering[,1])
colnames(clustering) <- "Gene_Clusters"
table(clustering)

write.csv(em_group_pathway_2,"./em_group_pathway_2.csv")

###Enrichment analysis
library(msigdbr)
library(clusterProfiler)
msigdbr_show_species()  
m_df <- msigdbr(species = "Homo sapiens") 

###Go/KEGG
m_t2g_C5 <- msigdbr(species = "Homo sapiens", category = "C5",subcategory = "BP") %>% 
  dplyr::select(gs_name, gene_symbol)

m_t2g_C7 <- msigdbr(species = "Homo sapiens", category = "C7") %>% 
  dplyr::select(gs_name, gene_symbol)
m_t2g_C2 <- msigdbr(species = "Homo sapiens", category = "C2",subcategory = c("CP:KEGG")) %>% 
  dplyr::select(gs_name, gene_symbol) ###signating pathways

###Cluster 1
gene1<-rownames(clustering)[clustering[,1]==1]
em_group1 <- enricher(gene1, TERM2GENE=m_t2g_C5)
em_group_pathway_1 <- enricher(gene1, TERM2GENE=m_t2g_C2)
###simplify(em_group1)
head(em_group1)
head(em_group_pathway_1)
barplot(em_group1,showCategory = 20)
barplot(em_group_pathway_1)
write.csv(em_group1,file="./GO_1.csv")
write.csv(em_group_pathway_1,file="./KEGG_1.csv")

###Cluster 2
gene2<-rownames(clustering)[clustering[,1]==2]
em_group2<- enricher(gene2, TERM2GENE=m_t2g_C5)
em_group_pathway_2 <- enricher(gene2, TERM2GENE=m_t2g_C2)
###simplify(em_group2)
head(em_group2)
head(em_group_pathway_2)
barplot(em_group2,showCategory = 20)
barplot(em_group_pathway_2)
write.csv(em_group2,file="./GO_2.csv")
write.csv(em_group_pathway_2,file="./KEGG_2.csv")

###Bar plot
KEGG_result <- read.csv("GGO_2.csv",header=T)###提取csv格式文件
##计算Rich Factor
KEGG_result <- mutate(KEGG_result,
                      RichFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))
##计算Fold Enrichment
KEGG_result$FoldEnrichment <- apply(KEGG_result,1,function(x){
  GeneRatio <- eval(parse(text = x["GeneRatio"]))
  BgRatio <- eval(parse(text = x["BgRatio"]))
  foldEnrichment <- round(GeneRatio/BgRatio,2)
  foldEnrichment
})
head(KEGG_result$RichFactor)
head(KEGG_result$FoldEnrichment)
colnames(KEGG_result)

#提取结果前Top15绘图(或自定义所需pathway绘图)：
top15 <- KEGG_result[1:15,]
#指定绘图顺序（转换为因子）
top15$pathway <- factor(top15$Description,levels = rev(top15$Description))

#自定义主题：
mytheme <- theme(axis.title = element_text(size = 13),
                 axis.text = element_text(size = 11,face = "bold"),
                 plot.title = element_text(size = 18,###表头字体大小
                                           hjust = 0.5,
                                           face = "bold"),
                 legend.title = element_text(size = 12),
                 legend.text = element_text(size = 12))

###可选绘图方案
#Top15富集数目条形图（横轴基因数，填充-log10（Pvalue））：
p <- ggplot(data = top15,
            aes(x = Count, y = pathway, fill = -log10(pvalue)))+
  scale_fill_distiller(palette = "RdPu",direction = 1) +
  geom_bar(stat = "identity", width = 0.8) +
  labs(x = "Number of Gene", y = "pathway", title = "Cluster1 GO enrichment") +
  theme_bw() + mytheme
p

######
#Top15富集因子条形图（挑选所需数据列绘图即可）（横轴-log10（Pvalue）,填充RichFactor)：
p1 <- ggplot(data = top15,
             aes(x = -log10(pvalue), y = pathway, fill = RichFactor)) +
  scale_fill_distiller(palette = "Blues",direction = 1) +
  geom_bar(stat = "identity", width = 0.8) +
  labs(x = "-log10(pvalue)", y = "pathway", title = "Cluster2 GO enrichment") +
  theme_bw() + mytheme
p1

