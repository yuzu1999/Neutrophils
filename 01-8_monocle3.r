library(monocle3)
library(Seurat)
library(tidyverse)
library(RColorBrewer)
#加载数据
#devtools::install_github('satijalab/seurat-data')

rm(list=ls())
gc()

setwd("/home/luchun/scRNA/Neutrophil/01_Human_BM/monocle2/")

sce <- dior::read_h5('../05_anno_data.h5')
sce
#17610 features across 21861 samples within 1 assay



#————————————————构建Monocle3所需数据——————————————————
expression_matrix = sce@assays$RNA@data
cell_metadata = data.frame(sce@meta.data)

gene_annotation = data.frame(expression_matrix[,1])
gene_annotation[,1] = row.names(gene_annotation)
colnames(gene_annotation)=c("gene_short_name")


##构建Monocle3 cds对象
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)


#预处理，相当于Seurat流程normalize过程
#counts数据选择norm_method = c("log")
#data数据选择norm_method = c("none")
cds <- preprocess_cds(cds, num_dim = 50,norm_method = c("none"))


#去除批次效应,多个样本时可以通过此方式去除批次
cds <- align_cds(cds, alignment_group = "orig.ident")


## 降维，默认是"Umap"方式
cds <- reduce_dimension(cds, preprocess_method = "PCA")

color <- c("#F46D43","#E6F598","#66C2A5","#5E4FA2")

jpeg("01_Neu_anno.jpg")
plot_cells(cds, reduction_method="UMAP", 
           show_trajectory_graph = FALSE,
           label_cell_groups = FALSE, 
           color_cells_by="Neu_anno")+ scale_color_manual(values = color)
dev.off()		   



##使用Seurat的UMAP信息，这样可以与Seurat对象的细胞分布保持一致
seurat_umap <- Embeddings(sce, reduction = "umap")[colnames(cds),]
cds@int_colData$reducedDims$UMAP <- seurat_umap

jpeg("02_Neu_anno.jpg")
plot_cells(cds, reduction_method="UMAP", 
           show_trajectory_graph = FALSE,
           label_cell_groups = FALSE, 
           color_cells_by="Neu_anno")+ scale_color_manual(values = color)
dev.off()




#——————————————————轨迹分析————————————————————————————
#分群(类似monocle的State)
cds <- cluster_cells(cds)

#预测轨迹
cds <- learn_graph(cds)

#交互式确定root节点，可以选择多个。我这里选择了一个
cds <- order_cells(cds)



##不同细胞类型拟时序数值
##拟时序值越高表示细胞分化程度越高，这里仅为演示，并非真实分化情况
jpeg("03_pseu.jpg")
plot_cells(cds, color_cells_by = "pseudotime",show_trajectory_graph=TRUE) 
dev.off()



saveRDS(cds,"cds_monocle3.rds")

cds <- readRDS('cds_monocle3.rds')




library(RColorBrewer)

# 将伪时间信息添加到 Seurat 对象的元数据中
pseudotime <- cds@principal_graph_aux@listData$UMAP$pseudotime


sce <- AddMetaData(
  object = sce,
  metadata = pseudotime,
  col.name = "Monocle3_Pseudotime"
)

mycolor2 <- rev(brewer.pal(10, "Blues"))[1:8]

tiff("04_pseu.tiff",width=500)
FeaturePlot(sce,features="Monocle3_Pseudotime",reduction='umap',pt.size=0.5)+ scale_color_gradientn(
	colors = mycolor2,
	limits = c(0, 50),  
    na.value = mycolor2[11])
dev.off()





# ——————————————————绘制拟合曲线图——————————————————————
# 提取伪时间信息
library(splines)
library(ggplot2)
library(cowplot)


pseudotime <- cds@principal_graph_aux@listData$UMAP$pseudotime

genes <- unique(c('CEBPB', 'S100A8', 'S100A9', 'IRF1','STAT1','STAT3','STAT6','ANXA1','LYZ2','CD14','TGFB1','RB1', 'ARG1','ARG2','NOS2','CXCR1', 'CYBB', 'VEGFA', 'CD274', 'TGFB1', 'IL4R','CD84', 'TNFRSF10B', 'OLR1','MKI67','TOP2A','PCNA','CD63','LGALS3','CEBPE'))

genes2 <- intersect(genes,rownames(assays(cds)[["counts"]]))
genes2

p_list <- list()


for(i in genes2){
	gene_expression <- assays(cds)[["counts"]][i, ]

	# 创建拟合数据框
	fit_data <- data.frame(pseudotime = pseudotime, gene_expression = gene_expression)

	# 使用 splines 包的 ns 函数拟合平滑样条
	fit <- lm(gene_expression ~ ns(pseudotime, df = 3), data = fit_data)

	# 生成拟合曲线的数据点
	new_pseudotime <- seq(min(pseudotime), max(pseudotime), length.out = 100)
	fit_values <- predict(fit, newdata = data.frame(pseudotime = new_pseudotime), interval = "confidence")

	# 创建拟合曲线的数据框
	fit_curve <- data.frame(
		pseudotime = new_pseudotime,
		fit = fit_values[, "fit"],
		lower = fit_values[, "lwr"],
		upper = fit_values[, "upr"]
	)



	# 绘制拟合曲线和误差区间
	p <- ggplot() +
		geom_ribbon(data = fit_curve, aes(x = pseudotime, ymin = lower, ymax = upper), fill = "#9ecae1", alpha = 0.4) +
		geom_line(data = fit_curve, aes(x = pseudotime, y = fit), color = "#223D6C", size = 2) +
		labs(title = i, x = "Pseudotime", y = "Expression") +
		theme_bw()+
		theme(
			aspect.ratio = 1/1,
			#panel.grid = element_blank(),
			axis.title.x = element_text(size = 16),  
			axis.title.y = element_text(size = 16),  
	
			axis.text.x = element_text(size = 14),  
			axis.text.y = element_text(size = 14),
			plot.title = element_text(size = 25, face = "bold")
		)
		
	p_list[[i]] <- p
}

tiff("05_expr.tiff",width=2800,height=1600)
plot_grid(plotlist=p_list,ncol=7)
dev.off()


jpeg("05_expr.jpg",width=2800,height=1600)
plot_grid(plotlist=p_list,ncol=7)
dev.off()




# ——————————————————绘制峰状图————————————————————————————————
library(ggplot2)
library(ggridges)

# 假设sce是你的Seurat对象
sce$leiden_res_1 <- factor(sce$leiden_res_1,levels=c(0,1,2,6,4,3,5,7,8))
meta_data <- sce@meta.data
sorted_meta_data <- meta_data[order(meta_data$Monocle3_Pseudotime), ]

color <- c("#F46D43","#E6F598","#66C2A5","#5E4FA2")
# 绘制平滑峰状图
p <- ggplot(sorted_meta_data, aes(x = Monocle3_Pseudotime, y = Neu_anno, fill = Neu_anno)) +
  geom_density_ridges() +
  scale_fill_manual(values = color)+ 
    labs(title = "Neu_anno",
		x = "Pseudotime",
		y = "",
		fill = "Neu_anno") +
	theme_bw()+
	theme(
		#panel.grid = element_blank(),
		aspect.ratio = 1/2,
		axis.title.x = element_text(size = 16),  
		axis.title.y = element_text(size = 16),  
	
		axis.text.x = element_text(size = 14),  
		axis.text.y = element_text(size = 14),
		plot.title = element_text(size = 20, face = "bold")
	)

tiff("06_ridge_neu_anno.tiff",width=600,height=300)
p
dev.off()


color <- c("#ffe788","#e8743c","#6a51a3","#CE9FCA","#9C9BE9","#76afda","#9e9ac8","#2873B3","#b0d45d")[c(1,2,3,7,5,4,6,8,9)]

p <- ggplot(sorted_meta_data, aes(x = Monocle3_Pseudotime, y = leiden_res_1, fill = leiden_res_1)) +
  geom_density_ridges() +
  scale_fill_manual(values = color)+ 
    labs(title = "leiden_res_1",
		x = "Pseudotime",
		y = "",
		fill = "leiden_res_1") +
	theme_bw()+
	theme(
		#panel.grid = element_blank(),
		aspect.ratio = 1/2,
		axis.title.x = element_text(size = 16),  
		axis.title.y = element_text(size = 16),  
	
		axis.text.x = element_text(size = 14),  
		axis.text.y = element_text(size = 14),
		plot.title = element_text(size = 20, face = "bold")
	)

tiff("06_ridge_leiden_res_1.tiff",width=600,height=300)
p
dev.off()




#——————————————基因聚类——————————————————————————

cds <- readRDS("cds_monocle3.rds")

#识别沿轨迹差异表达的基因
modulated_genes <- graph_test(cds, neighbor_graph = "principal_graph", cores = 4)
dim(modulated_genes)
#[1] 17610     6


grep("^HSP|^MT-|^RPS|^RPL",rownames(modulated_genes),value=TRUE)

index <- grep("^HSP|^MT-|^RPS|^RPL",rownames(modulated_genes))
modulated_genes2 <- modulated_genes[-index,]

Track_genes <- modulated_genes2[,c(5,2,3,4,1,6)] %>% filter(q_value < 1e-3)
write.csv(Track_genes, "Trajectory_genes.csv", row.names = F)



####寻找共表达基因模块
pr_deg_ids <- row.names(subset(modulated_genes2, q_value == 0))
length(pr_deg_ids)
#[1] 959

'''
deg <- read.csv("/home/luchun/scRNA/Neutrophil/01_Human_BM/leiden_res_1_deg_filter.csv",header=TRUE,row.names=1)
pr_deg_ids <- unique(deg$names)
length(pr_deg_ids)
#[1] 3993
'''


library(ClusterGVis)
mat <- pre_pseudotime_matrix(cds_obj = cds,
                             gene_list = pr_deg_ids)

# check
head(mat[1:5,1:5])
dim(mat)
#[1]   959 21861


# kmeans
ck <- clusterData(exp = mat,
                  cluster.method = "kmeans",
                  cluster.num = 5)

# add line annotation
pdf('07_monocle3.pdf',height = 10,width = 8,onefile = F)
visCluster(object = ck,
           plot.type = "both",
           add.sampleanno = F#,
           #markGenes = sample(rownames(mat),30,replace = F)
		   )
dev.off()