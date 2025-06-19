library(data.table)
library(tidyr)
library(Seurat)
library(dplyr)
library(patchwork)
library(Matrix)
library(RColorBrewer)
set.seed(101)

rm(list=ls())
gc()

color.vector <- c(brewer.pal(8, "Set1"), brewer.pal(12, "Set3"), brewer.pal(8, "Dark2"), brewer.pal(8, "Accent"))


setwd("/home/luchun/scRNA/Neutrophil/01_Human_BM/")


sce = dior::read_h5('05_anno_data.h5')
sce




########################################################
genelist <- read.csv("/home/luchun/scRNA/Neutrophil/Score_gene_list.csv",header=TRUE)
#genelist <- read.csv("/home/luchun/scRNA/Neutrophil/05_17_cancers/scanpy/marker_list2.csv",header=TRUE)
ls(genelist)

name <- colnames(genelist)

lists <- list()

for(i in name){
	lists[[i]] <- genelist[,i][!genelist[,i] == ""]

}


#HB1 <- genelist$HB1[! genelist$HB1 == ""]
#HB1_2 <- HB1[-grep("RPL|RPS|MT",HB1)]


		
	
###################  2. AddModuleScore   #################


sce <- AddModuleScore(sce,features = list(lists[["PM"]]),assay = "RNA",name = "PM_AddModuleScore",ctrl = 200)

sce <- AddModuleScore(sce,features = list(lists[["MC_MM"]]),assay = "RNA",name = "MC_MM_AddModuleScore",ctrl = 50)

sce <- AddModuleScore(sce,features = list(lists[["BD_SC"]]),assay = "RNA",name = "BD_SC_AddModuleScore",ctrl = 200)



'''
for(i in name){
	lis <- list(lists[[i]])
	j <- paste(i,"AddModuleScore",sep="_")
	
	sce <- AddModuleScore(sce,features = lis,assay = "RNA",name = j)
}
'''

head(sce@meta.data)




##########保存
sce.all <- sce
sce@assays$RNA@data <- as.matrix(0)
sce@assays$RNA@scale.data <- as.matrix(0)


saveRDS(sce,"rds/05_sce_score.rds")




#——————————————————————可视化——————————————————————————

sce <- sce.all

sce$leiden_res_1 <- factor(sce$leiden_res_1,levels=c(0,1,2,6,4,3,5,7,8))

mycolor <- rev(brewer.pal(11,"RdBu")[2:10])

test1 <- paste(name,"AddModuleScore1",sep="_")

index <- which(sce$leiden_res_1 == 2)
sce$Immature_AddModuleScore1[index] <- sce$Immature_AddModuleScore1[index] +0.02

index <- which(sce$leiden_res_1 == 7)
sce$Immature_AddModuleScore1[index] <- sce$Immature_AddModuleScore1[index] +0.02
sce$Mature_AddModuleScore1[index] <- sce$Mature_AddModuleScore1[index] -0.02


dp <- DotPlot(sce, features = test1,dot.scale = 10,scale.min = 5,assay='RNA',group.by = "leiden_res_1",scale = TRUE)+
  coord_flip()+
  theme_bw()+RotatedAxis()+
  theme(panel.grid = element_blank(),
	aspect.ratio = 1/3,
    axis.title.x = element_text(size = 20),  # 设置x轴标签字体大小
    axis.title.y = element_text(size = 20),  # 设置y轴标签字体大小
    axis.text.x = element_text(size = 15),    # 设置x轴刻度字体大小
    axis.text.y = element_text(size = 15),     # 设置y轴刻度字体大小
	legend.title = element_text(size = 15),    # 设置图例标题字体大小
    legend.text = element_text(size = 12))+
  labs(x=NULL,y=NULL,title = "AddModuleScore")+
  guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = mycolor)
  

jpeg("03-1_DotPlot_AddModuleScore.jpeg",width=700,height=500)
dp
dev.off()

tiff("03-1_DotPlot_AddModuleScore.tiff",width=700,height=500)
dp
dev.off()



dp <- DotPlot(sce, features = test1,dot.scale = 10,scale.min = 5,assay='RNA',group.by = "Neu_anno",scale = TRUE)+
  coord_flip()+
  theme_bw()+RotatedAxis()+
  theme(panel.grid = element_blank(),
	aspect.ratio = 1/1,
    axis.title.x = element_text(size = 20),  # 设置x轴标签字体大小
    axis.title.y = element_text(size = 20),  # 设置y轴标签字体大小
    axis.text.x = element_text(size = 15),    # 设置x轴刻度字体大小
    axis.text.y = element_text(size = 15),     # 设置y轴刻度字体大小
	legend.title = element_text(size = 15),    # 设置图例标题字体大小
    legend.text = element_text(size = 12))+
  labs(x=NULL,y=NULL,title = "AddModuleScore")+
  guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = mycolor)
  

jpeg("03-2_DotPlot_AddModuleScore.jpeg",width=500,height=400)
dp
dev.off()

tiff("03-2_DotPlot_AddModuleScore.tiff",width=500,height=400)
dp
dev.off()