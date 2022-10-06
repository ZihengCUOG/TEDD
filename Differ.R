library(Seurat)

merge.data <- readRDS("save.Rds")

library(future)
plan()
options(future.globals.maxSize = 8000 * 1024^2)
plan("multiprocess",workers = 16)
start = Sys.time()

merge.data <- PrepSCTFindMarkers(merge.data)

merge.data@meta.data$Tissue_Timepoint <- paste(merge.data@meta.data$Tissue, "_", "All","_",merge.data@meta.data$Timepoint, sep = "")
merge.data@active.ident <- as.factor(merge.data$Tissue_Timepoint)
c <- FindAllMarkers(merge.data,logfc.threshold = 0)

save(c, file = "c.Rda")


merge.data@meta.data$Tissue_Celltype_Timepoint <- paste(merge.data@meta.data$Tissue, "_",merge.data@meta.data$Celltype,"_",merge.data@meta.data$Timepoint, sep = "")

sum <- data.frame(table(merge.data@meta.data$Tissue_Celltype_Timepoint))

sum <- sum[sum$Freq > 10,]

cell <- merge.data@meta.data[merge.data@meta.data$Tissue_Celltype_Timepoint %in% sum$Var1,]

merge.data <- subset(merge.data, cells = rownames(cell))

merge.data@active.ident <- as.factor(merge.data$Tissue_Celltype_Timepoint)

merge.data <- PrepSCTFindMarkers(merge.data)

d <- FindAllMarkers(merge.data,logfc.threshold = 0)

save(d, file = "d.Rda")

a = list.files(,pattern = "c.Rda")
a
dir = paste("./",a,sep="")
n = length(dir)

load("Tedd.10_AdrenalGland_scrna_c.Rda")

data <- c

for (i in 2:n){
  load(dir[i])
  new.data <- c
  data <- rbind(data,new.data)
  #   head(merge.data)
}

#data <- data.frame(data)

a = list.files(,pattern = "d.Rda")
a
dir = paste("./",a,sep="")
#命令构建路径变量dir（方便更改），也可以不构建，后面示例
n = length(dir)

load("Tedd.10_AdrenalGland_scrna_d.Rda")

data2 <- d

for (i in 2:n){
  load(dir[i])
  new.data <- d
  data2 = rbind(data2,new.data)
  #   head(merge.data)
}

c <- data
d <- data2

save(c, file = "c.Rda")
save(d, file = "d.Rda")



load("c.Rda")
load("d.Rda")

c$Timepoint <- sub(".*_", "",c$cluster)

c$Tissue <- sub("_.*", "",c$cluster)

#c <- c[!c$Timepoint %in% c("19 h after routine fertilization","27 h after routine fertilization","4 h after IVF oocyte retrieval",
#"48 h after routine fertilization","E3","E4"),]

#c$cluster <- paste(c$Tissue, "_", "All","_",c$Timepoint, sep = "")

c <- c[,1:7]


meta <- read.csv("h.meta.csv")
meta$tissue <- paste(meta$Tissue, "_", meta$Celltype,"_",meta$Timepoint, sep = "")
ttt <- data.frame(table(meta$tissue))
tttt <- ttt[ttt$Freq > 50,]

id <- intersect(d$cluster, tttt$Var1)

d <- d[d$cluster %in% id, ]

gene <- rbind(c,d)

write.csv(gene, "h.differ.gene.csv",row.names = F)

#plot figure 3d

library(pheatmap)
mat <- read.csv("/data/Differh.multigene.csv",row.names = 1, check.names = F)
gene <- read.csv("/data/Differ/h.differ.gene.csv")
gene2 <- gene[gene$cluster == "Liver_All_GW13" & abs(gene$avg_log2FC) > 1, ]
mat <- mat[gene2$gene, c("Liver_All_21-30","Liver_All_51-60", "Liver_All_GW13","Liver_All_GW16", "Liver_All_GW17","Liver_All_GW26")]
pheatmap(mat, scale = "row")

