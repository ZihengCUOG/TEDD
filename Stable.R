library(Seurat)

setwd("/data/Stable")
load("Tedd.5.Rda")

meta <- merge.data@meta.data[,c("Tissue","Celltype","Timepoint","Sex")]

meta$ID <- rownames(meta)

counts <- merge.data@assays$RNA@counts

count <- counts[1,,drop = F]
merge <- sapply(split(meta$ID, meta$Tissue),function(cells) sum(count[,cells,drop = F] > 0)/ncol(count[,cells,drop = F]))

for (i in 2:nrow(counts)){
  count <- counts[i,,drop = F]
  tissue <- sapply(split(meta$ID, meta$Tissue),
                   function(cells) sum(count[,cells,drop = F] > 0)/ncol(count[,cells,drop = F]))
  merge <- cbind(merge,tissue)
}

colnames(merge) <- rownames(counts)

merge <- t(merge)

save(merge, file = "tissue.Rda")

count <- counts[1,,drop = F]
merge2 <- sapply(split(meta$ID, meta$Timepoint),
                 function(cells) sum(count[,cells,drop = F] > 0)/ncol(count[,cells,drop = F]))

for (i in 2:nrow(counts)){
  count <- counts[i,,drop = F]
  tissue <- sapply(split(meta$ID, meta$Timepoint),
                   function(cells) sum(count[,cells,drop = F] > 0)/ncol(count[,cells,drop = F]))
  merge2 <- cbind(merge2,tissue)
}

colnames(merge2) <- rownames(counts)

merge2 <- t(merge2)

save(merge2, file = "timepoint.Rda")

load("tissue.Rda")
load("timepoint.Rda")

t <- merge[,1,drop = F][merge[,1] > 0.6,]

t <- data.frame(t)

colnames(t) <- "Ratio"

t$Gene <- rownames(t)

t$Ident <- colnames(merge)[1]


for (i in 2:ncol(merge))try({
  t2 <- merge[,i,drop = F][merge[,i] > 0.6,,drop = F]
  t2 <- data.frame(t2)
  colnames(t2) <- "Ratio"
  t2$Gene <- rownames(t2)
  t2$Ident <- colnames(merge)[i]
  t <- rbind(t,t2)
},silent = T)

tissue <- t

t <- merge2[,1,drop = F][merge2[,1] > 0.6,]

t <- data.frame(t)

colnames(t) <- "Ratio"

t$Gene <- rownames(t)

t$Ident <- colnames(merge2)[1]


for (i in 2:ncol(merge2))try({
  t2 <- merge2[,i,drop = F][merge2[,i] > 0.6,,drop = F]
  t2 <- data.frame(t2)
  colnames(t2) <- "Ratio"
  t2$Gene <- rownames(t2)
  t2$Ident <- colnames(merge2)[i]
  t <- rbind(t,t2)
},silent = T)

timepoint <- t

t <- rbind(tissue, timepoint)

write.csv(t, file = "h.stable.csv",row.names = F)
