a = list.files(,pattern = ".tsv")
a
dir = paste("./",a,sep="")

n = length(dir)

merge.data <- read.table(file = dir[1], header=TRUE, sep="\t")

data <- c(a[1],ncol(merge.data))

for (i in 2:n){
   new.data = read.table(file = dir[i], header=T, sep="\t")
   new.data <- c(a[i],ncol(new.data))
   data = rbind(data,new.data)
   head(merge.data)
}

data <- data.frame(data)

a <- data[as.character(data$X2) %in% c(15,17),]

a <- as.character(a[,1])

a
dir = paste("./",a,sep="")
n = length(dir)
merge.data <- read.table(file = dir[1], header=TRUE, sep="\t")
merge.data <- merge.data[,c("gene_id","TPM"),drop = F]
col.name <- sub('.*/','',dir[1])
col.name <- substr(col.name,1,nchar(col.name)-4)
colnames(merge.data) <- c("id",col.name)
merge.data <- merge.data[!duplicated(merge.data$id),]

for (i in 2:n){
   new.data = read.table(file = dir[i], header=T, sep="\t")
   new.data <- new.data[,c("gene_id","TPM"),drop = F]
   new.data <- new.data[!duplicated(new.data$gene_id),]
   col.name <- sub('.*/','',dir[i])
   col.name <- substr(col.name,1,nchar(col.name)-4)
   colnames(new.data) <- c("id",col.name)
   merge.data = merge(merge.data,new.data, all=TRUE,by = "id")
   head(merge.data)
}

write.csv(merge.data,file = "./merge_all.csv",row.names=FALSE)
saveRDS(merge.data, file = "./merge_all.Rds")

a <- read.csv("merge_all.csv")
b <- read.table(file = "ENCFF002MXO.tsv", header=TRUE, sep="\t", row.names=1)

a$id <- sub('\\.[0-9]*$', '',a$id)

a <- a[!duplicated(a$id),]

rownames(a) <- a$id

a <- a[,-1]

httr::set_config(httr::config(ssl_verifypeer = FALSE))
library(biomaRt)
ensembl <- useMart('ensembl', dataset = 'hsapiens_gene_ensembl')
genes_ids <- sub('\\.[0-9]*$', '',rownames(a))
gs_heatdata <- getBM(
  attributes = c('external_gene_name', 'hgnc_symbol','ensembl_gene_id'),
  filters = 'ensembl_gene_id',
  values = genes_ids,
  mart = ensembl)

gs_heatdata

gs_heatdata <- gs_heatdata[!duplicated(gs_heatdata$ensembl_gene_id),]

rownames(gs_heatdata) <- gs_heatdata$ensembl_gene_id

genes_ids <- data.frame(genes_ids)

rownames(genes_ids) <- genes_ids$genes_ids

genes_ids$symbol <- gs_heatdata[rownames(genes_ids),"hgnc_symbol"]

a$symbol <- genes_ids$symbol

a <- a[!is.na(a$symbol),]

saveRDS(a, file = "merge_all.Rds")

a <- a[!duplicated(a$symbol),]
rownames(a) <- a$symbol

a <- a[,-ncol(a)]

write.csv(a,file = "./merge_all.csv", quote=FALSE)

a <- read.csv("merge_all.csv",row.names = 1)

##linux
awk 'BEGIN{FS=OFS="\t"} $3=="gene"{print $0}' gencode.v29.annotation.gff3 |sed 's/;/\t/g' |awk 'BEGIN{FS=OFS="\t"} {for(i=1; i<=NF; i++){split($i, a, "="); b[a[1]]=a[2]}} {print b["gene_name"],b["gene_type"]}{split("", b, ":")}' > test.txt

b$id <- sub('\\.[0-9]*$', '',rownames(b))

b <- b[!duplicated(b$id),]

rownames(b) <- b$id

httr::set_config(httr::config(ssl_verifypeer = FALSE))
library(biomaRt)
ensembl <- useMart('ensembl', dataset = 'hsapiens_gene_ensembl')
genes_ids <- sub('\\.[0-9]*$', '',rownames(b))
gs_heatdata <- getBM(
  attributes = c('external_gene_name', 'hgnc_symbol','ensembl_gene_id'),
  filters = 'ensembl_gene_id',
  values = genes_ids,
  mart = ensembl)

gs_heatdata

gs_heatdata <- gs_heatdata[!duplicated(gs_heatdata$ensembl_gene_id),]

rownames(gs_heatdata) <- gs_heatdata$ensembl_gene_id

genes_ids <- data.frame(genes_ids)

rownames(genes_ids) <- genes_ids$genes_ids

genes_ids$symbol <- gs_heatdata[rownames(genes_ids),"hgnc_symbol"]

b$symbol <- genes_ids$symbol

b <- b[!is.na(b$symbol),]

b <- b[,c("symbol","TPM")]

b <- b[!duplicated(b$symbol),]

rownames(b) <- b$symbol

b <- b[,"TPM",drop = F]

colnames(b) <- "ENCFF002MXO"

a$ENCFF002MXO <- b[rownames(a), "ENCFF002MXO"]

saveRDS(a, file = "merge_all.Rds")

write.csv(a,file = "./merge_all.csv", quote=FALSE)

load("meta.h.Rda")

a <- a[,intersect(rownames(meta.h),colnames(a))]

anno <- read.table("../../tsv-backup/test.txt",sep = "\t")

h <- rownames(a)

h <- data.frame(h)

rownames(h) <- h$h

anno <- anno[!duplicated(anno$V1),]

rownames(anno) <- anno$V1

h$gene_type <- anno[rownames(h),"V2"]

write.table(h, file = 'human.anno.txt',sep = "\t", quote=FALSE,row.names = F)

genetype.anno <- data.frame(table(h$gene_type))

genetype.anno <- genetype.anno[order(genetype.anno$Freq,decreasing = T),]

write.table(genetype.anno, file = 'human.anno.sum.txt',sep = "\t", quote=FALSE,row.names = F)

a <- readRDS("merge_all.Rds")
a[is.na(a)] <- 0
write.csv(a,file = "./merge_all.csv", quote=FALSE)

a <- read.csv("merge_all.csv",row.names = 1)
load("meta.h.Rda")
a <- a[,rownames(meta.h)]
write.csv(a,file = "./h.merge_all.csv", quote=FALSE)

meta.h <- meta.h[,-6]
meta.h <- meta.h[,-4]
meta.h$Second.term.name <- paste(meta.h$Biosample.term.name, meta.h$Timepoint, sep = "_")
meta.h$Second.term.name[grep("Unknown",meta.h$Second.term.name)] <- "NA"
sum <- data.frame(table(meta.h$Second.term.name))
sum <- sum[!sum$Var1 == "NA",]
data2 <- meta.h[,"Second.term.name",drop = F]
data2 <- data2[data2$"Second.term.name" %in% sum$Var1,,drop = F]
a <- a[,rownames(data2)]
mat <- t(a)
mat <- aggregate(mat, by=list(data2[,1]), FUN="mean")
rownames(mat) <- mat$Group.1
mat <- t(mat[,-1])
head(mat)
mat1 <- mat
a <- read.csv("h.merge_all.csv",row.names = 1)
data2 <- meta.h[meta.h$Timepoint == "Unknown","Biosample.term.name",drop = F]
a <- a[,rownames(data2)]
mat <- t(a)
mat <- aggregate(mat, by=list(data2[,1]), FUN="mean")
rownames(mat) <- mat$Group.1
mat <- t(mat[,-1])
head(mat)
mat2 <- mat
mat <- cbind(mat1,mat2)
meta.h <- meta.h[,c(1,2,3,4,7,6,5)]
meta.h <- meta.h[,c(1,3,4,5,6,7,2)]
meta.h$Timepoint[meta.h$Timepoint == "Unknown"] <- "NA"
meta.h$Sex[meta.h$Sex == "Unknown"] <- "NA"
colnames(meta.h) <- c("Experiment.accession","Biosample.organism","Biosample.term.name","Second.term.name","Timepoint","Sex","Biosample.type")
ident1 <- colnames(mat2)
b <- meta.h
b2 <- b[b$"Biosample.term.name" %in% ident1,c(7,3)]
ident2 <- colnames(mat1)
b1 <- b[b$"Second.term.name" %in% ident2,c(7,4)]
colnames(b1) <- c("Biosample.type", "Biosample.term.name")

b3 <- rbind(b1, b2)

b3 <- b3[!duplicated(b3$"Biosample.term.name"),]

b3 <- b3[order(as.character(b3$"Biosample.term.name")),]

b3 <- b3[order(as.character(b3$"Biosample.type"), decreasing= T),]

mat <- mat[,as.character(b3$"Biosample.term.name")]

write.csv(mat, "h.mat.csv")

write.csv(meta.h, "h.meta_all.csv")

####plot figure 3a

data <- read.csv("/Users/ddlove/Desktop/stLFR/数据库/datav8.2/bulk-expression/Single-gene/h.merge_all.csv",row.names = 1)
meta <- read.csv("/Users/ddlove/Desktop/stLFR/数据库/datav8.2/bulk-expression/Single-gene/h.meta_all.csv",row.names = 1)

meta <- meta[meta$Second.term.name %in% c("Hepatocyte_Fetal", "HepG2_Child","Liver_Fetal", "Liver_Adult"),]

data <- data[c("AFP"),rownames(meta)]

data2 <- data.frame(t(data))

data2$class <- meta$Second.term.name

ggplot(data2, aes(class, AFP)) + 
  geom_violin(aes(fill = class)) + coord_flip() +
  geom_boxplot(width = 0.1)





