setwd("/data/Chromosome")

b <- read.csv("m.differ.gene.csv")

gene <- read.table("mouse.ident.txt")

gene <- gene[!duplicated(gene$V5),]

rownames(gene) <- gene$V5

b$chromosome <- gene[b$gene, "V1"]
b$start <- gene[b$gene, "V2"]
b$end <- gene[b$gene, "V3"]
b$orientation <- gene[b$gene, "V4"]
b$gene_type <- gene[b$gene, "V6"]

b <- b[!is.na(b$chromosome),]

write.csv(b, file = "m.gene.csv",quote = F, row.names = F)
