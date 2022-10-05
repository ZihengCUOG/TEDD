library(biomaRt)
library(Seurat)

a = list.files(pattern = ".RDS")

count = readRDS(a[1])

httr::set_config(httr::config(ssl_verifypeer = FALSE))

ensembl <- useMart('ensembl', dataset = 'hsapiens_gene_ensembl')
genes_ids <- sub('\\.[0-9]*$', '',rownames(count))
gs_heatdata <- getBM(
  attributes = c('external_gene_name', 'hgnc_symbol','ensembl_gene_id'),
  filters = 'ensembl_gene_id',
  values = genes_ids,
  mart = ensembl)

gs_heatdata

gs_heatdata <- gs_heatdata[!is.na(gs_heatdata$hgnc_symbol),]

gs_heatdata <- gs_heatdata[!duplicated(gs_heatdata$ensembl_gene_id),]

rownames(gs_heatdata) <- gs_heatdata$ensembl_gene_id

genes_ids <- data.frame(genes_ids)

rownames(genes_ids) <- genes_ids$genes_ids

genes_ids$symbol <- gs_heatdata[rownames(genes_ids),"hgnc_symbol"]

Rcpp::sourceCpp(code='
#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
IntegerMatrix asMatrix(NumericVector rp,
                       NumericVector cp,
                       NumericVector z,
                       int nrows,
                       int ncols){

  int k = z.size() ;

  IntegerMatrix  mat(nrows, ncols);

  for (int i = 0; i < k; i++){
      mat(rp[i],cp[i]) = z[i];
  }

  return mat;
}
' )

as_matrix <- function(mat){

  row_pos <- mat@i
  col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])

  tmp <- asMatrix(rp = row_pos, cp = col_pos, z = mat@x,
                 nrows =  mat@Dim[1], ncols = mat@Dim[2])

  row.names(tmp) <- mat@Dimnames[[1]]
  colnames(tmp) <- mat@Dimnames[[2]]
  return(tmp)
}

count <- data.frame(as_matrix(count))

count$id <- genes_ids$symbol

count <- count[!is.na(count$id),]

count <- count[!duplicated(count$id),]

count <- count[!count$id=="",]

rownames(count) <- count$id

count <- count[,-ncol(count)]

meta <- readRDS("meta.Rds")

rownames(meta) <- gsub("-",".",rownames(meta))

colnames(count) <- gsub("-",".",colnames(count))

count <- count[,intersect(rownames(meta),colnames(count))]

meta <- meta[which(rownames(meta) %in% colnames(count)),]

library(Seurat)
Tedd.5 <- CreateSeuratObject(counts = count, project = "Tedd.5",meta.data = meta, min.cells = 1)
Tedd.5[["percent.mt"]] <- PercentageFeatureSet(Tedd.5, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200)
Tedd.5 <- NormalizeData(Tedd.5, normalization.method = "LogNormalize", scale.factor = 10000)
Tedd.5 <- FindVariableFeatures(Tedd.5, selection.method = "vst", nfeatures = 2000)

Tedd.5 <- ScaleData(Tedd.5)
Tedd.5 <- RunPCA(Tedd.5, features = VariableFeatures(object = Tedd.5))
Tedd.5 <- FindNeighbors(Tedd.5, dims = 1:30)
Tedd.5 <- FindClusters(Tedd.5, resolution = 0.5)
Tedd.5 <- RunUMAP(Tedd.5, dims = 1:30)

save(Tedd.5, file = "Tedd.5.Rda")

count <- GetAssayData(object = Tedd.5,slot = "counts")

count <- as_matrix(count)


write.table(count, 
            'counts.csv',sep = ',', row.names = T, col.names = T, quote = F)

embedding <- Tedd.5@reductions$umap@cell.embeddings

meta  <- Tedd.5@meta.data[,c("Tissue","Celltype","Timepoint","Sex")]

meta <- cbind(meta, embedding)

write.table(meta,'meta.csv',sep = ',', row.names = T, col.names = T, quote = F)

Tedd.5@active.ident <- as.factor(Tedd.5$Tissue)
Tissue.mat <- AverageExpression(Tedd.5)
Tissue.mat <- Tissue.mat$RNA

Tedd.5@active.ident <- as.factor(Tedd.5$Celltype)
Celltype.mat <- AverageExpression(Tedd.5)
Celltype.mat <- Celltype.mat$RNA

Tedd.5@active.ident <- as.factor(Tedd.5$Timepoint)
Timepoint.mat <- AverageExpression(Tedd.5)
Timepoint.mat <- Timepoint.mat$RNA

Tedd.5@active.ident <- as.factor(Tedd.5$Sex)
Sex.mat <- AverageExpression(Tedd.5)
Sex.mat <- Sex.mat$RNA

write.csv(Tissue.mat, "Tissue.mat.csv")
write.csv(Celltype.mat, "Celltype.mat.csv")
write.csv(Timepoint.mat, "Timepoint.mat.csv")
write.csv(Sex.mat, "Sex.mat.csv")

