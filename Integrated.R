library(Seurat)
library(glmGamPoi)
library(harmony)

load("Tedd.4_AdultAdrenalGland/Tedd.4.Rda")

Tedd.4.1 <- SCTransform(Tedd.4, verbose = FALSE)

Tedd.4$stim <- "Tedd.4"

load("Tedd.4_FetalAdrenalGland/Tedd.4.Rda")

Tedd.4.2 <- SCTransform(Tedd.4, verbose = FALSE)

Tedd.4.2$stim <- "Tedd.4"

load("Tedd.4_NeonatalAdrenalGland/Tedd.4.Rda")

Tedd.4.3 <- SCTransform(Tedd.4, verbose = FALSE,method = "glmGamPoi")

Tedd.4.3$stim <- "Tedd.4"

load("Tedd.5_FetalAdrenalGland/Tedd.5.Rda")

Tedd.5 <- SCTransform(Tedd.5, verbose = FALSE,method = "glmGamPoi")

Tedd.5$stim <- "Tedd.5"

gene <- Reduce(intersect,  list(Tedd.4.1 = rownames(Tedd.4.1), Tedd.4.2 = rownames(Tedd.4.2),Tedd.4.3 = rownames(Tedd.4.3), Tedd.5 = rownames(Tedd.5)))

save(gene, file = "gene.Rda")

object.list <- list(Tedd.4 = Tedd.4.1, Tedd.4.2 = Tedd.4.2, Tedd.4.3 = Tedd.4.3,Tedd.5 = Tedd.5)

features <- SelectIntegrationFeatures(object.list = object.list, nfeatures = 3000)

features

merged.data <- merge(Tedd.4.1, c(Tedd.4.2,Tedd.4.3,Tedd.5),project = "Implantation")

VariableFeatures(merged.data) <- features
merged.data <- RunPCA(object = merged.data, assay = "SCT", npcs = 30)
merged.data <- RunHarmony(object = merged.data,
                                    assay.use = "SCT",
                                    reduction = "pca",
                                    dims.use = 1:30,
                                    group.by.vars = "stim",
                                    plot_convergence = TRUE)


merged.data <- RunUMAP(object = merged.data, assay = "SCT", reduction = "harmony", dims = 1:30)
merged.data <- FindNeighbors(object = merged.data, assay = "SCT", reduction = "harmony", dims = 1:30)
merged.data <- FindClusters(object = merged.data, resolution = 0.5)

Tedd.10 <- merged.data

load("meta.Rda")

meta <- meta[rownames(Tedd.10@meta.data),]

Tedd.10@meta.data <- meta

save(Tedd.10, file = "Tedd.10.Rda")

DefaultAssay(Tedd.10) <- "SCT"

count <- GetAssayData(object = Tedd.10, slot = "counts",assay = "SCT")

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

count <- as_matrix(count)

write.table(count,
            'counts.csv',sep = ',', row.names = T, col.names = T, quote = F)

embeddings <- Tedd.10@reductions$umap@cell.embeddings

meta  <- Tedd.10@meta.data[,c("Tissue","Celltype","Timepoint","Sex")]

meta <- cbind(meta, embeddings)

write.table(meta,'meta.csv',sep = ',', row.names = T, col.names = T, quote = F)

Tedd.10@active.ident <- as.factor(Tedd.10$Tissue)
Tissue.mat <- AverageExpression(Tedd.10)
Tissue.mat <- Tissue.mat$SCT

Tedd.10@active.ident <- as.factor(Tedd.10$Celltype)
Celltype.mat <- AverageExpression(Tedd.10)
Celltype.mat <- Celltype.mat$SCT


Tedd.10@active.ident <- as.factor(Tedd.10$Timepoint)
Timepoint.mat <- AverageExpression(Tedd.10)
Timepoint.mat <- Timepoint.mat$SCT

Tedd.10@active.ident <- as.factor(Tedd.10$Sex)
Sex.mat <- AverageExpression(Tedd.10)
Sex.mat <- Sex.mat$SCT

write.csv(Tissue.mat, "Tissue.mat.csv")
write.csv(Celltype.mat, "Celltype.mat.csv")
write.csv(Timepoint.mat, "Timepoint.mat.csv")
write.csv(Sex.mat, "Sex.mat.csv")

###plot figure 3b
library(Seurat)
load("Tedd.10_Liver_scrna/Tedd.10.Rda")
DimPlot(Tedd.10, group.by = "Celltype")
DimPlot(Tedd.10, group.by = "Timepoint")
FeaturePlot(Tedd.10, features = "AFP")

library(Seurat)
load("Tedd.10_Liver_scatac/Tedd.10.Rda")
DimPlot(Tedd.10, group.by = "Celltype")
DimPlot(Tedd.10, group.by = "Timepoint")
DimPlot(Tedd.10, group.by = "Timepoint")
FeaturePlot(Tedd.10, features = "AFP")

