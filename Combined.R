library(Seurat)

Tedd.11 <- readRDS("../../Tedd.11/save2.Rds")
Tedd.11 <- SCTransform(Tedd.11, verbose = FALSE, method = "glmGamPoi")
Tedd.11$stim <- "Tedd.11"

Tedd.12 <- readRDS("../../Tedd.12_FetalLimb/save2.Rds")
Tedd.12 <- SCTransform(Tedd.12, verbose = FALSE,method = "glmGamPoi")
Tedd.12$stim <- "Tedd.12"

Tedd.13 <- readRDS("../../Tedd.13/save2.Rds")
Tedd.13 <- SCTransform(Tedd.13, verbose = FALSE,method = "glmGamPoi")
Tedd.13$stim <- "Tedd.13"

Tedd.14 <- readRDS("../../Tedd.14/save2.Rds")
Tedd.14 <- SCTransform(Tedd.14, verbose = FALSE,method = "glmGamPoi")
Tedd.14$stim <- "Tedd.14"

Tedd.15 <- readRDS("../../Tedd.15/save2.Rds")
Tedd.15 <- SCTransform(Tedd.15, verbose = FALSE, method = "glmGamPoi")
Tedd.15$stim <- "Tedd.15"

Tedd.17 <- readRDS("../../Tedd.17/save2.Rds")
Tedd.17 <- SCTransform(Tedd.17, verbose = FALSE,method = "glmGamPoi")
Tedd.17$stim <- "Tedd.17"

Tedd.18 <- readRDS("../../Tedd.18/save2.Rds")
Tedd.18 <- SCTransform(Tedd.18, verbose = FALSE,method = "glmGamPoi")
Tedd.18$stim <- "Tedd.18"

Tedd.19 <- readRDS("../../Tedd.19_Pancreas/save2.Rds")
Tedd.19 <- SCTransform(Tedd.19, verbose = FALSE,method = "glmGamPoi")
Tedd.19$stim <- "Tedd.19"

Tedd.23 <- readRDS("../../Tedd.23_Testis/save2.Rds")
Tedd.23 <- SCTransform(Tedd.23, verbose = FALSE,method = "glmGamPoi")
Tedd.23$stim <- "Tedd.23"

merge.data <- merge(Tedd.11, c(Tedd.12, Tedd.13, Tedd.14, Tedd.15, Tedd.17, Tedd.18, Tedd.19, Tedd.23),project = "mouse")

genes.to.keep <- Matrix::rowSums(merge.data@assays$RNA@counts > 0) >= 10

gene <- rownames(merge.data@assays$RNA@counts[genes.to.keep,])

merge.data <- subset(merge.data, features = gene)

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

saveRDS(merge.data, file = "save2.Rds")

merge.data@meta.data$Tissue_Timepoint <- paste(merge.data@meta.data$Tissue, "_","All","_",merge.data@meta.data$Timepoint, sep = "")
merge.data@active.ident <- as.factor(merge.data$Tissue_Timepoint)
Timepoint.mat <- AverageExpression(merge.data,slot = "data",group.by = "Tissue_Timepoint")
Timepoint.mat <- Timepoint.mat$SCT

merge.data@meta.data$Tissue_Celltype_Timepoint <- paste(merge.data@meta.data$Tissue, "_",merge.data@meta.data$Celltype, "_",merge.data@meta.data$Timepoint, sep = "")
merge.data@active.ident <- as.factor(merge.data$Tissue_Celltype_Timepoint)
Celltype.mat <- AverageExpression(merge.data,slot = "data",group.by = "Tissue_Celltype_Timepoint")
Celltype.mat <- Celltype.mat$SCT

data <- Reduce(cbind, list(Timepoint.mat, Celltype.mat))

write.csv(data, "Multi-gene/m.multigene.csv")

Tissue.mat <- AverageExpression(merge.data,slot = "data",group.by = "Tissue")
Tissue.mat <- Tissue.mat$SCT

Timepoint.mat <- AverageExpression(merge.data,slot = "data",group.by = "Timepoint")
Timepoint.mat <- Timepoint.mat$SCT

data <- cbind(Tissue.mat, Timepoint.mat)

write.csv(data, "Stable/m.mat.csv")

count <- GetAssayData(object = merge.data, slot = "data",assay = "SCT")

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

write.table(count,'Single-gene/m.counts.csv',sep = ',', row.names = T, col.names = T, quote = F)

