library(Seurat)

setwd("/data/Integrate")

load("Tedd.4_NeonatalAdrenalGland/Tedd.4.Rda")

Tedd.4.3 <- SCTransform(Tedd.4, verbose = FALSE,method = "glmGamPoi")

Tedd.4.3$stim <- "Tedd.4"

load("Tedd.5_FetalAdrenalGland/Tedd.5.Rda")

Tedd.5 <- SCTransform(Tedd.5, verbose = FALSE,method = "glmGamPoi")

Tedd.5$stim <- "Tedd.5"

merge.data <- merge(Tedd.4.3, c(Tedd.5))

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

write.csv(data, "Multi-gene/h.multigene.csv")

Tissue.mat <- AverageExpression(merge.data,slot = "data",group.by = "Tissue")
Tissue.mat <- Tissue.mat$SCT

Timepoint.mat <- AverageExpression(merge.data,slot = "data",group.by = "Timepoint")
Timepoint.mat <- Timepoint.mat$SCT

data <- cbind(Tissue.mat, Timepoint.mat)

write.csv(data, "Stable/h.mat.csv")

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

write.table(count,'Single-gene/h.counts.csv',sep = ',', row.names = T, col.names = T, quote = F)

