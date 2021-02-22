##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param sce_object
##' @return
##' @author dylanmr
##' @export

sceSCT <- function(sce_object, n_genes=3000) {

  # store sctransform output in a sce object
  scnorm <- sctransform::vst(counts(sce_object), return_cell_attr = T, n_cells = 20000, return_corrected_umi = T, method = "qpoisson")
  # log1p of corrected umi
  log_corrected <- log1p(scnorm$umi_corrected)
  # matrix of highly variable genes by pearson residuals
  feature.variance <- setNames(object = scnorm$gene_attr$residual_variance, nm = rownames(x = scnorm$gene_attr))
  feature.variance <- sort(x = feature.variance, decreasing = TRUE)
  top.features <- names(x = feature.variance)[1:min(n_genes, length(x = feature.variance))]
  scale.data <- scnorm$y[top.features, ]
  # clip the residuals
  clip.range <- c(-sqrt(ncol(sce_object) / 30), sqrt(ncol(sce_object) / 30))
  scale.data[scale.data < clip.range[1]] <- clip.range[1]
  scale.data[scale.data > clip.range[2]] <- clip.range[2]
  # 2nd regression
  scale.data <- Seurat::ScaleData(
    scale.data,
    features = NULL,
    model.use = 'linear',
    use.umi = FALSE,
    do.scale = F,
    do.center = T,
    scale.max = Inf,
    block.size = 750,
    min.cells.to.block = 3000
  )
  # filter sce object to have same rownames as corrected umi matrix
  sce_object <- sce_object[rownames(scnorm$umi_corrected),]
  # replace counts assay with corrected umi
  assay(sce_object, "counts") <- scnorm$umi_corrected
  # replace logcounts with log1p of corrected counts
  assay(sce_object, "logcounts") <- log_corrected
  # create an altExp with scaled highly variable genes
  altExp(sce_object, "hvg") <- SingleCellExperiment(list(counts=scale.data))
  # swap default experiment
  updated_exp <- swapAltExp(sce_object, name = "hvg", saved = "original")
  # make counts and logcounts identical
  assay(updated_exp, "logcounts") <- assay(updated_exp, "counts")
  # return object
  return(updated_exp)

}
