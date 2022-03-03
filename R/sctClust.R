##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title

##' @return
##' @author dylanmr
##' @export
##' 


sctClust <- function(x, dims=NULL, resolution=0.8, arguments=NULL) {
  
  future::plan(future::sequential)
  
  if(is(x, "SingleCellExperiment")) {
    x <- CreateSeuratObject(counts = counts(x), meta.data = colData(x) %>% as.data.frame)
  } 
  
  DefaultAssay(x) <- "RNA"
  
  if(is.null(arguments)) {
    x %>% 
      SCTransform(method = "qpoisson", vst.flavor = "v2") %>% 
      RunPCA(assay = "SCT", reduction.name = "pcsct", reduction.key = "pcsct_") -> x
    
    if(is.null(dims)) {
      dims <- round(as.numeric(intrinsicDimension::maxLikGlobalDimEst(data = x@reductions$pcsct[, 1:50], k = 20)))
    } else {
      dims <- dims
    }
    
    x %>% 
      RunUMAP(dims = seq(dims), reduction="pcsct", reduction.name = "umapsct", reduction.key = "umapsct_") %>%
      FindNeighbors(dims = seq(dims), reduction="pcsct", force.recalc=T) %>%
      FindClusters(resolution = resolution) -> x
    
    return(x)
    
  } else {
    
    x %>% 
      SCTransform(method = "qpoisson") %>% 
      .sc_correct_batch(arguments = arguments) %>% 
      RunPCA(assay="limma", reduction.name="pclimma", reduction.key = "pclimma_") %>% 
      RunPCA(assay = "SCT", reduction.name = "pcsct", reduction.key = "pcsct_") -> x
    
    if(is.null(dims)) {
      dims <- round(as.numeric(intrinsicDimension::maxLikGlobalDimEst(data = x@reductions$pclimma[, 1:50], k = 20)))
    } else {
      dims <- dims
    }
    
    x %>% 
      RunUMAP(dims = seq(dims), reduction="pclimma", reduction.name = "umaplimma", reduction.key = "umaplimma_") %>%
      RunUMAP(dims = seq(dims), reduction="pcsct", reduction.name = "umapsct", reduction.key = "umapsct_") %>%
      FindNeighbors(dims = seq(dims), reduction="pclimma", force.recalc=T) %>%
      FindClusters(resolution = resolution, graph.name = "limma_snn") -> x
    
    return(x)
  }
}

.sc_correct_batch <- function(seur, arguments) {
  design <- formula(paste("~", paste(arguments, collapse = "+"), collapse = " "))
  data <- seur@assays$SCT@scale.data
  mm <- model.matrix(design, seur[[]])
  mat <- limma::removeBatchEffect(data, batch=seur$kitv, design=mm)
  seur[["limma"]] <- seur@assays$SCT
  seur[["limma"]]@scale.data <- mat
  return(seur)
}

