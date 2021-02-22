##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param path
##' @param mit.cut
##' @return
##' @author dylanmr
##' @export
build_ds <- function(path, mit.cut = 10) {

  plan(sequential)
  dat <- fread(path)
  dat <- as.data.frame(dat) 
  rownames(dat) <- dat[,1]
  dat <- dat[,-1]
  dat <- CreateSeuratObject(dat)
  dat[["percent.mt"]] <- PercentageFeatureSet(dat, pattern = "^mt-") 
  dat <- 
    subset(dat, subset = nFeature_RNA > quantile(nFeature_RNA, .01) & 
             nFeature_RNA < quantile(nFeature_RNA, .99) & 
             percent.mt < mit.cut) %>%
    sctClust(dims=30)
  
  return(dat)

}
