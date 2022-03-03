library(targets)
library(tarchetypes)
library(future)
library(future.callr)

options(tidyverse.quiet = TRUE)
plan(callr)

tar_option_set(packages = c(
  "mclust", "Seurat", "tidyverse", "tidymodels", "yardstick", "R2HTML",
  "DropletUtils", "Matrix", "rjson", "stringr", "future", "doParallel",
  "org.Mm.eg.db", "scDblFinder", "AnnotationDbi", "SingleCellExperiment",
  "scran", "scater", "glmGamPoi", "miloR", "BiocNeighbors", "data.table",
  "magrittr")
)

dev <- function(obj, ageset) {
  CreateSeuratObject(counts = counts(obj), meta.data = colData(obj) %>% as.data.frame) %>% 
    subset(., subset = age < ageset) %>% 
    sctClust(dims=40) %>% 
    mapCamp() %>% 
    classifySex()
}

lep <- function(obj, ageset) {
  CreateSeuratObject(counts = counts(obj), meta.data = colData(obj) %>% as.data.frame) %>% 
    subset(., subset = age > ageset & treatment != "NA") %>% 
    sctClust(dims=40) %>% 
    mapCamp() %>% 
    classifySex()
}

values <- tibble::tibble(
  method_function = rlang::syms(c("dev", "lep"))
  )

## Load your R files
lapply(list.files("./R", full.names = TRUE), source)

preprocessed <- 
  list(
    tarchetypes::tar_files(
      hash,
      path_to_data(proj_num = "SCOP_37", type = "kite")
    ),
    tarchetypes::tar_files(
      exp,
      path_to_data(proj_num = "SCOP_37", type = "kallisto")
    ),
    tar_target(
      metapath, 
      path_to_meta(), 
      format="file"
    ),
    tar_target(
      meta, 
      readxl::read_excel(metapath) %>% 
        janitor::clean_names() %>% 
        mutate(kitv = factor(x10x_version))
    ),
    tar_target(
      web_summaries, 
      gen_summary(files = exp), 
      format = "file", 
      pattern = map(exp)
    ),
    tar_target(
      sce_objects, 
      build_sce_w_hto(exp, hash), 
      pattern = map(exp, hash), 
      iteration = "list"
    ),
    tar_target(
      param_df, 
      tibble::tibble(par = c(1:44)) %>% 
        dplyr::mutate(output_file = here::here(paste0("output/qc/hto_demux_qc", 1:n(), ".html")))
    ),  
    #tar_render_rep(
    #  rendered_report,
    #  "qc_report.Rmd",
    #  params = param_df
   # ),
    tar_target(
      mergedsce,
      comb_and_filter(sce_objects, meta, exp)
    )
  )
    
seuratobjects <- 
  tar_map(
    values = values,
    tar_target(
      fulldata,
      method_function(mergedsce, ageset = 35)
    ),
    tar_target(
      neuron, 
      subset(fulldata, subset = predicted.id == "neuron") %>% 
        sctClust(dims=40) %>% 
        mapCamp(class="neuron")
    ),
    tar_target(
      glia, 
      subset(fulldata, subset = predicted.id != "neuron") %>% 
        sctClust(dims=40) %>% 
        mapCamp(class="glia")
    ),
    tar_target(
      milo_neuron,
      gen_milo(neuron)
    ),
    tar_target(
      milo_glia,
      gen_milo(glia)
    )
  )

substates <- 
  list(
    tar_target(
      milo_agrp,
      subset(neuron_dev, idents = c(1,6)) %>% 
        sctClust(dims=20) %>% 
        gen_milo(., k=30, d=20)
    ),
    tar_target(
      milo_pomc,
      subset(neuron_dev, idents = c(4)) %>% 
        sctClust(dims=20) %>% 
        gen_milo(., k=30, d=20)
    )
  )

modeling <-
  list(
    tar_target(
      split_neuron,
      splitall(neuron_dev, grantdata = grantdata),
      iteration="list"
    ),
    tar_target(
      paramsformodel,
      paramsweep(split_neuron),
      pattern = map(split_neuron),
      iteration="list"
    ),
    tar_target(
      grantdata,
      qs::qread("data/grantdata.qs")
    ),
    tar_target(
      neuronpreds,
      predictage(devdata = split_neuron, lepdata = neuron_lep, grantdata = grantdata, lasso_params = paramsformodel),
      pattern = map(split_neuron, paramsformodel),
      iteration="list"
    )
  )


list(preprocessed, seuratobjects)#, modeling, substates)

  # tar_target(
  #   seurbycluster,
  #   splitall(neuron, glia),
  #   iteration="list"
  # ),
  # tar_target(
  #   lasso_params,
  #   paramsweep(seurbycluster),
  #   pattern = map(seurbycluster)
  # ),
  #  tar_target(
  #    predicted_age,
  #    predage(seur, lasso_params),
  #    pattern = map(seurbycluster, lasso_params)
  #  )


  # ),
  # 
  # tar_target(
  #   milo, 
  #   gen_milo(x = split_seur, k=30, d=30),
  #   pattern = map(split_seur), 
  #   iteration="list"
  # ),
  # tar_target(
  #   pseudobulk,
  #   Seurat::as.SingleCellExperiment(split_seur[[1]], assay="SCT") %>%
  #     scater::aggregateAcrossCells(., use.assay.type="counts", id=colData(.)[,c("predicted.id", "hash_id")]),
  #   pattern = map(split_seur), 
  #   iteration="list"
  # ),
  # tar_target(
  #   knn_graph,
  #   build_graph(split_seur[[1]]),
  #   pattern = map(split_seur), 
  #   iteration="list"
  # ),
  # tar_target(
  #   resampled_graph,
  #   knn_resampling(knn_graph, k = 30, meta = split_seur[[1]][[]]),
  #   pattern = map(knn_graph, split_seur),
  #   iteration = "list"
  # )
  # tar_target(
  #   integ_with_embr,
  #   build_ds("/data/pub-others/huisman-natcom-2019/GSE126480_Raw_counts.txt.gz") %>% 
  #     integrate_data(sc_comb, ., label = "huisman")
  # )
  #),
  # tar_map(
  #   values = clustering_params,
  #   names = "name", # Select columns from `values` for target names.
  #   tar_target(
  #     clusters,
  #     cluster_across_params(sce_final, dims = dims, k = k, weight = weight, clus_method = clus_method)
  #   )
  #)


