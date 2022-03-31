library(targets)
library(tarchetypes)
library(future)
library(future.callr)

options(tidyverse.quiet = TRUE)
plan(callr)

packages = c(
  "mclust", "Seurat", "tidyverse", "tidymodels", "yardstick", "R2HTML",
  "DropletUtils", "Matrix", "rjson", "stringr", "future", "doParallel",
  "org.Mm.eg.db", "scDblFinder", "AnnotationDbi", "SingleCellExperiment",
  "scran", "scater", "glmGamPoi", "miloR", "BiocNeighbors", "data.table",
  "magrittr", "purrr", "edgeR", "cacoa")

vapply(packages, library, logical(1L), character.only = TRUE, logical.return = TRUE)

values <- tibble::tibble(
  dataset = c("dev", "lepip", "lepan", "agrpfl"),
  kitv = c(list(1),list(c(1,2)),list(c(1,2)),list(c(1,2))),
  #method_function = rlang::syms(c("dev", "lepip", "lepan", "agrpfl")),
  trt_grouping = c("geno", "treatment", "treatment", "treatment"),
  sample_column = rep("hash",4), 
  cell_grouping = c(list(c("SCT_snn_res.0.8","geno","age")),list(c("SCT_snn_res.0.8","geno","time")), list("SCT_snn_res.0.8"), list("SCT_snn_res.0.8")),
  design = c(list(c("geno", "diet", "age")), list(c("geno", "time", "treatment")), list('treatment'), list('geno')),
  batch = c("kitv", "seq_pool", "hash_pool", "hash_pool"),
  names = c("dev", "lepip", "lepan", "agrpfl")
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
    names = "names",
    tar_target(
      fulldata,
      CreateSeuratObject(counts = counts(fulldata), meta.data = colData(fulldata) %>% as.data.frame) %>% 
        subset(., subset = dataset == dataset & kitv %in% kitv) %>% 
        sctClust(dims=40) %>% 
        mapCamp()
    ),
    tar_target(
      cacoa,
      gen_cacoa_obj(fulldata, trt_grouping = trt_grouping, sample_column = sample_column, cell_grouping = cell_grouping)
    ),
    tar_target(
        neuron,
        subset(fulldata, subset = predicted.id == "neuron") %>%
          sctClust(dims=40, quant=0.05) %>%
          mapCamp(class="neuron")
      ),
    tar_target(
      glia,
      subset(fulldata, subset = predicted.id != "neuron") %>%
        sctClust(dims=40, quant=0.05) %>%
        mapCamp(class="glia")
    ),
    tar_target(
      milo_neuron,
      gen_milo(neuron)
    ),
    tar_target(
      milo_glia,
      gen_milo(glia)
    ),
    tar_target(
      neuron_edger,
      splitwrapper(neuron, split.by="SCT_snn_res.0.8") %>% 
        map(~build_edger(.x, 10, design = c("geno", "time", "treatment"), batch = c("seq_pool")))
    ),
    tar_target(
      glia_edger,
      splitwrapper(glia) %>% 
        map(~build_edger(.x, 10, design = c("geno", "time", "treatment"), batch = c("seq_pool")))
    )
  )


# degenes <-
#   list(
#     #LEPIP
#     tar_target(
#       splitneur_lepip,
#       splitwrapper(neuron_lepip, split.by="SCT_snn_res.0.8"),
#       iteration="list"
#     ),
#     tar_target(
#       splitglia_lepip,
#       splitwrapper(glia_lepip),
#       iteration="list"
#     ),
#     tar_target(
#       neuron_lep_edger,
#       build_edger(splitneur_lepip, 10, design = c("geno", "time", "treatment"), batch = c("seq_pool")),
#       pattern = map(splitneur_lepip)
#     ),
#     tar_target(
#       glia_lep_edger,
#       build_edger(splitglia_lepip, 10, design = c("geno", "time", "treatment"), batch = c("seq_pool")),
#       pattern = map(splitglia_lepip)
#     ),
#     #LEP ANT
#     tar_target(
#       splitneur_lepan,
#       splitwrapper(neuron_lepan, split.by="SCT_snn_res.0.8"),
#       iteration="list"
#     ),
#     tar_target(
#       splitglia_lepan,
#       splitwrapper(glia_lepan),
#       iteration="list"
#     ),
#     tar_target(
#       neuron_lepan_edger,
#       build_edger(splitneur_lepan, 10, design = c("treatment"), batch = c("hash_pool")),
#       pattern = map(splitneur_lepan)
#     ),
#     tar_target(
#       glia_lepan_edger,
#       build_edger(splitglia_lepan, 10, design = c("treatment"), batch = c("hash_pool")),
#       pattern = map(splitglia_lepan)
#     ),
#     #AGRP FL
#     tar_target(
#       splitneur_agrpfl,
#       splitwrapper(neuron_agrpfl, split.by="SCT_snn_res.0.8"),
#       iteration="list"
#     ),
#     tar_target(
#       splitglia_agrpfl,
#       splitwrapper(glia_agrpfl),
#       iteration="list"
#     ),
#     tar_target(
#       neuron_agrpfl_edger,
#       build_edger(splitneur_agrpfl, 10, design = c("geno"), batch = c("hash_pool")),
#       pattern = map(splitneur_agrpfl)
#     ),
#     tar_target(
#       glia_agrpfl_edger,
#       build_edger(splitglia_agrpfl, 10, design = c("geno"), batch = c("hash_pool")),
#       pattern = map(splitglia_agrpfl)
#     ),
#     #DEV
#     tar_target(
#       splitneur_dev,
#       splitwrapper(neuron_dev),
#       iteration="list"
#     ),
#     tar_target(
#       splitglia_dev,
#       splitwrapper(glia_dev),
#       iteration="list"
#     ),
#     tar_target(
#       neur_dev_edger,
#       build_edger(splitneur_dev, ncell = 10, design = c("geno", "diet", "age"), batch = c("kitv")),
#       pattern = map(splitneur_dev)
#     ),
#     tar_target(
#       glia_dev_edger,
#       build_edger(splitglia_dev, ncell = 10, design = c("geno", "diet", "age"), batch = c("kitv")),
#       pattern = map(splitglia_dev)
#     )
#   )
# 
# cacoa_obj <-
#   list(
#     tar_target(
#       cacoa_neuron_lepip, 
#       gen_cacoa_obj(neuron_lepip, trt_grouping = "treatment", sample_column = "hash", cell_grouping = c("SCT_snn_res.0.8","geno","time"))
#     ),
#     tar_target(
#       cacoa_neuron_lepan, 
#       gen_cacoa_obj(neuron_lepan, trt_grouping = "treatment", sample_column = "hash", cell_grouping = c("SCT_snn_res.0.8"))
#     ),
#     tar_target(
#       cacoa_neuron_agrpfl,
#       gen_cacoa_obj(neuron_agrpfl, trt_grouping = "geno", sample_column = "hash", cell_grouping = c("SCT_snn_res.0.8"))
#     ),
#     tar_target(
#       cacoa_neuron_dev,
#       gen_cacoa_obj(neuron_dev, trt_grouping = "geno", sample_column = "hash", cell_grouping = c("SCT_snn_res.0.8"))
#     )
#   )


substates <- 
  list(
    tar_target(
      milo_agrp,
      subset(neuron_dev, subset = predicted.id == "Agrp") %>% 
        sctClust(dims=20) %>% 
        gen_milo(., k=30, d=20)
    ),
    tar_target(
      milo_pomc,
      subset(neuron_dev, subset = predicted.id == "Pomc_Lepr") %>% 
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


list(preprocessed, seuratobjects)#, degenes, substates, cacoa_obj, cacoa_testing)#, modeling, substates)

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


# dev <- function(obj, ageset) {
#   CreateSeuratObject(counts = counts(obj), meta.data = colData(obj) %>% as.data.frame) %>% 
#     subset(., subset = age <= ageset & kitv == 1) %>% 
#     sctClust(dims=40) %>% 
#     mapCamp() %>% 
#     classifySex()
# }
# 
# lepip <- function(obj, ageset) {
#   CreateSeuratObject(counts = counts(obj), meta.data = colData(obj) %>% as.data.frame) %>% 
#     subset(., subset = age > ageset & treatment %in% c("Lep","Sal")) %>% 
#     sctClust(dims=40) %>% 
#     mapCamp() %>% 
#     classifySex()
# }
# 
# lepan <- function(obj, ageset) {
#   CreateSeuratObject(counts = counts(obj), meta.data = colData(obj) %>% as.data.frame) %>% 
#     subset(., subset = age > ageset & treatment %in% c("ANT", "CON")) %>% 
#     sctClust(dims=40) %>% 
#     mapCamp() 
# }
# 
# agrpfl <- function(obj, ageset) {
#   CreateSeuratObject(counts = counts(obj), meta.data = colData(obj) %>% as.data.frame) %>% 
#     subset(., subset = age == ageset & kitv == 2) %>% 
#     sctClust(dims=40) %>% 
#     mapCamp() %>% 
#     classifySex()
# }


