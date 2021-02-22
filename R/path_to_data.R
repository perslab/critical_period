##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param proj_num project number to find paths
##' @param type which type of data is being located
##' @return paths to data
##' @author dylanmr
##' @export

path_to_data <- function(proj_num = "SCOP_37", type = NULL) {

  # get type (gene exp or hto data)
  type <- if_else(is.null(type), true = "kallisto", false = "kite")
  # glob all file names in the project folders (do not end in h5ad)
  x <- Sys.glob(paste0("/", file.path("data", "sc-10x", "data-kallisto", proj_num, type, "*", "counts_unfiltered", "*"))) %>%
    str_subset("h5ad", negate = T)
  if (type == "kallisto") {
    # there are 4 files that go with each hash pool, we need to batch each of these together
    # get the hash pool name
    y <- stringr::str_match_all(x, pattern = ".*kallisto\\/(.*)\\/counts.*") %>%
      do.call(rbind, .) %>%
      .[, 2]
    # create a dataframe where each row is an individual file, column 2 has the hash pool that file is associated with0
    data <- cbind(x, y)
    # create a list based on column 2, has the hash pool identity
    tapply(data[, 1], INDEX = data[, 2], identity, simplify = FALSE) %>% unname()
  } else {
    # there are 3 files that go with each hash pool, we need to batch each of these together
    # get the hash pool name
    y <- stringr::str_match_all(x, pattern = ".*kite\\/(.*)\\/counts.*") %>%
      do.call(rbind, .) %>%
      .[, 2]
    # create a dataframe where each row is an individual file, column 2 has the hash pool that file is associated with
    data <- cbind(x, y)
    # create a list based on column 2, has the hash pool identity
    tapply(data[, 1], INDEX = data[, 2], identity, simplify = FALSE) %>% unname()
  }
}
