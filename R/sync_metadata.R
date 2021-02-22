##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param path path to googlesheet
##' @return path to filtered file
##' @author dylanmr
##' @export
sync_metadata <- function(path = "https://docs.google.com/spreadsheets/d/1JiIFeY3HomWQy0Jc-yL5JVMmshRUn4Dqxn6mE250mHM/edit#gid=0") {

  # google sheet is public, do not require authorization to access the file
  gs4_deauth()
  # read in the file, specify range to ensure we capture necessary rows x columns
  # read in columns as all characters
  dat <- read_sheet(path, range = "A1:Z1000", col_types = "c")
  # filter data frame and convert columns to correct formats
  dat %>%
    janitor::clean_names() %>%
    dplyr::transmute(
      cage = (cage),
      pup_number = pup_number, 
      sequencing_pool = factor(sequencing_pool),
      age = as.numeric(age),
      date = str_replace_all(collection_date, "\\.", "-") %>%
        lubridate::dmy(),
      weight = str_replace_all(weight_g, ",", ".") %>%
        as.numeric(),
      length = str_replace_all(weight_g, ",", ".") %>%
        as.numeric(length_cm),
      littersize = as.numeric(pre_culled_litter_size),
      leptin = str_replace_all(weight_g, ",", ".") %>%
        as.numeric(elisa_leptin_ng_ml),
      diet = factor(dplyr::if_else(diet == "ob/ob", true = "CHOW", false = diet)),
      geno = factor(dplyr::if_else(grepl("Ob", cage), true = "ob", false = "wt"))
    ) %>% 
    unite("hash_id", cage, pup_number, sep = "-") -> dat

  # return metadata
  return(dat)
}

