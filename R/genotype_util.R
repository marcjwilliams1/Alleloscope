#' Map genotypes to allele specific states
#'
#' @param Obj_filtered An Alleloscope object with a n cell by (m region * 4) genotype_values matrix.
#' Every 2 columns in the genotype_table matrix are (rho_hat, theta_hat) of each region.

#' @return An alleloscope object with an additional slot "genotype_states" that has the from M|m where M is the 
#' major allele copy number and m is the minor allele copy number
#' 
#' @export
map_genotypes <- function(Obj_filtered){
  # check parameters
  if(is.null(Obj_filtered)){
    stop("Please provide a valid Alleloscope object for Obj_filtered.")
  }
  
  mapping <- c("1" = "0|1",
               "2" = "1|0",
               "3" = "0|2",
               "4" = "1|1",
               "5" = "2|0",
               "6" = "0|3",
               "7" = "1|2",
               "8" = "2|1",
               "9" = "3|0",
               "10" = "0|4",
               "11" = "1|3",
               "12" = "2|2",
               "13" = "3|1",
               "14" = "4|0",
               "15" = "0|5",
               "16" = "1|5",
               "17" = "2|3",
               "18" = "3|2",
               "19" = "4|1",
               "20" = "5|0",
               "21" = "0|6",
               "22" = "1|5",
               "23" = "2|4",
               "24" = "3|3",
               "25" = "4|2",
               "26" = "5|1",
               "27" = "6|0")
  
  genotypes_states <- unlist(lapply(paste0(Obj_filtered$genotypes), function(x) mapping[[x]]))
  genotypes_states_mat <- matrix(genotypes_states, 
                                 nrow = nrow(Obj_filtered$genotypes), 
                                 ncol = ncol(Obj_filtered$genotypes))
  colnames(genotypes_states_mat) <- colnames(Obj_filtered$genotypes)
  rownames(genotypes_states_mat) <- rownames(Obj_filtered$genotypes)
  Obj_filtered$genotypes_states <- genotypes_states_mat
  
  return(Obj_filtered)
}

long_data_frame_helper <- function(x, colname1 = NULL, colname2 = NULL){
  df <- as.data.frame(x)
  df$cell_id <- row.names(df)
  df <- tidyr::pivot_longer(df, -cell_id, names_to = colname1, values_to = colname2)
  return(df)
}

#' Creat long dataframe with allele specific states and rho and theta values
#'
#' @param Obj_filtered An Alleloscope object
#' 
#' @return dataframe with allele specific states and rho and theta values extracted from the alleloscope object
#' 
#' @export
create_long_data_frame <- function(Obj_filtered){
  
  if (!"genotypes_states" %in% names(Obj_filtered)){
    Obj_filtered <- map_genotypes(Obj_filtered)
  }
  
  #get states dataframe
  states <- long_data_frame_helper(Obj_filtered$genotypes_states, colname1 = "segid", colname2 = "cn_AS")
  #add total and Maj, Min columns
  states <- tidyr::separate(states, "cn_AS", c("Maj", "Min"), "\\|", remove = FALSE)
  states$Maj <- as.numeric(states$Maj)
  states$Min <- as.numeric(states$Min)
  states$cn_total <- states$Maj + states$Min
  
  
  rawdata <- long_data_frame_helper(Obj_filtered$genotype_values, colname1 = "segid_param", colname2 = "value")
  rawdata <- tidyr::separate(rawdata, segid_param, c("param", "segid"), "_")
  rawdata <- dplyr::filter(rawdata, param %in% c("rho", "theta"))
  rawdata <- tidyr::pivot_wider(rawdata, names_from = "param", values_from = "value")
  
  df <- inner_join(states, rawdata, by = c("cell_id", "segid"))
  df <- tidyr::separate(df, "segid", c("chr", "start"), ":")
  df <- dplyr::left_join(df, dplyr::select(Obj_filtered$seg_table_filtered, chr, start, end), by = c("chr", "start"))
  df <- dplyr::select(df, cell_id, chr, start, end, cn_AS, cn_total, Maj, Min, rho, theta)
  df <- dplyr::arrange(df, cell_id, chr, start)
  
  df$start <- as.numeric(df$start)
  df$end <- as.numeric(df$end)
  
  return(as.data.frame(df))
}



