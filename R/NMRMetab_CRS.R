###############################################################################################
#This function will take in either a object of class data.frame and a using a second dataframe#
#with a list of all the metabolite identifed, will return the CRS score for all the bins.     #
#this could work with other metabolite identifiers; for example bins with the name of the     #
#metabolite, as  log as they are unique. mistakes could happen when greple take NAD and NADH  #
#INPUTS:                                                                                      #
# df = the dataframe with the column names as HMDBs and no other columns descriptors such as  #
#      group identifier or sample names. column names have to be 'Metab' for the metabolite names
# 	and 'HMDB' with the corresponding HMDB name                                           #
# metabolite_list = a new dataframe with 2 columns. A column called 'HMDB' with the HMDB      #
#      identifier and another colunn called 'Match' with the names of the metabolites         #
# identifier = this is set to 'Match' as a default. This argument will use the 'Match' column #
#      to identify the metabolites in the output                                              #
#
#OUTPUTS:
# metabdf = filtered dataframe fof the specific metabolite
# corrmatrix = correlation matrices for the corresponding filtered dataframe
# averages = the crs score calculated from the corrmatrix
# CRS_pass = the crs pass score for the whole df. follows rules set out by Rudi in its thesis
# frequency_metabs = table of frequency of the metabolites in the original dataframe
# frequency_plot = plot of the frequency table.
###############################################################################################


#' @title CRS function
#' @name NMRMetab_CRS
#' @export
#' @description This function will take in either a object of class data.frame and a using a second dataframe with a list of all the metabolite identifed, will return the CRS score for all the bins. this could work with other metabolite identifiers; for example bins with the name of the     metabolite, as  log as they are unique. mistakes could happen when greple take NAD and NADH
#' @author Michele Fresneda Alarcon
#' @param dat a data.frame. Column as variable and rows as sample
#' @param metabolite_list a data.frame. with rows as metabolites.. column names are 'HMDB' and 'Metab'
#' @param save_csv boolean. save all into result into file
#' @param correlation_type character string. type of correlation to use for CRS score 'preason'(default) or 'spearman'



NMRMetabCRS <- function(dat, metabolite_list, correlation_type = 'pearson', save_csv = F) {
  HMDB <- metabolite_list$HMDB # a vector of HMDB from the metabolite_list data
  matches <- metabolite_list$Metab # a second vector of the matching names of HMDBs

  #### DEFINING THE VARIABLE LOCALLY TO THE FUNCTION #####
  CRS_score <- Metabolite <- median_score <- Standard_dev <- sigLevel <- NULL


  corr_data <- list()
  # using grepl to loop over the dataframe
  for (i in 1:nrow(metabolite_list)) {
    single_metab_df <- as.data.frame(dat[, grepl(pattern = HMDB[[i]], x = colnames(dat))])
    if(ncol(single_metab_df) == 1){
      nam = colnames(dat)[grepl(pattern = HMDB[[i]], x = colnames(dat))]
      colnames(single_metab_df) = nam
    }
    corr_data[[i]] <- single_metab_df
  }
  # give the names of the dataframes created in the list to the names of the metabolites for easier reading
  names(corr_data) <- matches

  corr_matrix <- lapply(corr_data, function(x) {
    if (ncol(x) > 0) {
      # correlation matrix
      corr_matrix <- Hmisc::rcorr(as.matrix(x), type = "pearson")$r
      corr_matrix <- as.data.frame(corr_matrix)
    } else {
      # print to the console in case a metabolite in our list is not actually present in our data
      print("0 columns, skipping")
      NULL
    }
  })

  names(corr_matrix) <- matches
  # remove unused (NULL) metabolites from the list
  corr_matrix <- purrr::compact(corr_matrix)

  CRS_scores <- lapply(corr_matrix, function(x) {
    if (ncol(x) > 0) {
      # CRS_score percentace is calculated as the mean of the rows *100
      x %>% dplyr::transmute(metab_bins = rownames(x), CRS_score = rowMeans(x) * 100) %>% tibble::tibble()
    } else {
      # print to the console in case a metabolite in our list is not actually present in our data
      print("0 columns, skipping")
      NULL
    }
  })
  names(CRS_scores) <- names(corr_matrix)

  CRS_pass <- do.call(dplyr::bind_rows, list(CRS_scores, .id = "Metabolite")) %>%
    dplyr::filter(CRS_score != 100, Metabolite != 'UNKNOWN') %>%
    dplyr::summarise(median_score = median(CRS_score), Standard_dev = sd(CRS_score)) %>%
    dplyr::mutate(CRS_pass = round(median_score - Standard_dev)) %>%
    round(digits = 2)

  frequency_metabs <- do.call(dplyr::bind_rows, list(CRS_scores, .id = "Metabolite")) %>%
    count(Metabolite)

  plot_bar <- frequency_metabs %>%
    dplyr::filter(Metabolite != "UNKNOWN") %>%
    ggplot(aes(x = Metabolite, y = n)) +
    geom_bar(stat = "identity") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90)) +
    labs(y = "frequency", title = "Metabolite frequency", subtitle = "number of times each metabilite was identified in binned data")


  summary_data <- dplyr::bind_rows(CRS_scores, .id = 'Metabolite') %>%
    #dplyr::mutate(bin_IDs = gsub("^.*_", "", metab_bins)) %>%
    tidyr::nest(data = c(metab_bins,CRS_score)) %>%
    dplyr::right_join(frequency_metabs, by =c('Metabolite')) %>%
    dplyr::rename(frequency_bin = n) %>%
    dplyr::mutate(
      highest_corr = purrr::map(data, function(x){
        max_score <- x %>% dplyr::pull('CRS_score') %>% max()
      }),
      smallest_corr = purrr::map(data, function(x){
        max_score <- x %>% dplyr::pull('CRS_score') %>% min()
      }),
      bin_IDs = purrr::map(data,function(x){
        x$bin <- gsub("^.*_", "", x$metab_bins)
        paste0(x$bin, collapse = ', ')
      }),
      best_bin = purrr::map(data, function(x){
        best_bin <- x %>% dplyr::slice_max(CRS_score) %>% dplyr::pull(metab_bins)

      })) %>%
    tidyr::unnest(cols = c(highest_corr, smallest_corr, bin_IDs, best_bin))




  output = list('corr_data'= corr_data,
                'corrmatrix'= corr_matrix,
                'CRS_score'= CRS_scores,
                'CRS_pass' = CRS_pass,
                'frequency_table' = frequency_metabs,
                'frequency_plot' = plot_bar,
                'summary_data' = summary_data)

  ##### save crs dataframes into csv

  #
  if (save_csv == T) {
    dplyr::bind_rows(CRS_scores, .id = 'Metabolite') %>% write.csv(file = 'CRS_scores.csv',row.names = F)
  }
  #if (save_excell == T) {

    #wb = xlsx::createWorkbook()

    #lapply(names(CRS_scores),
     #      function(df) {
      #       sheet = xlsx::createSheet(wb, df)
       #      xlsx::addDataFrame(CRS_scores[[df]], sheet = sheet, row.names = FALSE)
        #   } )

    #xlsx::saveWorkbook(wb, "CRS_score.xlsx")
  #}


  return(output)
}




#' @title swap metabolite names to HMDB
#' @name NMRMetab_names_to_HMDB
#' @export
#' @description switches the names to HMDB from the original NMR dataframe
#' @author Michele Fresneda Alarcon
#' @param dat a data.frame. Column as variable and rows as sample
#' @param metabolite_list a data.frame. with rows as metabolites.. column names are 'HMDB' and 'Metab' (each eanty in metab must be equal to the acual names assigend in the dat vatiable)
#' @param index_col column number with first metabolite


NMRMetab_names_to_HMDB <- function(dat, metabolite_list, index_col = 3){

  patt <- metab_list$Metab
  repl <- metab_list$HMDB
  group_columns <- colnames(dat)[1:index_col-1]
  new_names <- lapply(str_split(colnames(dat[index_col:ncol(dat)]), pattern = "_"), function(x) {
    number <- x[length(x)]
    find_names <- grepl(paste0("^", x, "$", collapse = "|"), patt)
    HMDBs <- repl[find_names]
    HMDBs <- paste(HMDBs, collapse = "_")
    HMDBs <- paste(HMDBs, number, sep = "_")
  })


  new_names <- unlist(new_names)

  colnames(dat) <- c(group_columns, new_names)
  return(dat)
}
