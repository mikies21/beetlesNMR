
#' @title foldchange function
#' @name NMRMetab_foldchange
#' @export
#' @description This function will take in either a object of class data.frame and a using a second dataframe with a list of all the metabolite identifed, will return the CRS score for all the bins. this could work with other metabolite identifiers; for example bins with the name of the     metabolite, as  log as they are unique. mistakes could happen when greple take NAD and NADH
#' @author Michele Fresneda Alarcon
#' @param data a data.frame. Column as variable and rows as sample
#' @param groupID string. name off the variable containing your grouping details must only have 2 facors eg Healthy vs Diseases
#' @param dividendID string. top of the division
#' @param index_col index colum oif first mobservation
#' @param divisorID string. what group you're dividing by


NMRMetab_foldchange = function(data, groupID, index_col = 3, dividendID, divisorID){

  coln = colnames(data)[index_col:ncol(data)]

  fold_df = data %>%
    dplyr::select(.data[[groupID]], c(index_col:ncol(data))) %>%
    dplyr::group_by(dplyr::across(.data[[groupID]])) %>%
    dplyr::summarise(dplyr::across(all_of(coln),mean)) %>%
    t() %>%
    as.data.frame() %>%
    janitor::row_to_names(1) %>%
    dplyr::mutate(across(everything(), as.numeric)) %>%
    tibble::rownames_to_column('metabolite') %>%
    tibble::tibble() %>%
    dplyr::mutate(FC = .data[[dividendID]]/.data[[divisorID]])

  return(fold_df)
}

