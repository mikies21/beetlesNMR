
#' @title VIP extraction
#' @name NMRMetab_PLSDA_VIP
#' @export
#' @description This function will take in either a object of class data.frame and a using a second dataframe with a list of all the metabolite identifed, will return the CRS score for all the bins. this could work with other metabolite identifiers; for example bins with the name of the     metabolite, as  log as they are unique. mistakes could happen when greple take NAD and NADH
#' @author Michele Fresneda Alarcon
#' @param dat a data.frame. Column as variable and rows as sample
#' @param groupID string. name off the fvariable containing your grouping details
#' @param ncomp number of components
#' @param index_col index colum oif first mobservation

NMRMetab_get_vips = function(dat, groupID = 1, index_col = 3, ncomp = 2){

  . <- metabolite <- comp <-

  datgrp <- dat %>% dplyr::pull(var = groupID)
  datmatrix <- dat[, index_col:ncol(dat)]
  plsda_model <- mixOmics::plsda(
    X = dat[, index_col:ncol(dat)],
    Y = factor(dat[, groupID]),
    ncomp = ncomp
  ) ### change number of components appropiately

  vips = mixOmics::vip(plsda_model) %>%
    as.data.frame() %>%
    tibble::rownames_to_column('metabolite') %>%
    tibble::tibble() %>%
    dplyr::select(c(1,ncomp+1))%>%
    #dplyr::filter(comp2 > 1) %>%
    dplyr::arrange(desc(.[[2]])) %>%
    dplyr::filter(.[[2]] > 1)

  colnames(vips) = c('metabolite', 'comp')

  met = vips$metabolite

  plot_vips = ggplot2::ggplot(vips, aes(x = reorder(metabolite, comp), y = comp)) +
    ggplot2::geom_point(size = 3, col = "#B31B21") +
    ggplot2::coord_flip() +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::labs(x = 'metabolite',title = paste0('VIP for comp ', ncomp))

  plot(plot_vips)
  output = list('VIP_df' = vips, 'plot_vips' = plot_vips)
} ##gettting component 1
