
#' @title PCA analysis. plot
#' @name NMRMetab_PCA_plot
#' @export
#' @description This function will take in either a object of class data.frame and a using a second dataframe with a list of all the metabolite identifed, will return the CRS score for all the bins. this could work with other metabolite identifiers; for example bins with the name of the     metabolite, as  log as they are unique. mistakes could happen when greple take NAD and NADH
#' @author Michele Fresneda Alarcon
#' @param data a data.frame. Column as variable and rows as sample
#' @param groupID string. name off the fvariable containing your grouping details
#' @param elipses number of components
#' @param index_col index colum oif first mobservation
#' @param pcs principal component to plot

NMRMetab_PCA_plot = function(data, groupID, index_col = 2, elipses = F, pcs = c(1,2)) {

  PCx <- PCy <- NULL

  PCA <- prcomp(data[index_col:ncol(data)], center = F, scale. = F)

  drugs_scores <- as.data.frame(PCA$x[, pcs]) ## just taking the first 20 comp. no need to go higher
  prop_var <- round(summary(PCA)$importance[2, ] * 100, digits = 2)
  prop_var <- prop_var[pcs]
  PCnames <- colnames(drugs_scores)[pcs]
  PCnamex <- paste(PCnames[1], " (", prop_var[1], "%)", sep = "")
  PCnamey <- paste(PCnames[2], " (", prop_var[2], "%)", sep = "")
  col_group <- data[, groupID]
  colnames(drugs_scores) <- c('PCx','PCy')
  #drugs_scores = cbind.data.frame(group = )
  plot1 <- ggplot2::ggplot(drugs_scores, aes(x = PCx,y = PCy,col = col_group)) +
    ggplot2::geom_point(size = 3)+
    ggplot2::theme_bw(base_size = 16) +
    ggplot2::labs(col = groupID, x = PCnamex,y = PCnamey) +
    ggplot2::scale_color_brewer(palette = 'Dark2')

  if (elipses == T) {
    plot1 = plot1 + ggplot2::stat_ellipse(aes(x = PCx, y = PCy, colour = col_group))
  }
  return(plot1)

}


#' @title PCA analysis loading plot
#' @name NMRMetab_PCA_loading_plot
#' @export
#' @description This function will take in either a object of class data.frame and a using a second dataframe with a list of all the metabolite identifed, will return the CRS score for all the bins. this could work with other metabolite identifiers; for example bins with the name of the     metabolite, as  log as they are unique. mistakes could happen when greple take NAD and NADH
#' @author Michele Fresneda Alarcon
#' @param data a data.frame. Column as variable and rows as sample
#' @param pcs string. name off the fvariable containing your grouping details
#' @param index_col index colum oif first mobservation

NMRMetab_PCA_loading_plot = function(data, index_col = 2, pcs = c(1,2)) {

  PCx <- PCy <- NULL

  PCA <- prcomp(data[index_col:ncol(data)], center = F, scale. = F)

  drugs_scores <- as.data.frame(PCA$x[, pcs]) ## just taking the first 20 comp. no need to go higher
  prop_var <- round(summary(PCA)$importance[2, ] * 100, digits = 2)
  prop_var <- prop_var[pcs]
  PCnames <- colnames(drugs_scores)[pcs]
  PCnamex <- paste(PCnames[1], " (", prop_var[1], "%)", sep = "")
  PCnamey <- paste(PCnames[2], " (", prop_var[2], "%)", sep = "")


  loadings <- as.data.frame(PCA$rotation[, pcs])
  colnames(loadings) <- c('PCx', 'PCy')
  plot2 <- ggplot2::ggplot(loadings, aes(x = PCx, y = PCy, label = rownames(loadings))) +
    ggplot2::geom_point() +
    ggplot2::theme_bw(base_size = 16) +
    ggplot2::labs(title = "PCA loadings", x = PCnamex, y = PCnamey)
  return(plot2)
}
