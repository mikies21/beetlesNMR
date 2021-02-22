
#' @title PCA analysis. plot
#' @name NMRMetab_PCA_plot
#' @export
#' @description This function will take in either a object of class data.frame and a using a second dataframe with a list of all the metabolite identifed, will return the CRS score for all the bins. this could work with other metabolite identifiers; for example bins with the name of the     metabolite, as  log as they are unique. mistakes could happen when greple take NAD and NADH
#' @author Michele Fresneda Alarcon
#' @param dat a data.frame. Column as variable and rows as sample
#' @param groupID string. name off the fvariable containing your grouping details
#' @param elipses number of components
#' @param index_col index colum oif first mobservation

NMRMetab_PCA_plot = function(data, groupID, index_col = 2, elipses = F) {
  PCA <- prcomp(data[index_col:ncol(data)], center = F, scale. = F)


  drugs_scores <- as.data.frame(PCA$x[, c(1:3)]) ## just taking the first 20 comp. no need to go higher
  prop_var <- round(summary(PCA)$importance[2, ] * 100, digits = 2)
  prop_var <- prop_var[c(1,2)]
  PCnames <- colnames(drugs_scores)[c(1, 2)]
  PCnamex <- paste(PCnames[1], " (", prop_var[1], "%)", sep = "")
  PCnamey <- paste(PCnames[2], " (", prop_var[2], "%)", sep = "")
  col_group <- data[, groupID]
  #drugs_scores = cbind.data.frame(group = )
  plot1 <- ggplot(
    drugs_scores, aes(
      x = PC1,
      y = PC2,
      col = col_group
    )
  ) +
    geom_point(size = 3)+
    theme_bw(base_size = 16) +
    labs(
      col = grp,
      x = PCnamex,
      y = PCnamey
    )+
    scale_color_brewer(palette = 'Dark2')

  if (elipses == T) {
    plot1 = plot1 + stat_ellipse(aes(x = PC1, y = PC2, colour = col_group))
  }
  return(plot1)

}


loadings_plot = function(x, index_column = 2) {

  PCA <- prcomp(x[index_column:ncol(x)], center = F, scale. = F)
  loadings <- as.data.frame(PCA$rotation[, 1:3])
  plot2 <- ggplot(loadings, aes(x = PC1, y = PC2)) +
    geom_point() +
    theme_bw(base_size = 16) +
    labs(title = "PCA loadings")
  return(plot2)
}
