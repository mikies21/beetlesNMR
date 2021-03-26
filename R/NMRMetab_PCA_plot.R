
#' @title PCA analysis. plot
#' @name NMRMetab_PCA_plot
#' @export
#' @description This function will take in either a object of class data.frame and a using a second dataframe with a list of all the metabolite identifed, will return the CRS score for all the bins. this could work with other metabolite identifiers; for example bins with the name of the     metabolite, as  log as they are unique. mistakes could happen when greple take NAD and NADH
#' @author Michele Fresneda Alarcon
#' @param data a data.frame. Column as variable and rows as sample
#' @param groupID string. name off the fvariable containing your grouping details
#' @param elipses number of components
#' @param index_col index colum oif first mobservation
#' @param pcs vector with integers only. principal components to plot default 1, 2
#' @param size_point double. control the size of the points in the plot

NMRMetab_PCA_plot = function(data, groupID, index_col = 2, elipses = F, pcs = c(1,2), size_point = 3) {

  PCx <- PCy <- PCz <- NULL

  PCA <- prcomp(data[index_col:ncol(data)], center = F, scale. = F)

  drugs_scores <- as.data.frame(PCA$x[, pcs]) ## just taking the first 20 comp. no need to go higher
  prop_var <- round(summary(PCA)$importance[2, ] * 100, digits = 2)
  prop_var <- prop_var[pcs]
  PCnames <- colnames(drugs_scores)[pcs]
  PCnamex <- paste(PCnames[1], " (", prop_var[1], "%)", sep = "")
  PCnamey <- paste(PCnames[2], " (", prop_var[2], "%)", sep = "")
  col_group <- data[, groupID]

  if (length(pcs)==3) {
    PCnamez <- paste(PCnames[3], " (", prop_var[3], "%)", sep = "")
    colnames(drugs_scores) <- c('PCx','PCy','PCz')

    drugs_scores$col_group3d <- col_group
    mycolors <- RColorBrewer::brewer.pal(n = 8,name = 'Dark2')[1:length(unique(col_group))]
    drugs_scores$color <- mycolors[as.numeric(as.factor(drugs_scores$col_group3d))]

    plot1 <- rgl::plot3d(
      x=drugs_scores$PCx, y=drugs_scores$PCy, z=drugs_scores$PCz,
      col = drugs_scores$color,
      type = 'p',
      size = size_point,
      radius = size_point,
      xlab= PCnamex,
      ylab= PCnamey,
      zlab= PCnamez,
      box = F)

  } else{
    colnames(drugs_scores) <- c('PCx','PCy')
    #drugs_scores = cbind.data.frame(group = )
    plot1 <- ggplot2::ggplot(drugs_scores, aes(x = PCx,y = PCy,col = col_group)) +
      ggplot2::geom_point(size = size_point)+
      ggplot2::theme_bw(base_size = 16) +
      ggplot2::labs(col = groupID, x = PCnamex,y = PCnamey)

    if (elipses == T) {
      plot1 = plot1 + ggplot2::stat_ellipse(aes(x = PCx, y = PCy, colour = col_group))
    }
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
