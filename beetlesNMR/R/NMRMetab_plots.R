

#' @title plot NMR_data
#' @name NMRMetab_plot_binned
#' @export
#' @description plot binned spectra
#' @author Michele Fresneda Alarcon
#' @param binned_data a data.frame. Column as variable and rows as sample
#' @param index_col integer. index number of the colum with the firsst metabllite measurement
#' @param group_var string. the name of the grouping vairable. defails is sampleID. change appropiately



NMRMetab_plot_binned <- function(binned_data, index_col = 2, group_var = 'sampleID') {

  x <- grp <- y <- NULL
    #if (is_null(group_var)) {
    # reshape2::melt(binned_data) %>%
    #   ggplot2::ggplot(aes(x = variable, y = value, col = groupID, group = groupID)) +
    #   ggplot2::geom_line(show.legend = F) +
    #   ggplot2::theme_bw(base_size = 7) +
    #   ggplot2::theme(axis.text.x  = element_text(angle = 90))
    # }

    datgrp <- binned_data %>% dplyr::pull(var = group_var)
    datmatrix <- binned_data[, index_col:ncol(binned_data)]
    dataTemp <- data.frame(
      x = rep(1:ncol(datmatrix), nrow(datmatrix)),
      y = as.vector(t(datmatrix)),
      grp = factor(rep(datgrp, each = ncol(datmatrix)))
    )

    sample_sum <- dataTemp %>%
      dplyr::group_by(x, grp) %>%
      dplyr::summarize(
        mean = mean(y),
        sd = sd(y),
        mean_p2sd = mean + 2 * sd,
        mean_m2sd = mean - 2 * sd
      ) %>%
      dplyr::ungroup()

    p <- ggplot2::ggplot(data = sample_sum, aes(x = x, y = mean, group = grp, col = grp, stat = "identity")) +
      geom_line()

    p <- p +
      theme_bw() +
      ggtitle("NMR bins") +
      xlab("Bin") +
      ylab("Intensity") +
      theme(plot.title = element_text(hjust = 0.5), axis.text.x  = element_text(angle = 90),legend.position = 'none')

    plot(p)
    return(p)


  }


#' @title plot NMR_data with bins
#' @name NMRMetab_plot_raw_with_bins
#' @export
#' @description plot binned spectra
#' @author Michele Fresneda Alarcon
#' @param raw_data a data.frame. Column as variable and rows as sample
#' @param pattern_file pattern file with bins
#' @param min_x double
#' @param max_x double
#' @param ymax int height of the rectangle and of the text for the name of the bin


NMRMetab_plot_raw_with_bins <- function(raw_data, pattern_file = NULL, min_x, max_x, ymax = 1000000) {

  ppm <- min_ppm <- max_ppm <- bin <- value <- sampleID <- change <- NULL



  binnend_raw = raw_data %>%
    dplyr::filter(between(ppm,left = min_x, right = max_x))


  if (is.null(pattern_file)) {
    plot1 = binnend_raw %>%
      tidyr::pivot_longer(!ppm, names_to = 'sampleID',values_to = 'value') %>%
      ggplot2::ggplot() +
      ggplot2::geom_line(aes(x = ppm, y = value, col = sampleID),show.legend = F) +
      ggplot2::theme_bw(base_size = 7)+
      ggplot2::scale_x_reverse()
  }
  else{
    binnend_pattern = pattern_file %>%
      dplyr::filter(min_ppm > min_x & max_ppm < max_x)


    binnend_pattern = binnend_pattern %>% dplyr::mutate(change = max_ppm - min_ppm)

    plot1 = binnend_raw %>%
      tidyr::pivot_longer(!ppm, names_to = 'sampleID',values_to = 'value') %>%
      ggplot2::ggplot() +
      ggplot2::geom_line(aes(x = ppm, y = value, col = sampleID),show.legend = F) +
      ggplot2::theme_bw(base_size = 7) +
      geom_rect(data = binnend_pattern, aes(xmin = min_ppm, xmax =max_ppm, ymin = 0, ymax = ymax), alpha = 0.2, col = 'black')+
      ggrepel::geom_text_repel(data= binnend_pattern, aes(x=min_ppm+change/2, y=ymax, label=bin,angle = 45), size=3)+
      ggplot2::scale_x_reverse()
  }
  plot1
  return(plot1)

}
