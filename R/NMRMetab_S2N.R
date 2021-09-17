
#' @title Signal to Noise calculation
#' @name NMRmetab_S2N
#' @export
#' @description this function will calculate the signal to noise ration for each sample and each bin.
#' @author Michele Fresneda Alarcon
#' @param data a data.frame. Column as samples and rows as ppm. fits columm 'ppm'. this come from bruker topspin
#' @param pattern_file a data.frame with column as 'min_ppm', max_ppm' and 'bin'. this is the pattern file the you created


# signal to noise for each individual spectra

NMRmetab_S2N <-  function(data, index_column, pattern_file, noise_column){

  pattern <- pattern_file %>% mutate(width = max_ppm-min_ppm)
  data_new <- data[, -noise_column]
  S2N_calc <- data %>%
    tidyr::pivot_longer(cols = index_column:ncol(data), names_to = 'bin', values_to = 'value') %>%
    dplyr::left_join(pattern, by =  'bin') %>%
    dplyr::mutate(intensity_over_width = value/width, S2N = intensity_over_width/NOISE) %>%
    tidyr::pivot_wider(id_cols = 1:(index_column-1), names_from = 'bin', values_from = 'S2N')
  return(S2N_calc)

}
