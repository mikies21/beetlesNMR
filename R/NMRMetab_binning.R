<<<<<<< HEAD
#' @title Binning of NMR spectra
#' @name NMRMetab_binning
#' @export
#' @author Michele Fresneda Alarcon
#' @param data a data.frame. Column as samples and rows as ppm. fits columm 'ppm'. this come from bruker topspin
#' @param pattern_file a data.frame with column as 'min_ppm', max_ppm' and 'bin'. this is the pattern file the you created

NMRMetab_binning <- function(data, pattern_file){

  ppm <- NULL

  list_of_bins = list()
  for (i in 1:nrow(pattern_file)) {
    df = data %>%
      dplyr::filter(dplyr::between(ppm,pattern_file[i, 'min_ppm'], pattern_file[i, "max_ppm"]))
    if (nrow(df) == 0) {
      print(paste0('0 rows ---- ',pattern_file[i,'bin'],'  min=', pattern_file[i, 'min_ppm'],'  max=',pattern_file[i, 'max_ppm']))
    }
    list_of_bins[[i]] = df
    names(list_of_bins)[[i]] = pattern_file[i,'bin']
  }

  df = list_of_bins %>% lapply(function(x){
    ppm = x$ppm
    #print(ppm)
    apply(x %>% dplyr::select(-ppm), 2, function(y) {
      areaucurve = bayestestR::area_under_curve(x = ppm, y = y, method = 'trapezoid') * 1000
    })
  }) %>%
    dplyr::bind_rows(.id = 'bin') %>%
    #dplyr::select(-ppm) %>%
    t() %>%
    as.data.frame() %>%
    janitor::row_to_names(1) %>%
    tibble::rownames_to_column('sampleID') %>%
    tibble::tibble()

  df[2:ncol(df)] = apply(df[2:ncol(df)], 2 , as.numeric)

  return(df)
}



#list_of_bins = apply(pattern_file, 1, function(row){
#  patt_file_ppm = suppressWarnings(row %>% as.numeric() %>% na.omit() %>% sort())
#  min_ppm = patt_file_ppm[1]
#  max_ppm = patt_file_ppm[2]
#  df = data %>%
#    dplyr::filter(dplyr::between(ppm,min_ppm,max_ppm))
#  if (nrow(df) == 0) {
#    print(paste0('0 rows ---- ',pattern_file[i,'bin'],'  ',pattern_file[i, 'min_ppm'],pattern_file[i, 'max_ppm']))
#  }
#})
#names(list_of_bins) = pattern_file$bin
=======
#' @title Binning of NMR spectra
#' @name NMRMetab_binning
#' @export
#' @author Michele Fresneda Alarcon
#' @param data a data.frame. Column as samples and rows as ppm. fits columm 'ppm'. this come from bruker topspin
#' @param pattern_file a data.frame with column as 'min_ppm', max_ppm' and 'bin'. this is the pattern file the you created

NMRMetab_binning <- function(data, pattern_file){

  ppm <- NULL

  list_of_bins = list()
  for (i in 1:nrow(pattern_file)) {
    df = data %>%
      dplyr::filter(dplyr::between(ppm,pattern_file[i, 'min_ppm'], pattern_file[i, "max_ppm"]))
    if (nrow(df) == 0) {
      print(paste0('0 rows ---- ',pattern_file[i,'bin'],'  min=', pattern_file[i, 'min_ppm'],'  max=',pattern_file[i, 'max_ppm']))
    }
    list_of_bins[[i]] = df
    names(list_of_bins)[[i]] = pattern_file[i,'bin']
  }

  df = list_of_bins %>% lapply(function(x){
    ppm = x$ppm
    #print(ppm)
    apply(x %>% dplyr::select(-ppm), 2, function(y) {
      areaucurve = bayestestR::area_under_curve(x = ppm, y = y, method = 'trapezoid') * 1000
    })
  }) %>%
    dplyr::bind_rows(.id = 'bin') %>%
    #dplyr::select(-ppm) %>%
    t() %>%
    as.data.frame() %>%
    janitor::row_to_names(1) %>%
    tibble::rownames_to_column('sampleID') %>%
    tibble::tibble()

  df[2:ncol(df)] = apply(df[2:ncol(df)], 2 , as.numeric)

  return(df)
}



#list_of_bins = apply(pattern_file, 1, function(row){
#  patt_file_ppm = suppressWarnings(row %>% as.numeric() %>% na.omit() %>% sort())
#  min_ppm = patt_file_ppm[1]
#  max_ppm = patt_file_ppm[2]
#  df = data %>%
#    dplyr::filter(dplyr::between(ppm,min_ppm,max_ppm))
#  if (nrow(df) == 0) {
#    print(paste0('0 rows ---- ',pattern_file[i,'bin'],'  ',pattern_file[i, 'min_ppm'],pattern_file[i, 'max_ppm']))
#  }
#})
#names(list_of_bins) = pattern_file$bin
>>>>>>> b1f0e94cf480a7b6e5e9136fbc6f6428acf4b709
