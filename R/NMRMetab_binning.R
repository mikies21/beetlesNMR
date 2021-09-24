#' @title Binning of NMR spectra
#' @name NMRMetab_binning
#' @export
#' @description using a pattern file and the data extracted from bruker experiment files. this function sues the bins min and max to get a value, area under the curve to attribute to tat particular bin
#' @author Michele Fresneda Alarcon
#' @param data a data.frame. Column as samples and rows as ppm. fits columm 'ppm'. this come from bruker topspin
#' @param pattern a data.frame with column as 'min_ppm', max_ppm' and 'bin'. this is the pattern file the you created

#NMRMetab_binning <- function(data, pattern){

#  ppm <- NULL

#  list_of_bins = list()
#  for (i in 1:nrow(pattern_file)) {
#    df = data %>%
#      dplyr::filter(dplyr::between(ppm,pattern_file[i, 'min_ppm'], pattern_file[i, "max_ppm"]))
#    if (nrow(df) == 0) {
#      print(paste0('0 rows ---- ',pattern_file[i,'bin'],'  min=', pattern_file[i, 'min_ppm'],'  max=',pattern_file[i, 'max_ppm']))
#    }
#    list_of_bins[[i]] = df
#    names(list_of_bins)[[i]] = pattern_file[i,'bin']
#  }

#  df = list_of_bins %>% lapply(function(x){
#    ppm = x$ppm
#    #print(ppm)
#    apply(x %>% dplyr::select(-ppm), 2, function(y) {
##      areaucurve = sum(y)
#    })
#  }) %>%
#    dplyr::bind_rows(.id = 'bin') %>%
#    #dplyr::select(-ppm) %>%
#    t() %>%
#    as.data.frame() %>%
#    janitor::row_to_names(1) %>%
#    tibble::rownames_to_column('sampleID') %>%
#    tibble::tibble()
#
#  df[2:ncol(df)] = apply(df[2:ncol(df)], 2 , as.numeric)

#  return(df)
#}



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


# from tame NMR -----------------------------------------------------------

NMRMetab_binning <- function(data, pattern){
  data <- data[2:ncol(data)]
  ppms <- data[,1]

  ppmInterval2Pos <- function(ppms, interval){
    bin = which(ppms>=min(interval) & ppms<=max(interval))
    c(min(bin),max(bin))
  }

  fixDupes = function(labels, dupes){
    unis = unique(labels[dupes])
    for(lab in unis){
      pos = which(labels == lab)
      labels[pos] <- paste(lab, '_', 1:length(pos), sep='')
    }
    labels
  }

  # find and modify duplicated bin labels
  dupes = duplicated(pattern[,'bin'])
  if (any(dupes))
    pattern[,'bin'] = fixDupes(pattern[,'bin'], which(dupes))

  data = data * 1.0
  #convert ppms to positions in the data matrix
  bins = do.call('rbind', lapply(1:nrow(pattern), function(i) ppmInterval2Pos(ppms, pattern[i,c('max_ppm','min_ppm')])))
  binSize = abs(bins[,2] - bins[,1]) + 1
  print(bins)
  print(dim(bins))
  dataInt = do.call('cbind', lapply(1:nrow(bins), function(i) apply(data[bins[i,1]:bins[i,2],], 2, sum)/binSize[i]))
  dataInt = as.data.frame(dataInt, stringsAsFactors=F)

  dupes = duplicated(colnames(data))
  if (any(dupes)){
    tempNames = fixDupes(colnames(data), which(dupes))
  } else {
    tempNames = colnames(data)
  }

  rownames(dataInt) = tempNames
  names(dataInt) = pattern[,'bin']
  dataInt
}

