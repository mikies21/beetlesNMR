#' @title read Bruker files
#' @name NMRMmetab_readBruker
#' @export
#' @description takes the path where al the experiment files are, returns a dataframe with ppm as first column and samples as subsequent column. not working properly yet as evn tho data was previously aligned in topspin this was not reflected in this.
#' @author Dr. Jie Hao
#' @author Michele Fresneda Alarcon
#' @param path_to_file a string. path to the NMR experiment folder


NMRMmetab_readBruker<- function(path_to_file){
  path <- paste(system.file(package = 'beetlesNMR'), 'python/read_bruker_pyt.py', sep = '/')

  reticulate::source_python(file = path)
  dat <- read_data(path_to_file)
  dat <- tibble::rownames_to_column(.data = dat, var = 'ppm')
  return(dat)
}

#NMRMmetab_readBruker(path_to_file = '../../../Documents/8. RAW_NMR/all_batches_QCd/')
