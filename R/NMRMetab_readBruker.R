#' @title read Bruker files
#' @name NMRMmetab_readBruker
#' @export
#' @description takes the path where al the experiment files are, returns a dataframe with ppm as first column and samples as subsequent column. not working properly yet as evn tho data was previously aligned in topspin this was not reflected in this.
#' @author Dr. Jie Hao
#' @author Michele Fresneda Alarcon
#' @param BrukerDataDir a string. path to the NMR experiment folder


NMRMmetab_readBruker<- function(path_to_file){
  #path <- paste(system.file(package = 'beetlesNMR'), 'read_bruker_pyt.py', sep = '/')
  #command <- paste('python', path, path_to_file)
  #try(suppressWarnings(response <- system(command, intern = T)), silent = T)

  #if (!is.null(attr(response, 'status'))){
  #  if(attr(response, 'status') == 1){
  #    response <- ""
  #    cat('something went wrong. not sure what')
  #  }
  #}

  #reticulate::source_python(file = 'C:/Users/micfres/Documents/R/win-library/4.1/beetlesNMR/.inst/read_bruker_pyt.py')
  #read_read_data(path_to_file)
}

