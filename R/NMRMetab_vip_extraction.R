
#' @title VIP extraction
#' @name NMRMetab_PLSDA_VIP
#' @export
#' @description This function will take in either a object of class data.frame and a using a second dataframe with a list of all the metabolite identifed, will return the CRS score for all the bins. this could work with other metabolite identifiers; for example bins with the name of the     metabolite, as  log as they are unique. mistakes could happen when greple take NAD and NADH
#' @author Michele Fresneda Alarcon
#' @param dat a data.frame. Column as variable and rows as sample
#' @param metabolite_list a data.frame. with rows as metabolites.. column names are 'HMDB' and 'Metab'
#' @param save_excell boolean. save all into result into file


get_vips_dat = function(dat){
  plsda_model <- plsda(
    X = dat[, 9:ncol(dat)],
    Y = factor(dat[, 'disease_status']),
    ncomp = 5
  ) ### change number of components appropiately

  vips = vip(plsda_model) %>%
    as.data.frame() %>%
    rownames_to_column('metabolite') %>%
    tibble %>%
    #dplyr::filter(comp2 > 1) %>%
    arrange(desc(comp1)) %>%
    dplyr::filter(comp1 > 1)

  met =vips$metabolite

  plot_vips = ggplot2::ggplot(vips, aes(x = reorder(metabolite, comp1), y = comp1)) +
    geom_point(size = 3, col = "#B31B21") +
    coord_flip() +
    theme_bw(base_size = 13) +
    labs(x = 'metabolite',title = 'VIP scores of Healthy vs Disease neutrophils',
         subtitle = 'neutrophils treated with JAK inhibitor')

  output = list('VIP_df' = vips, 'plot_vips' = plot_vips)
} ##gettting component 1
