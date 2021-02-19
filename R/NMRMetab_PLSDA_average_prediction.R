

#' @title NMRMetab_average_prediction
#' @name NMRMetab_PLSDA_average_prediction
#' @export
#' @description This function will take in either a object of class data.frame and a using a second dataframe with a list of all the metabolite identifed, will return the CRS score for all the bins. this could work with other metabolite identifiers; for example bins with the name of the     metabolite, as  log as they are unique. mistakes could happen when greple take NAD and NADH
#' @author Michele Fresneda Alarcon
#' @param dat a data.frame. Column as variable and rows as sample
#' @param groupID a data.frame. with rows as metabolites.. column names are 'HMDB' and 'Metab'
#' @param components_model boolean. save all into result into file
#' @param iterations something



average_prediction_metrics = function(dat,
                                      groupID = 'disease_status',
                                      index_col = 9,
                                      components_model = 5,
                                      iterations = 5){
  predictions_list <- list()
  pred = list()
  vips = list()

  for (i in 1:components_model) {
    set.seed(21)
    for(j in 1:iterations){
      splits = rsample::initial_split(data = dat, strata = 'disease_status')

      dat_train = rsample::training(splits) %>%
        NMRMetab_norm_scale(index_col = index_col, scaling = 'Pareto') %>%
        invisible()
      dat_test = rsample::testing(splits) %>%
        NMRMetab_norm_scale(index_col = index_col, scaling = 'Pareto') %>%
        invisible()



      plsda_model <- mixOmics::plsda(
        X = dat_train[, index_col:ncol(dat_train)],
        Y = factor(dat_train[, groupID]),
        ncomp = components_model
      ) ### change number of components appropiately

      vips[[j]] = mixOmics::vip(plsda_model) %>%
        as.data.frame() %>% tibble::rownames_to_column('metabolite')

      prediction <- tibble::tibble(
        "prediction" = factor(predict(
          plsda_model,
          dat_test[index_col:ncol(dat_test)]
        )$class$mahalanobis.dist[, i]), ### change number of components appropiately
        "true_val" = factor(dat_test[, groupID])
      )
      #     prediction$prediction = relevel(prediction$prediction, ref = 'Healthy')#
      #     prediction$true_val = relevel(prediction$true_val, ref = 'Healthy')
      pred[[j]] <- prediction %>%
        yardstick::conf_mat(truth = "true_val", estimate = "prediction") %>%
        summary() %>%
        dplyr::filter(.metric %in% c("accuracy", "bal_accuracy", "precision", "recall", "f_meas")) %>%
        dplyr::select(-.estimator)

    }

    predictions_list[[i]] = pred %>%
      dplyr::bind_rows(.id = "run") %>%
      dplyr::group_by(.metric) %>%
      dplyr::summarise(estimate = round(mean(.estimate) * 100, digits = 1), sd = round(sd(.estimate) * 100, digits = 1))

  }

  predictions_list = predictions_list %>% dplyr::bind_rows(.id ='component')
  print(predictions_list)

  names(vips) = c(paste(rep('run',iterations), 1:iterations, sep = ''))
  vips = vips %>%
    dplyr::bind_rows(.id ='run')

  vips = vips %>% dplyr::select(run, metabolite, comp2) %>%
    tidyr::pivot_wider(names_from = 'metabolite',values_from = 'comp2') %>%
    t() %>% janitor::row_to_names(1) %>%as.data.frame() %>% tibble::rownames_to_column('metabolite') %>% tibble::as_tibble()
  vips[1:iterations+1] <- sapply(vips[1:iterations+1], as.numeric)
  return(list('predictions' = predictions_list,'vips' = vips))
}
