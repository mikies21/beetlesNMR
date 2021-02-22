<<<<<<< HEAD
#' @title compute PLS-DA prediction metrics, VIPs
#' @name NMRmetab_average_prediction_metrics
#' @description alot of stuff
#' @param dat a matrix or a data.frame. rows are samples and columns are grouping and metabolites or bins
#' @param groupID a vercors with the column name contining the grouping of the data. This parameter is to be use to mae sure that the 2 new data.frames created have all groups mentioned in the column we identifiy
#' @param index_col name of the column containing the groups
#' @param components_model number of components in the model
#' @param iterations number of iterations. ie number of times you are running the random splits -> plsda -> test prediction
#' @param run_CV boolean. default is false. intensive computer work. performs the CV on each train split for the maximum number of comomponent.
#' @param CV_validation string. either 'loo' default of 'Mfold'
#' @param CV_folds integer. number of folds. only used for Mfold validation
#' @param CV_repeats integer. number of repeats only used for Mfold validation
#' @export


NMRmetab_average_prediction_metrics = function(dat,
                                      groupID = 'disease_status',
                                      index_col = 9,
                                      components_model = 5,
                                      iterations = 5,
                                      run_CV = F,
                                      CV_validation = 'loo',
                                      CV_folds = 5,
                                      CV_repeats = 100){
  CV <- repeats <- .metric <- .estimator <- .estimate <- run <- comp <- type <- mahalanobis.dist <- metabolite <- NULL
  predictions_list <- list()
  pred = list()
  vips = list()
  cvs = list()


  for (i in 1:components_model) {
    set.seed(21)
    for(j in 1:iterations){
      splits = rsample::initial_split(data = dat, strata = 'disease_status')

      dat_train = rsample::training(splits) %>%
        NMRMetab_norm_scale(index_col = index_col, scaling = 'Pareto')

      dat_test = rsample::testing(splits) %>%
        NMRMetab_norm_scale(index_col = index_col, scaling = 'Pareto')



      plsda_model <- mixOmics::plsda(
        X = dat_train[, index_col:ncol(dat_train)],
        Y = factor(dat_train[, groupID]),
        ncomp = components_model
      ) ### change number of components appropiately

      if (run_CV & components_model == i) {
        cv_df = mixOmics::perf(object = plsda_model, validation = CV_validation, folds = CV_folds, nrepeats = CV-repeats, progressBar = T)
        cv_df = cv_df[['error.rate']] %>%
          lapply(data.frame) %>%
          dplyr::bind_rows(.id = 'type')
        cv_df$comp = c(1:(nrow(cv_df)/2))
          #dplyr::mutate(comp = factor(c(1:(nrow(error.dfs)/2))))
        #cv_df$comp = factor(c(1:(nrow(error.dfs)/2)))
        #error.dfs = reshape2::melt(error.dfs)
        print(cv_df)
        cvs[[j]] = cv_df
      }

      vips[[j]] = mixOmics::vip(plsda_model) %>%
        as.data.frame() %>%
        tibble::rownames_to_column('metabolite')

      prediction <- tibble::tibble(
        "prediction" = factor(stats::predict(
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
  median_vips = vips %>%
    dplyr::select(-run) %>%
    dplyr::group_by(metabolite) %>%
    dplyr::summarise(dplyr::across(.cols = dplyr::starts_with('comp'),median))
  mean_vips = vips %>%
    dplyr::select(-run) %>%
    dplyr::group_by(metabolite) %>%
    dplyr::summarise(dplyr::across(.cols = dplyr::starts_with('comp'),mean))

  # GET MEAN AND MEDIAN FOR VIPS FOR EACH COMPONENT

  #vips = vips %>%
  #  dplyr::select(run, metabolite, comp2) %>%
  #  tidyr::pivot_wider(names_from = 'metabolite',values_from = 'comp2') %>%
  #  t() %>%
  #  janitor::row_to_names(1) %>%
  #  as.data.frame() %>%
  #  tibble::rownames_to_column('metabolite') %>%
  #  tibble::as_tibble()

#  vips[1:iterations+1] <- sapply(vips[1:iterations+1], as.numeric)
  names(cvs) = c(paste(rep('run',iterations), 1:iterations, sep = ''))
  cvs = dplyr::bind_rows(cvs, .id = 'run')
  p = cvs %>%
    dplyr::select(run, comp, type, mahalanobis.dist) %>%
    dplyr::filter(type == 'BER') %>%
    ggplot2::ggplot(ggplot2::aes(x = comp,
                                 y = mahalanobis.dist,
                                 linetype = type,
                                 colour = run,
                                 )) +
    ggplot2::geom_line(stat = 'identity')+
    ggplot2::geom_point()+
    ggplot2::theme_bw()

  return(list('predictions' = predictions_list,
              'vips' = vips,
              'CVs' = cvs,
              'cv_plot' = p,
              'median_vips' = median_vips,
              'mean_vips' = mean_vips))
}

=======
#' @title compute PLS-DA prediction metrics, VIPs
#' @name NMRmetab_average_prediction_metrics
#' @description alot of stuff
#' @param dat a matrix or a data.frame. rows are samples and columns are grouping and metabolites or bins
#' @param groupID a vercors with the column name contining the grouping of the data. This parameter is to be use to mae sure that the 2 new data.frames created have all groups mentioned in the column we identifiy
#' @param index_col name of the column containing the groups
#' @param components_model number of components in the model
#' @param iterations number of iterations. ie number of times you are running the random splits -> plsda -> test prediction
#' @param run_CV boolean. default is false. intensive computer work. performs the CV on each train split for the maximum number of comomponent.
#' @param CV_validation string. either 'loo' default of 'Mfold'
#' @param CV_folds integer. number of folds. only used for Mfold validation
#' @param CV_repeats integer. number of repeats only used for Mfold validation
#' @export


average_prediction_metrics = function(dat,
                                      groupID = 'disease_status',
                                      index_col = 9,
                                      components_model = 5,
                                      iterations = 5,
                                      run_CV = F,
                                      CV_validation = 'loo',
                                      CV_folds = 5,
                                      CV_repeats = 100){
  CV <- repeats <- .metric <- .estimator <- .estimate <- run <- comp <- type <- mahalanobis.dist <- metabolite <- NULL
  predictions_list <- list()
  pred = list()
  vips = list()
  cvs = list()


  for (i in 1:components_model) {
    set.seed(21)
    for(j in 1:iterations){
      splits = rsample::initial_split(data = dat, strata = 'disease_status')

      dat_train = rsample::training(splits) %>%
        NMRMetab_norm_scale(index_col = index_col, scaling = 'Pareto')

      dat_test = rsample::testing(splits) %>%
        NMRMetab_norm_scale(index_col = index_col, scaling = 'Pareto')



      plsda_model <- mixOmics::plsda(
        X = dat_train[, index_col:ncol(dat_train)],
        Y = factor(dat_train[, groupID]),
        ncomp = components_model
      ) ### change number of components appropiately

      if (run_CV & components_model == i) {
        cv_df = mixOmics::perf(object = plsda_model, validation = CV_validation, folds = CV_folds, nrepeats = CV-repeats, progressBar = T)
        cv_df = cv_df[['error.rate']] %>%
          lapply(data.frame) %>%
          dplyr::bind_rows(.id = 'type')
        cv_df$comp = c(1:(nrow(cv_df)/2))
          #dplyr::mutate(comp = factor(c(1:(nrow(error.dfs)/2))))
        #cv_df$comp = factor(c(1:(nrow(error.dfs)/2)))
        #error.dfs = reshape2::melt(error.dfs)
        print(cv_df)
        cvs[[j]] = cv_df
      }

      vips[[j]] = mixOmics::vip(plsda_model) %>%
        as.data.frame() %>%
        tibble::rownames_to_column('metabolite')

      prediction <- tibble::tibble(
        "prediction" = factor(stats::predict(
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
  median_vips = vips %>%
    dplyr::select(-run) %>%
    dplyr::group_by(metabolite) %>%
    dplyr::summarise(dplyr::across(.cols = dplyr::starts_with('comp'),median))
  mean_vips = vips %>%
    dplyr::select(-run) %>%
    dplyr::group_by(metabolite) %>%
    dplyr::summarise(dplyr::across(.cols = dplyr::starts_with('comp'),mean))

  # GET MEAN AND MEDIAN FOR VIPS FOR EACH COMPONENT

  #vips = vips %>%
  #  dplyr::select(run, metabolite, comp2) %>%
  #  tidyr::pivot_wider(names_from = 'metabolite',values_from = 'comp2') %>%
  #  t() %>%
  #  janitor::row_to_names(1) %>%
  #  as.data.frame() %>%
  #  tibble::rownames_to_column('metabolite') %>%
  #  tibble::as_tibble()

#  vips[1:iterations+1] <- sapply(vips[1:iterations+1], as.numeric)
  names(cvs) = c(paste(rep('run',iterations), 1:iterations, sep = ''))
  cvs = dplyr::bind_rows(cvs, .id = 'run')
  p = cvs %>%
    dplyr::select(run, comp, type, mahalanobis.dist) %>%
    dplyr::filter(type == 'BER') %>%
    ggplot2::ggplot(ggplot2::aes(x = comp,
                                 y = mahalanobis.dist,
                                 linetype = type,
                                 colour = run,
                                 )) +
    ggplot2::geom_line(stat = 'identity')+
    ggplot2::geom_point()+
    ggplot2::theme_bw()

  return(list('predictions' = predictions_list,
              'vips' = vips,
              'CVs' = cvs,
              'cv_plot' = p,
              'median_vips' = median_vips,
              'mean_vips' = mean_vips))
}

>>>>>>> b1f0e94cf480a7b6e5e9136fbc6f6428acf4b709
