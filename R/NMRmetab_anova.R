#' @title Normalisation and Scaling
#' @name NMRMetab_anova
#' @export
#' @author Eva Caamano Gutierrez
#' @author Arturas Grauslys
#' @param data a data.frame. Column as variable and rows as sample
#' @param index_col integer value. The column number with first bin/matabolite
#' @param group_test charachter string. name of the group column to test
#' @param sigLevel double. set to 0.05, significance level to detect differences between groups,
#' @param adjMethod charachter string. multiple comparison adjustment method. Available methods: "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none". "BH" stands for Benjamini-Hochberg; "BY" for Benjamini-Yakutieli; "fdr" is false discovery rate (same as "BH")
#' @importFrom utils write.csv
#' @importFrom graphics abline text
#' @importFrom stats TukeyHSD anova aov p.adjust



NMRMetab_anova = function(data, index_col = 3, group_test = 'group', sigLevel = 0.05, adjMethod = 'fdr'){

  # PREP
  # separate data from groups
  data_ = as.matrix(data[,index_col:ncol(data)])
  grp = as.factor(data[,group_test])
  outFolder = getwd()


  cat(sprintf('\nPerforming ANOVA with %s multiple testing correction and %.2f significance level\n', adjMethod, sigLevel))
  cat(sprintf('and Tukey post-hoc analysis with %d pairwise comparisons.\n\n', 2^length(levels(grp))))
  res = do_anova_m(data_, grp, adjustMethod = adjMethod, thresh = sigLevel)
  if (ncol(res$tukey_pvals) != 0){
    print_anova_res(res, sigLvl = sigLevel)
  }



  p = plot.Pvals(res$anova_pvals,sigLvl = sigLevel)
  p
  #ggsave(path = outDir, filename = 'P_values.pdf', plot = p)

  return(list('results' = res, 'p.vals_plot' = p))
}


#' @noRd

plot.Pvals = function(res, showLabels = F, sigLvl = 0.05, main='P_values (ANOVA)'){

  labels = rownames(res)
  X = factor(1:nrow(res))
  sig = res[,"adj.p_val"] <= sigLvl
  cols = ifelse(sig, "steelblue","navy")
  #cols = factor(ifelse(sig, "significant","non-significant"))
  Y = -log10(res[,"adj.p_val"])
  pltTemp = data.frame(X=X, Y=Y, cols=cols)
  rownames(pltTemp) = labels

  p = ggplot(data=pltTemp, aes(x=X, y=Y))+
    geom_point(col=cols)+
    geom_abline(intercept = -log10(sigLvl), slope = 0, col='red')

  if(showLabels) p = p + ggrepel::geom_text_repel(aes(x = X, y=Y, label=rownames(res)))

  p = p+
    theme_bw()+
    ggtitle('P-values')+
    xlab('Bins')+
    ylab('-log10( p-value )')+
    theme(axis.ticks.x=element_blank(),
          axis.text.x = element_blank(),
          plot.title = element_text(hjust = 0.5))
  p
}

#' @noRd
plot_anova = function(res, main='Significant Bins (ANOVA)', showLabels = F){
  sigLevel <- NULL
  labels = rownames(res)
  X = 1:nrow(res)
  sig = res[,2] <= sigLevel
  cols = ifelse(sig, "steelblue","navy")

  Y = -log10(res[,2])
  plot(X,Y, col=cols, pch=16, main=main, xlab='Bins', ylab="-log10(p-val)")
  if(showLabels) {
    Xsig = X[sig]
    Ysig = Y[sig]
    text(Xsig+0.02*max(Xsig),Ysig+0.02*max(Ysig), labels = rownames(res)[sig])
  }
  abline(h=-log10(sigLevel))
}

#' @noRd
do_anova_m = function(data, groups, adjustMethod='fdr', thresh=0.05){
  # Performs multiple anova tests (one on each variable in the data)
  aov.res = apply(data,2,function(x) aov(x~groups))
  anova.res = lapply(aov.res, anova)
  res<-do.call('rbind', lapply(anova.res, function(x) { c(x["F value"][1,], x["Pr(>F)"][1,])}))
  res = cbind(res,p.adjust(res[,2], method = adjustMethod))
  colnames(res) = c('F-stat','p-val','adj.p_val')

  res = as.data.frame(res)
  sgnf_filter = which(res['adj.p_val'] <= thresh)
  posthoc.res = lapply(aov.res[sgnf_filter], TukeyHSD, conf.level=1-thresh)
  tukey_pvals = extract.pVals.Tukey(posthoc.res, thresh)

  return(list(anova_pvals = res, posthoc_res = posthoc.res, tukey_pvals = tukey_pvals, anova_res = anova.res))
}

#' @noRd
extract.pVals.Tukey = function(tukey.res, thresh=0.05){
  do.call('rbind', lapply(tukey.res, function(x) x[1][[1]][,'p adj']))
}


#' @noRd
print_anova_res = function(res, sigLvl = 0.05){
  anovaSig = sum(res$anova_pvals[,3] <= sigLvl)
  cat(sprintf('Total number of bins:\t%s \n', nrow(res$anova_pvals)))
  cat(sprintf('Significant bins:\t%s\n\n', anovaSig))

  cat(sprintf('Tukey post-hoc analysis:\n'))
  cat(sprintf('Comp. groups \t sig. bins\n'))
  if (ncol(res$tukey_pvals) != 0) {
    for(i in 1:ncol(res$tukey_pvals)){
      cat(sprintf('\t %s \t\t %d\n', colnames(res$tukey_pvals)[i], sum(res$tukey_pvals[i,]<=sigLvl)))
    }
  }
}

