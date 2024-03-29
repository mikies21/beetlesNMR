% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NMRMetab_anova.R
\name{NMRMetab_anova}
\alias{NMRMetab_anova}
\title{Normalisation and Scaling}
\usage{
NMRMetab_anova(
  data,
  index_col = 3,
  group_test = "group",
  sigLevel = 0.05,
  adjMethod = "fdr"
)
}
\arguments{
\item{data}{a data.frame. Column as variable and rows as sample}

\item{index_col}{integer value. The column number with first bin/matabolite}

\item{group_test}{charachter string. name of the group column to test}

\item{sigLevel}{double. set to 0.05, significance level to detect differences between groups,}

\item{adjMethod}{charachter string. multiple comparison adjustment method. Available methods: "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none". "BH" stands for Benjamini-Hochberg; "BY" for Benjamini-Yakutieli; "fdr" is false discovery rate (same as "BH")}
}
\description{
Normalisation and Scaling
}
\author{
Eva Caamano Gutierrez

Arturas Grauslys
}
