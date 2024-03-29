% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NMRMetab_CRS.R
\name{NMRMetab_CRS}
\alias{NMRMetab_CRS}
\alias{NMRMetabCRS}
\title{CRS function}
\usage{
NMRMetabCRS(dat, metabolite_list, correlation_type = "pearson", save_csv = F)
}
\arguments{
\item{dat}{a data.frame. Column as variable and rows as sample}

\item{metabolite_list}{a data.frame. with rows as metabolites.. column names are 'HMDB' and 'Metab'}

\item{correlation_type}{character string. type of correlation to use for CRS score 'preason'(default) or 'spearman'}

\item{save_csv}{boolean. save all into result into file}
}
\description{
This function will take in either a object of class data.frame and a using a second dataframe with a list of all the metabolite identifed, will return the CRS score for all the bins. this could work with other metabolite identifiers; for example bins with the name of the     metabolite, as  log as they are unique. mistakes could happen when greple take NAD and NADH
}
\author{
Michele Fresneda Alarcon
}
