#--------------------------------------------------------------------------------#
#---------- Univariate 2 group statistical test for NMR metabolomics data -------#
#--------------------------------------------------------------------------------#
#
# Dr Eva Caamano Gutierrez & Dr Arturas Grauslys, 2017
# This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.
# https://creativecommons.org/licenses/by-nc-sa/4.0/

# The script's aim is to allow unexperienced R users to perform univariate tests in NMR metabolomics
# data prepared such as rows contain samples and columns variables
# furthermore the first two columns are data describers: row names and 2 groups
# There are different options that the user must specify including:
# (a) whether to assume normality (T/F)
# (b) whether to assume equal variance (T/F)
# (c) whether a paired test should be implemented (T/F)
# Outputs a file with confidence intervals and adjusted p-values by Benjamini-Hochberg and Bonferroni
# In order to use this script follow this steps:
# (1) Source this file
# (2) Call the NMRMetab_UnivariateTest() function with adequate parameters
# e.g. NMRMetab_UnivariateTest(myBucketTable,normality=T,equal.variance=F,paired=F)
# These options can be changed and/or expanded but any aditions should be
# documented and dated in this preface.
#--------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

#' @title Univariate test
#' @name NMRMetab_foldchange
#' @export
#' @description This function will take in either a object of class data.frame and a using a second dataframe with a list of all the metabolite identifed, will return the CRS score for all the bins. this could work with other metabolite identifiers; for example bins with the name of the     metabolite, as  log as they are unique. mistakes could happen when greple take NAD and NADH
#' @author Eva Caamano Gutierrez
#' @author Arturas Grauslys
#' @param data a data.frame. Column as variable and rows as sample
#' @param groupID string. name off the variable containing your grouping details must only have 2 facors eg Healthy vs Diseases
#' @param index_col index colum oif first mobservation
#' @param divisorID string. what group you're dividing by
#' @param paired boolean
#' @param normality boolean
#' @param equal.variance boolean

NMRMetab_UnivarTest<-function(data, groupID,index_col = 3, paired=F,normality=T,equal.variance=F){


    #This function is created for only 2 groups
    dat<-data[,3:ncol(data)]
    groups<-unique(data[,2])

    pvals<-c()
    metabolite<-c()
    ConfInt<-c()

    for(i in 1:ncol(dat)){
      g1<-dat[data[,2]==groups[1],i]
      g2<-dat[data[,2]==groups[2],i]

      if(normality==T){
        if(equal.variance==F){
          if(paired==F){
            a<-t.test(g1,g2,var.equal=F,paired=F,na.action=na.rm)
            if(i==1){
              print("Doing a two-sided independent samples Welch Test")
            }
          }else{
            a<-t.test(g1,g2,var.equal=F,paired=T,na.action=na.rm)
            if(i==1){
              print("Doing a two-sided paired Welch Test")
            }
          }
        }else{
          if(paired==F){
            a<-t.test(g1,g2,var.equal=T,paired=F)
            if(i==1){
              print("Doing a two-sided independent samples Student's t-test")
            }
          }else{
            a<-t.test(g1,g2,var.equal=T,paired=T)
            if(i==1){
              print("Doing a two-sided paired Student's t-test")
            }
          }
        }
      }else{
        if(equal.variance==F){
          if(paired==F){
            a<-ks.test(g1,g2)
            #aa<-t.test(g1,g2,var.equal = F,paired=F)
            if(i==1){
              print("Doing a Kolmogorov-Smirnov Test")
            }
          }else{
            a<-wilcox.test(g1,g2,paired=T,conf.int = T)
            if(i==1){
              print("Doing a Wilcoxon's matched pairs test")
            }
          }
        }else{
          if(paired==F){
            a<-wilcox.test(g1,g2,var.equal=T,paired=F,conf.int = T)
            print("Doing a two-sided Mann-Whitney U-test")
          }else{
            a<-wilcox.test(g1,g2,paired=T,conf.int = T)
            print("Doing a Wilcoxon's matched pairs test test")
          }
        }
      }

      pvals<-c(pvals,a$p.value);
      metabolite<-c(metabolite,colnames(dat)[i])
      if(is.null(a$conf.int)){
        #ConfInt<-c(ConfInt,paste("(",round(aa$conf.int[1],5),",",round(aa$conf.int[2],5),")",sep=""))
        ConfInt<-c(ConfInt,"Not in this test")
        #ConfInt<-c(ConfInt,NA)
      }else{
        ConfInt<-c(ConfInt,paste("(",round(a$conf.int[1],5),",",round(a$conf.int[2],5),")",sep=""))
      }
    }
    pvalsNA<-round((pvals),4)
    pvalsBH<-round(p.adjust(pvals,method="BH"),5)
    pvalsBonf<-round(p.adjust(pvals,method="bonferroni"),5)

    out<-data.frame(metabolite, ConfInt, pvalsNA,pvalsBH, pvalsBonf)
    colnames(out)<-c("Metabolite/Bucket","95% Conf. Int. Mean.","Unadjusted pvals","BH pvals","Bonf pvals")
    return(out)
  }


