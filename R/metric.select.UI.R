#' Indicator metric selection
#'
#' Determines which indicator metrics which best differentiate the test site from its nearest-neighbbour reference sites. Metrics that indicate impairment will be
#' used preferentially.
#'
#' A interative selection algorithm is used as follows:
#'
#' 1. The first metric selected for the final set is the one which displayes the greatest distance from the Reference condition mean
#'
#' 2. Metrics with a pearson correlation greater than 0.7 to (any of) the selected metric(s) are excluded from further steps
#'
#' 3. The ranked departure of remaining metrics is divided by the (maximum) correlation with the metric(s) previously included in the analysis
#'
#' 4. The metric with the greatest score is selected for inclusion in the final set
#'
#' 5. Return to step 2 until the number of selected metrics is equal to the greater of 4 or 1/5 the number of Reference sites
#'
#' If no metrics or too few metrics demonstrate impairment, the following metrics are included until the maximum is reached:
#' Richness, Percent Dominance, HBI, Percent EPT.
#'
#' @param Test Vector containing metric scores at the test site. Should be a single row from \code{benth.met} or \code{add.met}.
#' @param Reference Data frame of metric scores at the reference sites. Should be output from \code{benth.met} or \code{add.met}.
#' @param Rank Use rank differences in metric selection
#' @param outbound Used if outlier.rem=T A numeric value between 0 and 1 indicating the outlier boundary for defining values as final outliers (default to 0.1)
#' @return $Best.Metrics - Vector containing the final selected indicator metrics
#' @return $Indicative.Metrics - Vector containing all metrics that indicate impairment
#' @return $raw.data - Data frame containing only selected best metrics
#' @return $ref.sites - Vector containing input reference site names
#' @return $outlier.ref.sites - Vector containing sites removed as potential outliers
#' @keywords Benthic Metrics
#' @export
#' @examples
#' data(YKBioData,envir = environment())
#' bio.data<-benth.met(YKBioData,2,2)$Summary.Metrics
#' nn.refsites<- c("075-T-1", "019-T-1","003-T-1","076-T-1","071-T-1","022-T-1","074-T-1",
#' "002-T-1","004-T-1","073-T-1","186-T-1","062-T-1","005-T-1","025-T-1",
#' "187-T-1","023-T-1","193-T-1","192-T-1","196-T-1","194-T-1")
#' metric.select(bio.data[201,],bio.data[nn.refsites,])

metric.select.UI <- function (Test,Reference,outlier.rem=T,rank=F, outbound=0.1) {
  Reference<-Reference[rowSums(is.na(Reference))!=ncol(Reference), !is.na(Test)]
  Test<-Test[,!is.na(Test)]
  
  raw.data1<-rbind(Reference,Test)
  raw.data<-tsa.zscore(Test=Test,Reference=Reference)
  #raw.data[is.nan(as.matrix(raw.data))]<-0
  
  if (outlier.rem==T) {
    Reference<-Reference[,which(apply(Reference, 2, mad)!=0 & !is.na(apply(Reference, 2, mad)!=0))]
    Reference<-Reference[c(which(pcout(Reference,outbound=outbound)$wfinal01==1)),]
    raw.data<-tsa.zscore(Test[,colnames(Test)%in%colnames(Reference)],Reference[,colnames(Reference)%in%colnames(Test)])
  }
  
  nRef<-nrow(raw.data)-1
  nInd<-ncol(raw.data)
  test.var<-NULL
  
  #This loops through each metric and removes metrics that have fewer than 25% unique variables
  restricted.metrics<-NULL
  for (i in 1:nInd) {
    if ((length(unique(raw.data[1:nRef,i]))>1) & !(max(table(raw.data[1:nRef,i]))>((1/3)*(nRef))) & (IQR(raw.data[1:nRef,i])>0)){ #& !has_warning(!has_error(invisible(covMcd(data[,i]))))){ #ceiling(nrow(tsa3)*0.1)
      restricted.metrics[i]<-paste0(colnames(raw.data)[i])
    } else {
      next
    }
  }
  data<-raw.data[,colnames(raw.data) %in% restricted.metrics]
  
  #indicative.metrics<-c(hb[which(hb%in%colnames(data)[which(data[nrow(data),]>0)])],lb[which(lb%in%colnames(data)[which(data[nrow(data),]<0)])])
  indicative.metrics<-colnames(data)
  indicative.metrics<-indicative.metrics[!duplicated(indicative.metrics)]
  reduced.data<-data[,colnames(data)%in%indicative.metrics]
  
  if (length(indicative.metrics)>1) {
    
    if (rank==T){
      diff<-length(reduced.data):1
    } else {
      diff<-sort(as.vector(reduced.data[nrow(reduced.data),]),decreasing=T)
    }
    names(diff)<-names(sort(abs(reduced.data[nrow(reduced.data),]),decreasing=T))
    test.var<-names(diff[1])
    cors<-abs(cor(reduced.data[1:nRef,],method="p")[,test.var])+0.001
    
    if (any(cors<0.75)) {
      cors<-cors[which(cors<0.75)]
      cors<-cors[(names(diff))]
      cors<-cors[which(!is.na(cors))]
      
      test.var[2]<-names(sort((diff[names(cors)]/cors),decreasing=T)[1])
      #[which(names(sort(abs(diff/cor(reduced.data[1:(nrow(reduced.data)-1),],method="k")[,names(test.var)]),decreasing=T))%in%(names(cor(reduced.data[1:(nrow(reduced.data)-1),],method="k")[,names(test.var)])[which(cor(reduced.data[1:(nrow(reduced.data)-1),],method="k")[,names(test.var)]<0.7)]))][1]
      
      for (var in 1:(min((length(indicative.metrics)-2),ceiling(1/5*nrow(data))-2))) {
        cors<-apply((abs(cor(reduced.data[1:nRef,],method="p")[,test.var])+0.001),1,max)
        if (any(cors<0.75)) {
          cors<-cors[which(cors<0.75)]
          cors<-cors[(names(diff))]
          cors<-cors[which(!is.na(cors))]
          test.var[var+2]<-names(sort((diff[names(cors)]/cors),decreasing=T)[1])
        } else {
          break
        }
      }
    }
  }
  
  if ((!is.null(test.var)) & (length(test.var)<max(4,ceiling(1/5*nrow(data))))){
    if (("O:E" %in% colnames(raw.data)) & !("O:E" %in% test.var)) {
      test.var<-c(test.var,"O:E")
    }
    if (("Percent.Dominance" %in% colnames(raw.data)) & !("Percent.Dominance" %in% test.var)) {
      test.var<-c(test.var,"Percent.Dominance")
    }
    if (("Richness" %in% colnames(raw.data)) & !("Richness" %in% test.var)) {
      test.var<-c(test.var,"Richness")
    }
    if (("Percent.mEPT" %in% colnames(raw.data)) & !("Percent.mEPT" %in% test.var)) {
      test.var<-c(test.var,"Percent.mEPT")
    }
    if (("Percent.Chironomidae" %in% colnames(raw.data)) & !("Percent.Chironomidae" %in% test.var)) {
      test.var<-c(test.var,"Percent.Chironomidae")
    }
  }
  
  if (is.null(test.var)){
    if (length(indicative.metrics)==1) {
      test.var<-indicative.metrics
    } else {test.var<-NULL}
    if (("O:E" %in% colnames(data)) & !("O:E" %in% test.var)) {
      test.var<-c(test.var,"O:E")
    }
    if (("Percent.Dominance" %in% colnames(data)) & !("Percent.Dominance" %in% test.var)) {
      test.var<-c(test.var,"Percent.Dominance")
    }
    if (("Richness" %in% colnames(data)) & !("Richness" %in% test.var)) {
      test.var<-c(test.var,"Richness")
    }
    if (("Percent.mEPT" %in% colnames(data)) & !("Percent.mEPT" %in% test.var)) {
      test.var<-c(test.var,"Percent.mEPT")
    }
    if (("Percent.Chironomidae" %in% colnames(data)) & !("Percent.Chironomidae" %in% test.var)) {
      test.var<-c(test.var,"Percent.Chironomidae")
    }
  }
  
  test.var<-test.var[1:max(3,ceiling((1/5)*nRef))]
  test.var<-test.var[!is.na(test.var)]
  
  metric.auto<-NULL
  metric.auto$Best.Metrics<-test.var
  #metric.auto$Indicative.Metrics<-indicative.metrics
  metric.auto$raw.data<-raw.data1[rownames(raw.data),colnames(raw.data1) %in% test.var]
  metric.auto$ref.sites<-rownames(raw.data1)[-c(nrow(raw.data1))]
  if (outlier.rem){
    metric.auto$outlier.ref.sites<-rownames(raw.data1)[-c(nrow(raw.data1))][!rownames(raw.data1)[-c(nrow(raw.data1))]%in%rownames(Reference)]
    metric.auto$raw.data.with.outliers<-raw.data1[,colnames(raw.data1) %in% test.var]
  }
  class(metric.auto)<-"met.sel"
  return(metric.auto)
}

