#' Z-Score Calculation
#'
#' This calculates z-scores centered on Reference Sites. A certain number of NAs can be permitted in the reference set for each metric but
#' metrics with NAs at the test site are automatically removed
#' @param Test Data frame of metric scores at the test site.
#' @param Reference Data frame of metric scores at the reference sites
#' @param na.cutoff A value between 0-1 indicating the percent of the reference set that can contain NAs for any metric. NAs are replaced by the mean value.
#' @keywords z-score
#' @export
#' @examples
#' tsa.zscore()

tsa.zscore<-function(Test,Reference,na.cutoff=0.7) {
  #first check if metrics match between test site and reference set
  if (any((colnames(Test)==colnames(Reference))==F)|
      (ncol(Test)!=ncol(Reference))) {
    stop("Metric mismatch between test site and reference set")
  }
  tsa<-rbind(Reference,Test)
  tsa<-t(apply(tsa,1,function(x) (x-colMeans(Reference,na.rm=T))/apply(Reference,2, function(x) sd(x,na.rm = T))))
  if (any(is.na(tsa))) {
    tsa<-tsa[,!is.na(tsa[nrow(tsa),])]
    tsa<-tsa[,apply(apply(tsa,2,is.na),2,function(x) length(which(x==T))<(na.cutoff*nrow(tsa)-1))]
    if (any(apply(tsa,2,function(x) any(is.na(x))))){
      for (i in names(which(apply(tsa,2,function(x) any(is.na(x)))))){
        tsa[is.na(tsa[,i]),i]<-mean(Reference[!is.na(Reference[,i]),i])
      }
    }
  }
  
  return(tsa)
}

