#' Reference specific Benthic Metric Calculation
#'
#' Calculation of a variety of metrics for determing impairment of Benthic Macroinvertebrate Communities, based on available Reference Condition information.
#' @param Test.taxa Output of \code{\link[BenthicAnalysis]{benth.met}}$Raw.Data for the test site
#' @param Reference.taxa Output of \code{\link[BenthicAnalysis]{benth.met}}$Raw.Data for the Reference site
#' @return Indicator metrics caluclated from \code{\link[BenthicAnalysis]{benth.met}} in addition to O:E ratios and Bray-Curtis Distance.
#' @keywords Benthic Metrics
#' @export
#' @examples
#' #load datasets
#' data(YKEnvData,envir = environment()) #Biological dataset
#' data(YKBioData,envir = environment()) #Environmental dataset
#'
#' #Calculate indicator metrics from raw biological data
#' bio.data.test<-benth.met(YKBioData,2,2)
#'
#' #standardize row names between datasets
#' rownames(YKEnvData)<-rownames(bio.data.test$Site.List)
#'
#' #Match a test site (#201) to the nearest neighbour reference set
#' nn.sites<-site.match(YKEnvData[201,],YKEnvData[1:118,],k=NULL,adaptive=T)
#'
#' #Extract the raw taxa data for additional metric calculation
#' taxa.data<-rbind(bio.data.test$Raw.Data[names(nn.sites$final.dist),],bio.data.test$Raw.Data[rownames(bio.data.test$Raw.Data[201,]),])
#'
#' #Calculate additional metrics based on nearest neighbour reference sites
#' additional.metrics<-add.met(Test=bio.data.test$Raw.Data[rownames(bio.data.test$Raw.Data[201,]),],Reference=bio.data.test$Raw.Data[names(nn.sites$final.dist),])
#' additional.metrics

add.met<-function (Test,Reference,original=F,tax.fields=2) {
  if (any(colnames(Test)%in%colnames(Reference)==F)){
    stop("Column name mismatch between Test and Reference Set")
  }
  
  raw.data<-rbind(Reference,Test)
  
  nRef<-nrow(Reference)
  
  pRef<-colSums(vegan::decostand(raw.data[rownames(Reference),],"pa"))/nrow(Reference)
  e<-adapt.sum1(pRef[names(which(pRef>=0.5))])
  e.var<-apply(vegan::decostand(raw.data[rownames(Reference),names(which(pRef>=0.5))],"pa"),1,function(x) adapt.sum1(x)/e)
  o<-adapt.sum1(vegan::decostand(raw.data[rownames(Test),names(which(pRef>=0.5))],"pa"))/e
  e.o<-c(e.var,o)
  e.o.stand<-(e.o-mean(e.o[1:nRef]))/sd(e.o[1:nRef])
  
  ref.bray<-rowMeans(braydist(raw.data[rownames(Reference),]))
  test.bray<-rowMeans(braydist(raw.data))[nrow(raw.data)]
  bray<-c(ref.bray,test.bray)
  bray.stand<-(bray-mean(bray[1:nRef]))/sd(bray[1:nRef])
  
  tax.names<-as.data.frame(t(matrix(c(lapply(colnames(raw.data), function(x) substr(x, start=1,stop=(gregexpr(pattern =';',x)[[1]][1]-1))),
                                      lapply(colnames(raw.data), function(x) substr(x, start=(gregexpr(pattern =';',x)[[1]][1]+1),stop=nchar(x)))),ncol=2)))
  colnames(tax.names)<-colnames(raw.data)
  raw.data<-rbind(tax.names,raw.data)
  raw.data<-cbind(rownames(raw.data),raw.data)
  
  ca.ord<-cca(log(rbind(Reference,Test)[,names(which(pRef>=0.1))]+1))
  #ca1<-c(ca.ord$CA$u[,1],predict(ca.ord,log(Test[,names(which(pRef>=0.1))]+1),type="wa")[1])
  #names(ca1)[nRef+1]<-rownames(Test)
  #ca2<-c(ca.ord$CA$u[,2],predict(ca.ord,log(Test[,names(which(pRef>=0.1))]+1),type="wa")[2])
  #names(ca2)[nRef+1]<-rownames(Test)
  ca1<-ca.ord$CA$u[,1]
  ca2<-ca.ord$CA$u[,2]
  
  if (original==T){
    raw.data<-data.frame(cbind(benth.met(x=raw.data,tax.fields=2,site.fields=1)$Summary.Metrics,c(e.var,o),c(ref.bray,test.bray),ca1,ca2))
    colnames(raw.data)[(ncol(raw.data)-3):ncol(raw.data)]<-c("O:E","Bray-Curtis","CA1","CA2")
  }
  if (original==F){
    raw.data<-data.frame(cbind(c(e.var,o),c(ref.bray,test.bray),ca1,ca2))
    colnames(raw.data)<-c("O:E","Bray-Curtis","CA1","CA2")
  }
  
  return(raw.data)
}
