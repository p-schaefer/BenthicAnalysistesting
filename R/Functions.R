
#' Test Site Analysis Function
#'
#' Implimentation of Bowman and Somers (2006) Test Site Analysis (\url{http://goo.gl/h4JAGP}).
#' For implimenation as in Bowman and Somers (2006) use distance=NULL, outlier.rem=F, m.select=F. This function allows optional funcationality including automatic selection of
#' indicator metrics using \code{\link[BenthicAnalysis]{metric.select}}, as well as outlier removal from the reference set using \code{\link[mvoutlier]{arw}} or \code{\link[mvoutlier]{pcout}}.
#' An alternative mahalnobis distance calculation is also implimented using weighted means and covariance matrix. The weights are supplied from
#' the ecological distance between test sites and reference sites using \code{\link[BenthicAnalysis]{site.match}}.
#'
#' @param Test Data frame of metric scores at the test site.
#' @param Reference Data frame of metric scores at the reference sites.
#' @param distance Vector of weights to use for optional weighted Mahalanobis Distance calculation. Can be output of \code{\link[BenthicAnalysis]{site.match}}$final.dist,
#' or NULL for unweighted Mahalanobis Distance Calculation
#' @param outlier.rem Logical argument indicating whether \code{\link[mvoutlier]{pcout}} should be used to remove outliers from the Reference set.
#' @param m.select Logical argument indicating whether the best subset of metrics should be automatically selected from the data.
#' @param na.cutoff A value between 0-1 indicating the percent of the reference set that can contain NAs for any metric. NAs are replaced by the mean value.
#' @param outbound Used if outlier.rem=T A numeric value between 0 and 1 indicating the outlier boundary for defining values as final outliers (default to 0.1)
#' @return $general.results - Table containing the number and names of reference sites used and number of indicator metrics used
#' @return $tsa.results - Table containing the status rank of the test site, numerical interval and equivalence test,
#' test site's mahalanobis sistance, the upper and lower critical mahalanobis distance values, the non-centrallity parameter lambda and F-value of the test site
#' @return $jacknife - Table conataining the jacknifed consistency of the status rank and  95% Jacknife interval for test site's mahalanobis distance
#' @return $partial.tsa - Only present if the test site is ranked impaired or possibly impaired. Table containing p and F values for each metric's
#' contribution to the overall mahalanobis distance score
#' @return $mahalanobis.distance - vector of reference sites and test site mahalanobis distance scores
#' @return $z.scores - Table containing test and reference sites metrics standrdized to z-scores of the reference sites only
#' @keywords Test Site Analysis, Mahalanobis Distance
#' @references \url{http://goo.gl/h4JAGP}
#' @export
#' @examples
#' #load datasets
#' data(YKEnvData,envir = environment()) #Biological dataset
#' data(YKBioData,envir = environment()) #Environmental dataset
#'
#' #Calculate indicator metrics from raw biological data
#' bio.data.test<-benth.met(YKBioData,2,2)
#'
#' #Extract just the summary metrics
#' bio.data<-bio.data.test$Summary.Metrics
#'
#' #standardize row names between datasets
#' rownames(YKEnvData)<-bio.data.test$Site.List
#'
#' #Match a test site (#201) to the nearest neighbour reference set
#' nn.sites<-site.match(YKEnvData[201,-c(1)],YKEnvData[1:118,-c(1)],k=NULL,adaptive=T)
#'
#' #Calculate additional metrics based on selected Reference sites
#' taxa.data<-add.met(Test=bio.data.test$Raw.Data[201,],Reference=bio.data.test$Raw.Data[names(nn.sites$final.dist),])
#'
#' #TSA test of indicator metrics at test site and reference sites selected used site.match()
#' tsa.results<-tsa.test(Test=taxa.data[nrow(taxa.data),],Reference=taxa.data[names(nn.sites$final.dist),],distance=nn.sites$final.dist, outlier.rem=T, m.select=T)
#' tsa.results
#'
#' #Evaluate Results
#' boxplot(tsa.results)
#' plot(tsa.results)
#' pcoa.tsa(tsa.results)


tsa.test<- function(Test, Reference, distance=NULL, outlier.rem=F, m.select=F,rank=F,na.cutoff=0.7,outbound=0.1) {
  Reference<-Reference[rowSums(is.na(Reference))!=ncol(Reference), ]
  raw.data1<-rbind(Reference,Test)
  
  
  if (any((colnames(Test)%in%colnames(Reference))==F)|
      (ncol(Test)!=ncol(Reference))) {
    stop("Metric mismatch between test site and reference set")
  }
  if(!is.null(distance)) {
    if (any(names(distance)%in%rownames(Reference))==F) {
      stop("Missing one or more ecological distances of Reference sites")
    } else {
      distance<-distance[names(distance)%in%rownames(Reference)]
    }
  }
  
  Reference<-Reference[,!is.na(Test)]
  Reference<-Reference[,apply(apply(Reference,2,is.na),2,function(x) length(which(x==T))<(na.cutoff*nrow(Reference)))]
  if (any(apply(Reference,2,function(x) any(is.na(x))))){
    for (i in names(which(apply(Reference,2,function(x) any(is.na(x)))))){
      Reference[is.na(Reference[,i]),i]<-mean(Reference[!is.na(Reference[,i]),i])
    }
  }

  #Reference<-Reference[,!is.na(colSums(Reference))]
  #Reference<-Reference[,!is.infinite(colSums(Reference))]
  #Reference<-Reference[,!is.na(colSums(Reference))]
  Test<-Test[,colnames(Reference)]
  #if (!is.null(add.metric) & rownames(Test)%in%rownames(add.metric) & any(rownames(Reference)%in%rownames(add.metric)==F)){
  #  stop("Missing add.metrics data for one or more test and/or reference sites")
  #}

  data.raw<-rbind(Reference,Test)

  nRef<-nrow(Reference)
  nInd<-ncol(data.raw)

  if (m.select==F & outlier.rem==F){
    if (ncol(data.raw)>(2*nrow(data.raw))) {
      stop("Too many metrics selected for number of reference sites")
    }
    cc <- try(mahalanobis(Reference,colMeans(Reference),cov(Reference),inverted=F), silent=T)
    if(is(cc,"try-error")) {
      stop("Matrix singularity found. One or more selected metrics are linear combinations of one another. Remove redundant metrics or use m.select=T")
    } else {
      data<-rbind(Reference,Test)
    }
  }

  if (m.select==T){
    metric.select.object<-metric.select(Test=data.raw[nrow(data.raw),],Reference=data.raw[1:nRef,],outlier.rem=outlier.rem,rank=rank,outbound=outbound)
    data<-metric.select.object$raw.data
  }

  if (m.select==F & outlier.rem==T) {
    data<-data.raw
    if (((nrow(data))-1-(ncol(data)*2))<1){
      stop("Too few reference sites for number of indicators")
    }
    if (!any(apply(Reference, 2, mad)!=0 & !is.na(apply(Reference, 2, mad)!=0))){
      stop("More than 50% equal values in one or more variables or too many NAs")
    }
    Reference1<-Reference[c(which(pcout(Reference,outbound=outbound)$wfinal01==1)),]
    data<-rbind(Reference1,Test)
    #if (has_warning(covMcd(data[1:nRef,])) & any(findCorrelation(cor(data[1:nRef,]),cutoff=0.8))) {
    #  data<-data[,-c(findCorrelation(cor(data[1:nRef,]),cutoff=0.8))]
    #}
    #if (has_warning(covMcd(data[1:nRef,])) & any(findCorrelation(cor(data[1:nRef,]),cutoff=0.7))) {
    #  data<-data[,-c(findCorrelation(cor(data[1:nRef,]),cutoff=0.7))]
    #}
    #if (has_warning(covMcd(data[1:nRef,])) & any(findCorrelation(cor(data[1:nRef,]),cutoff=0.6))) {
    #  data<-data[,-c(findCorrelation(cor(data[1:nRef,]),cutoff=0.6))]
    #}
    #if (has_warning(covMcd(data[1:nRef,])) & any(findCorrelation(cor(data[1:nRef,]),cutoff=0.5))) {
    #  data<-data[,-c(findCorrelation(cor(data[1:nRef,]),cutoff=0.5))]
    #}
    #data<-data[arw(data[1:nRef,],suppressWarnings(covMcd(data[1:nRef,])$center),suppressWarnings(covMcd(data[1:nRef,])$cov))$w,]
    #if (rownames(data)[nrow(data)]!=rownames(Test)){
    #  data[(nrow(data)+1),]<-Test[,colnames(Test)%in%colnames(data)]
    #}
  } 

  nRef<-nrow(data)-1
  nInd<-ncol(data)

  if ((nRef-(nInd*2))<1){
    stop("Too few reference sites for number of indicators")
  }

  if (is.null(distance)) {
    tsa.ref.cent<-colMeans(data[1:nRef,])
    tsa.cov<-cov(data[1:nRef,])
  } else {
    tsa.ref.cent<-cov.wt(data[1:nRef,],wt=as.numeric(1/distance[rownames(data[1:nRef,])]),center=T)$center
    tsa.cov<-cov.wt(data[1:nRef,],wt=as.numeric(1/distance[rownames(data[1:nRef,])]),center=T)$cov
  }

  tsa.dist<-mahalanobis(data,tsa.ref.cent,(tsa.cov),inverted=F)
  tsa.lambda<-qchisq(0.05,nInd, ncp = 0, lower.tail = FALSE, log.p = FALSE)*nRef
  tsa.F<-((nRef-nInd)*nRef*tsa.dist[nRef+1])/(nInd*(nRef-1))
  #tsa.F<-(nRef*tsa.dist[nRef+1]*(nRef-nInd))/(nInd*(nRef-1)*(nRef+1))
  tsa.NCPinterval<-1-pf(tsa.F, nInd, (nRef-nInd), tsa.lambda, log=FALSE)
  tsa.NCPequivalence<-1-tsa.NCPinterval
  tsa.impairment<-if(tsa.NCPinterval<0.05){
    "Impaired"
  } else {
    if(tsa.NCPinterval>0.95){
      "Not Impaired"
    } else {
      "Possibly Impaired"
    }
  }

  if (tsa.impairment!="Not Impaired"){
    part.tsa<-data.frame(matrix(nrow=nInd,ncol=2))
    rownames(part.tsa)<-colnames(data)
    colnames(part.tsa)<-c("p","F")

    for (n in 1:(nInd)) {
      part.tsa.stand<-data[,-c(n)]
      part.tsa.cov<-tsa.cov[-c(n),-c(n)]
      part.tsa.ref.cent<-tsa.ref.cent[-c(n)]

      part.tsa.dist<-mahalanobis(part.tsa.stand[nrow(part.tsa.stand),],part.tsa.ref.cent,(part.tsa.cov),inverted=F)
      part.tsa.F<-(nRef - nInd)*(tsa.dist[nrow(part.tsa.stand)] - part.tsa.dist)/(nRef - 1 + part.tsa.dist)
      part.tsa.NCPinterval<-1-pf(part.tsa.F, 1, (nRef-nInd))
      part.tsa[n,]<-c(part.tsa.NCPinterval,part.tsa.F)
    }
  }
  
  tsa.dist.r<-NULL
  for (i in 1:nRef) {
    data.r<-data[-c(i),]
    nRef.r<-nrow(data.r)-1
    nInd.r<-ncol(data.r)
    if (is.null(distance)) {
      tsa.ref.cent.r<-colMeans(data.r[1:nRef.r,])
      tsa.cov.r<-cov((data.r[1:nRef.r,]))
    } else {
      tsa.ref.cent.r<-cov.wt(data.r[1:nRef.r,],wt=as.numeric(1/distance[rownames(data.r[1:nRef.r,])]),center=T)$center
      tsa.cov.r<-cov.wt(data.r[1:nRef.r,],wt=as.numeric(1/distance[rownames(data.r[1:nRef.r,])][1:nRef.r]),center=T)$cov
    }

    tsa.dist.r[i]<-mahalanobis(data.r,tsa.ref.cent.r,(tsa.cov.r),inverted=F)[nrow(data.r)]

    if (i==nRef){
      tsa.F.r<-((nRef.r-nInd.r)*nRef.r*tsa.dist.r)/(nInd.r*(nRef.r-1))
      tsa.lambda.r<-qchisq(0.05,nInd.r, ncp = 0, lower.tail = FALSE, log.p = FALSE)*nRef.r
      tsa.NCPinterval.r<-1-pf(tsa.F.r, nInd.r, (nRef.r-nInd.r), tsa.lambda.r, log=FALSE)
      tsa.NCPequivalence.r<-1-tsa.NCPinterval.r
      tsa.CI<-quantile(tsa.dist.r,probs=c(0.025,0.975))

      if (tsa.impairment=="Not Impaired") {
        jacknife.consistency<-length(which(tsa.NCPinterval.r>=0.95))/nRef
      }
      if ((tsa.impairment=="Impaired")){
        jacknife.consistency<-length(which(tsa.NCPinterval.r<=0.05))/nRef
      }
      if ((tsa.impairment=="Possibly Impaired")){
        jacknife.consistency<-(length(which(tsa.NCPinterval.r>=0.05 & tsa.NCPinterval.r<=0.95)))/nRef
      }
    }
  }
  
  rand.F<-NULL
  rand.mah<-NULL
  for (i in 1:nrow(data)){
    rand.F[i]<-((nRef-nInd)*nRef*mahalanobis(data,colMeans(data[-c(i),]),(cov(data[-c(i),])),inverted=F)[i])/(nInd*(nRef-1))
    rand.mah[i]<-mahalanobis(data,colMeans(data[-c(i),]),(cov(data[-c(i),])),inverted=F)[i]
    if (i==nrow(data)){
      rand.p<-1-(length(rand.F[rand.F<=rand.F[i]])/length(rand.F))
    }
  }

  general.results<-data.frame(matrix(nrow=6,ncol=1))
  rownames(general.results)<-c("Test Site","Reference Set","Selected Indicator Metrics","Significant Indicator Metrics","Number of Metrics","Number of Reference Sites")
  colnames(general.results)<-""
  general.results[,1]<-c(paste0(rownames(data)[(nRef+1)]),
                         paste0(rownames(data)[1:nRef],collapse="",sep=", "),
                         paste0(colnames(data),collapse="",sep=", "),
                         if (tsa.impairment!="Not Impaired"){
                           paste0(rownames(part.tsa[part.tsa[,1]<0.05,]),collapse="",sep=", ")
                         } else {"None"},
                         ncol(data),
                         nRef)

  tsa.results<-data.frame(matrix(nrow=9,ncol=1))
  rownames(tsa.results)<-c("TSA Impairment","Interval Test","Equivalence Test","Randomization p value","Test Site D2",
                           "Lower Critical","Upper Critical",
                           "TSA Lambda","TSA F Value")
  colnames(tsa.results)<-""
  tsa.results[,1]<-c(tsa.impairment,
                     round(tsa.NCPinterval,3),
                     round(tsa.NCPequivalence,3),
                     round(rand.p,3),
                     round(tsa.dist[nRef+1],3),
                     round(((nInd*(nRef-1))*qf(0.05, df1=nInd, df2=(nRef-nInd), ncp=tsa.lambda, log=FALSE))/((nRef-nInd)*nRef),3),
                     round(((nInd*(nRef-1))*qf(0.95, df1=nInd, df2=(nRef-nInd), ncp=tsa.lambda, log=FALSE))/((nRef-nInd)*nRef),3),
                     round(tsa.lambda,3),round(tsa.F,3))

  jacknife<-data.frame(matrix(nrow=3,ncol=1))
  colnames(jacknife)<-" "
  rownames(jacknife)<-c("Jacknife Consistency","Jacknife 2.5%","Jacknife 97.5%")
  jacknife[,1]<-c(round(jacknife.consistency*100,2),round(tsa.CI[1],3),round(tsa.CI[2],3))


  output<-NULL
  output$tsa.results<-tsa.results
  output$jacknife<-jacknife
  if (tsa.impairment!="Not Impaired"){
    output$partial.tsa<-part.tsa
  }
  output$mahalanobis.distance<-tsa.dist
  if (m.select) {
    output$selected.metrics<-colnames(data)
  }
  output$randomization.mahalanobis<-rand.mah

  output$general.results<-general.results
  if (!is.null(distance)){
    output$ecological.distance<-distance[names(tsa.dist)[1:nRef]]
  } 
  
  output$z.scores<-tsa.zscore(Test=Test,Reference=Reference[rownames(data)[-c(nrow(data))],])
  output$outlier.rem<-outlier.rem
  output$m.select<-m.select
  output$ref.sites<-rownames(Reference)
  output$test.site<-rownames(Test)
  if (outlier.rem) {
    if (m.select) {
      output$raw.data<-metric.select.object$raw.data
      output$outlier.ref.sites<-metric.select.object$outlier.ref.sites
      output$raw.data.with.outliers<-raw.data1
    } else {
      output$raw.data<-raw.data1[rownames(Reference)[-c(nrow(data))]%in%rownames(data),]
      output$outlier.ref.sites<-rownames(Reference)[-c(nrow(data))][!rownames(Reference)[-c(nrow(data))]%in%rownames(data)]
      output$raw.data.with.outliers<-raw.data1
    }
    
  } else {
    if (m.select) {
      output$raw.data<-metric.select.object$raw.data
    } else {
      output$raw.data<-raw.data1[rownames(Reference)[-c(nrow(data))]%in%rownames(data),]
    }
  }
  class(output)<-"tsa.object"
  return(output)
}

#' @export
print.tsa.object<-function(tsa.object) {
  cat("TSA Results:")
  print(tsa.object$tsa.results)
  cat("\n")
  if (!is.null(tsa.object$partial.tsa)) {
    cat("Partial TSA:\n")
    print(tsa.object$partial.tsa)
    cat("\n")
  }
  cat("Jacknife Results:")
  print(tsa.object$jacknife)
  cat("\n")
  cat("Mahalanobis Distance:\n")
  print(tsa.object$mahalanobis.distance)
}

#' @export
summary.tsa.object<-function(tsa.object){
  cat("General Results:")
  print(tsa.object$general.results)
  cat("\n")
  cat("TSA Results:")
  print(tsa.object$tsa.results)
  cat("\n")
  if (!is.null(tsa.object$partial.tsa)) {
    cat("Partial TSA Results:")
    print(tsa.object$partial.tsa)
    cat("\n")
  }
  cat("Jacknife Results:")
  print(tsa.object$jacknife)
  cat("\n")
  cat("Mahalanobis Distance:")
  print(tsa.object$mahalanobis.distance)
  cat("\n")
  if (!is.null(tsa.object$ecological.distance)) {
    cat("Ecological Distance:")
    print(tsa.object$ecological.distance)
    cat("\n")
  }
  if (!is.null(tsa.object$selected.metrics)) {
    cat("Selected Metrics:")
    print(tsa.object$selected.metrics)
    cat("\n")
  }
  cat("z-scores:")
  print(tsa.object$z.scores)
  cat("\n")
}

##############################################################################################################
#' Test Site Analysis confidence 
#'
#' Assess the confidence of TSA results based on jacknifing entire reference set.
#'
#' @param match.object An object from \code{\link[BenthicAnalysis]{site.match}}
#' @param tsa.object An object from \code{\link[BenthicAnalysis]{tsa.test}}
#' @param Reference Metric scores of Reference Set
#' @param benth.metric An object from \code{\link[BenthicAnalysis]{benth.met}}
#' @param levels Percent of Reference sites to remove
#' @param reps Replicate jacknife samples to assess
#' @return $summary.results
#' @keywords Test Site Analysis, Jacknife
#' @export
#' @examples
#' tsa.conf()

tsa.conf<-function(match.object,tsa.object,Reference=NULL,benth.metric=NULL,levels=c(1,5,10,20,30,40),reps=100){
  if (!is.null(Reference)){
    if (any((rownames(Reference)%in%names(match.object$all.dist))==F)) {
      stop("Site mismatch between reference set and match.object")
    }
    if (nrow(Reference)<(length(match.object$all.dist))){
      stop("Reference set must be greater than number of reference sites selected in match.object")
    }
  }
  
  if (is.null(Reference) & is.null(benth.metric)){
    stop("Either Reference or benth.metric must be defined")
  }
  
  if (!is.null(Reference) & is.null(benth.metric)){
    if (any(!rownames(tsa.object$raw.data[1:(nrow(tsa.object$raw.data)-1),])%in%rownames(Reference)==T)){
      Reference<-rbind(Reference,tsa.object$raw.data[!rownames(tsa.object$raw.data[1:(nrow(tsa.object$raw.data)-1),])%in%rownames(Reference),])
    }
  }
  
  full.results<-data.frame(matrix(nrow=(length(levels)*reps),ncol=10))
  colnames(full.results)<-c("Level","Rep","Impairment Rank","Test Site D2","TSA F-value","Randomization p-value",
                            "Number of Reference Sites","Number of Metrics","Selected Metrics","Significant Metrics")
  full.results[,"Level"]<-rep(levels,each=reps)
  full.results[,"Rep"]<-rep(1:reps,length(levels))
  
  summary.results<-data.frame(matrix(nrow=length(levels),ncol=6))
  rownames(summary.results)<-as.character(levels)
  colnames(summary.results)<-c("Rank Consistancy","Percent worse rank","Percent better rank","Mean D2 difference",
                               "Mean TSA F difference","Mean Randomization p difference")
  
  detailed.results<-data.frame(matrix(nrow=length(levels),ncol=20))
  rownames(detailed.results)<-levels
  colnames(detailed.results)<-c("Most selected metrics","Most significant metrics",
                                "Mean D2", "Median D2","25p D2","75p D2","5p D2","95p D2",
                                "Mean F", "Median F","25p F","75p F","5p F","95p F",
                                "Mean p", "Median p","25p p","75p p","5p p","95p p")
  
  if (!is.null(benth.metric)){
    taxa.avail<-TRUE
    ref.taxa.data<-benth.metric$Raw.Data[names(match.object$all.dist),]
    test.taxa.data<-benth.metric$Raw.Data[tsa.object$general.results[1,],]
  } else {
    taxa.avail<-FALSE
  }
  
  
  for (x in 1:length(levels)){
    for (y in 1:reps){
      samp<-sample(1:length(match.object$all.dist), size=(100-levels[x]),replace = FALSE)
      
      if (taxa.avail==TRUE){
        match.samp<-site.match(Test=match.object$env.data[nrow(match.object$env.data),],
                               Reference=match.object$env.data[samp,],
                               k= if (!is.null(match.object$k)){match.object$k} else {NULL},
                               adaptive=match.object$adaptive,
                               RDA.reference= if (match.object$method=="RDA-ANNA") {ref.taxa.data[rownames(match.object$env.data[samp,]),]} else {NULL})
        taxa.data<-add.met(Test=test.taxa.data,Reference=ref.taxa.data[names(match.samp$final.dist),])
        tsa.samp<-tsa.test(Test=taxa.data[nrow(taxa.data),],
                           Reference=taxa.data[names(match.samp$final.dist),],
                           outlier.rem=tsa.object$outlier.rem,
                           m.select=tsa.object$m.select,
                           distance= if (!is.null(tsa.object$ecological.distance)){match.samp$final.distance} else {NULL})
      }
      
      if (taxa.avail==FALSE){
        match.samp<-site.match(Test=match.object$env.data[nrow(match.object$env.data),],
                               Reference=match.object$env.data[samp,],
                               k= if (!is.null(match.object$k)){match.object$k} else NULL,
                               adaptive=match.object$adaptive,
                               RDA.reference= if (match.object$method=="RDA-ANNA") {Reference[rownames(match.object$env.data[samp,]),]} else {NULL})
        tsa.samp<-tsa.test(Test=tsa.object$raw.data[nrow(tsa.object$raw.data),colnames(tsa.object$raw.data)%in%colnames(Reference)],
                           Reference=Reference[names(match.samp$final.dist),],
                           outlier.rem=tsa.object$outlier.rem,
                           m.select=tsa.object$m.select,
                           distance= if (!is.null(tsa.object$ecological.distance)){match.samp$final.distance} else {NULL})
      }
      
      full.results[(full.results[,"Level"]==levels[x] & full.results[,"Rep"]==y),3:10]<-c(tsa.samp$tsa.results[1,],tsa.samp$tsa.results[5,],
                                                                                 tsa.samp$tsa.results[9,],tsa.samp$tsa.results[4,],
                                                                                 tsa.samp$general.results[6,],tsa.samp$general.results[5,],
                                                                                 tsa.samp$general.results[3,],tsa.samp$general.results[4,])
    }
    temp.results<-data.frame(full.results[full.results[,"Level"]==levels[x],])
    summary.results[x,c(1,4:6)]<-c(length(which(temp.results[,3]==tsa.object$tsa.results[1,]))/reps,
                                   mean(as.numeric(temp.results[,4])-as.numeric(tsa.object$tsa.results[5,])),
                                   mean(as.numeric(temp.results[,5])-as.numeric(tsa.object$tsa.results[9,])),
                                   mean(as.numeric(temp.results[,6])-as.numeric(tsa.object$tsa.results[4,])))
                                   
    if (tsa.object$tsa.results[1,]=="Not Impaired"){
      summary.results[x,2]<-length(which(temp.results[,3]!="Not Impaired"))/reps
      summary.results[x,3]<-0
    }
    if  (tsa.object$tsa.results[1,]=="Possibly Impaired"){
      summary.results[x,2]<-length(which(temp.results[,3]=="Impaired"))/reps
      summary.results[x,3]<-length(which(temp.results[,3]=="Not Impaired"))/reps
    }
    if (tsa.object$tsa.results[1,]=="Impaired") {
      summary.results[x,2]<-0
      summary.results[x,3]<-length(which(temp.results[,3]!="Impaired"))/reps
    }
  }
  
  
  
  output<-NULL
  output$summary.results<-summary.results
  class(output)<-"tsa.conf"
  return(output)
}

#' @export
grep.paste<-function (data) {#paste0 the vector seperated by |
  substr(paste0(data,collapse="",sep="|"),1,(nchar(paste0(data,collapse="",sep="|"))-1))
}

#' @export
is.nan.data.frame <- function(x) #is.nan function applied to data.frames
  do.call(cbind, lapply(x, is.nan))

#' @export
is.inf.data.frame <- function(x) #is.nan function applied to data.frames
  do.call(cbind, lapply(x, is.infinite))

#' @export
my.sum <- function(x) ifelse( !all(is.na(x)), sum(x, na.rm=T), NA)

#' @export
range01 <- function(x){(x - min(x)) / (max(x) - min(x))}

#' @export
"%w/o%" <- function(x, y) x[!x %in% y] #--  x without y

#' @export
adapt.sum<-function (data){
  if (!is.vector(data)) {
    rowSums(data)
  } else {data}
}

#' @export
adapt.sum1<-function (data){
  if (length(data)!=1) {
    sum(data)
  } else {data}
}

#' @export
braydist<-function(x) {
  y<-nrow(x)
  dist<-data.frame(matrix(nrow=y,ncol=(y-1)))
  rownames(dist)<-rownames(x)
  for (i in 1:y){
    dist[i,]<-as.matrix(vegdist(rbind(x[rownames(x)[-c(i)],],x[i,]),"bray"))[y,1:(y-1)]
  }
  return(dist)
}

#' @export
richness.calc<-function(x){
  if (length(x)==0){
    return(rep(0,nrow(x)))
  } else {
    rich<-NULL
    for (i in 1:nrow(x)){
      tally<-NULL
      taxa.present<-colnames(x)[which(x[i,]>0)]
      tally<-length(grep("idae",taxa.present))
      taxa.remain<-taxa.present[-c(grep("idae",taxa.present))]
      if (length(grep(grep.paste(c("Chironominae","Tanypodinae","Diamesinae","Orthocladiinae","Podonominae","Prodiamesinae","Telmatogetoninae")),taxa.remain))>0){
        tally<-tally+length(grep(grep.paste(c("Chironominae","Tanypodinae")),taxa.remain))
        if (length(grep(grep.paste(c("Diamesinae","Orthocladiinae","Podonominae","Prodiamesinae","Telmatogetoninae")),taxa.remain))){
          tally<-tally+1
        }
        taxa.remain<-taxa.remain[-c(grep(grep.paste(c("Chironominae","Tanypodinae","Diamesinae","Orthocladiinae","Podonominae","Prodiamesinae","Telmatogetoninae")),taxa.remain))]
      }
      if (length(taxa.remain)>0){
        taxa.orig<-taxa.present[which(!taxa.present%in%taxa.remain)]
        t1<-unlist(lapply(taxa.orig, function(x) substr(x, start=1,stop=(gregexpr(pattern =';',x)[[1]][1]-1))))
        t2<-unlist(lapply(taxa.remain, function(x) substr(x, start=1,stop=(gregexpr(pattern =';',x)[[1]][1]-1))))
        t3<-t1[which(duplicated(t1))]
        
        tally<-tally+length(which(!t2%in%t3))
      }
      rich[i]<-tally
    }
    return(rich)
  }
}

#' @export
number.ftrait<-function(x,y){
  data(trait.feeding,envir = environment())
  output<-NULL
  for (i in 1:nrow(x)){
    trait<-trait.feeding[which(trait.feeding$TAXON%in%toupper(
      unlist(
        lapply(colnames(x)[which(x[i,]>0)], function(x) substr(x, start=(gregexpr(pattern =';',x)[[1]][1]+1),stop=nchar(x)))))),]
    
    trait<-trait[trait$TRAITVAL==y,]
    if (nrow(trait)>0){
      output[i]<-sum(x[i,which(toupper(
        unlist(
          lapply(colnames(x), function(x) substr(x, start=(gregexpr(pattern =';',x)[[1]][1]+1),stop=nchar(x)))))%in%trait$TAXON)])
    } else {
      output[i]<-0
    }
  }
  return(output)
}

#' @export
number.htrait<-function(x,y){
  data(trait.habit,envir = environment())
  output<-NULL
  for (i in 1:nrow(x)){
    trait<-trait.habit[which(trait.habit$TAXON%in%toupper(
      unlist(
        lapply(colnames(x)[which(x[i,]>0)], function(x) substr(x, start=(gregexpr(pattern =';',x)[[1]][1]+1),stop=nchar(x)))))),]
    
    trait<-trait[trait$TRAITVAL==y,]
    if (nrow(trait)>0){
      output[i]<-sum(x[i,which(toupper(
        unlist(
          lapply(colnames(x), function(x) substr(x, start=(gregexpr(pattern =';',x)[[1]][1]+1),stop=nchar(x)))))%in%trait$TAXON)])
    } else {
      output[i]<-0
    }
  }
  return(output)
}

#' @export
panel.hist <- function(x, ...) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}

#' @export
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}


##########################################################################################################

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

