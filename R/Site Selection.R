#############################################################################
#Nearest Neighbour Reference Site selection
##############################################################################
#' Match test sites with reference sites
#'
#' Methods for matching test sites with nearest neighbour reference sites based on ecological similarity. The Assessment by Nearest-Neighbour Analysis (ANNA)
#' and Redundancy Analysis ANNA (RDA-ANNA) are implimented as in (\url{http://www.bioone.org/doi/abs/10.1086/678702}). Ecological distance
#' are calculated as Euchlidean distance across significant ordination axis. Axis are standardized and weighted according the proportion
#' of variance explained. Significance of axis is determined by the broken-stick method for ANNA and a permutation test for RDA-ANNA.
#'
#' An option for adaptive selection of k of based on dynamic 1-dimensional k-means clustering implimented in \code{\link[Ckmeans.1d.dp]{Ckmeans.1d.dp}}
#' is available, otherwise k must be supplied.
#'
#' @param Test Vector of environmental variables at Test site
#' @param Reference Data frame of environmental Variables at Reference sites
#' @param k Numeric, the numer of nearest neighbour reference sites to find. Ignored if NULL
#' @param adaptive Logical, whether to use adaptive breaks to select nearest neighbour reference sites.
#' @param ad.factor Power factor used in the inverse-power function if adaptive==T
#' @param ad.constant Constant factor used in the inverse-power function if adaptive==T
#' @param RDA.reference Data frame containing indicator metric data at reference sites. Row names must match with Reference. If it is supplied the function will perform RDA-ANNA, otherwise ANNA.
#' @return $final.dist - Vector containing nearest neighbour reference sites to test site as well as the ecological distance of the test site to each reference site.
#' @return $method - Whether ANNA or RDA-ANNA was used
#' @return $k - The number of reference sites requested or "Adaptive"
#' @return $all.dist - The ecological distance of all reference sites to the test site
#' @return $ref.scores - Ordination axis scores of reference sites
#' @return $test.scores - Ordination axis scores of test site
#' @keywords Test Site Analysis, ANNA
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
#' #Evaluate nearest neighbour selections
#' plot(nn.sites)
#' sitematch.plot(nn.sites)

site.match<-function(Test, Reference, k=NULL, adaptive=T, ad.factor=2, ad.constant=1, RDA.reference=NULL) {
  if (any(colnames(Test)%in%colnames(Reference)==F)|
      (ncol(Test)!=ncol(Reference))) {
    stop("Variable name mismatch between test site and reference set")
  }
  
  if (adaptive==F & is.null(k)){
    k=nrow(Reference)
  }
  Reference<-as.matrix(Reference)
  Reference<-Reference[,apply(Reference,2, min)-apply(Reference,2, max)!=0]
  Reference<-Reference[,!is.na(colnames(Reference))]
  Test<-Test[,colnames(Reference)]
  
  if (is.null(RDA.reference)){
    anna.ref<-prcomp(Reference,center = TRUE, scale = T)
    sig<-length(which(eigenvals(anna.ref)>bstick(anna.ref)))
    anna.test.x<-predict(anna.ref,Test)
    anna.test<-(anna.test.x - apply(anna.ref$x,2,min)) / (apply(anna.ref$x,2,max) - apply(anna.ref$x,2,min))
    anna.ref.x<-decostand(anna.ref$x,method="range")
    sig<-min(sig, ncol(anna.ref.x))
    temp<-rbind(sweep(sweep(anna.ref.x,2,anna.test,"-"),2, eigenvals(anna.ref)/sum(eigenvals(anna.ref)),"*")[,1:sig],rep(0,sig))
    anna.dist<- sort(as.matrix(dist(temp))[,nrow(Reference)+1])
    anna.dist<-anna.dist[-c(which(names(anna.dist)==""))]
  }

  if (!is.null(RDA.reference)){
    if (any((rownames(RDA.reference)==rownames(Reference))==F)) {
      stop("site name mismatch between RDA.reference and Reference")
    }

    if (ncol(RDA.reference)>(ncol(Reference)-1)){
      stop("Too many metrics for number of environmental variables")
      #RDA.reference<-RDA.reference[,colSums(RDA.reference[1:nrow(Reference),])!=0]
      #rda.pca<-prcomp(RDA.reference,center = TRUE, scale = T)
      #sig<-length(which(cumsum(eigenvals(rda.pca)/sum(eigenvals(rda.pca)))<0.95))
      #RDA.reference<-rda.pca$x[,1:sig]

    }
    
    if (ncol(RDA.reference)>(nrow(Reference-1))){
      stop("Too many metrics selected for number of Reference Sites")
    }
    
    if (ncol(RDA.reference)>((1/2)*nrow(Reference-1))){
      warning("Number of indicator metrics greater than 1/2 number of reference sites")
    }
    
    RDA.reference<-scale(RDA.reference,T,T)[,]
    RDA.reference[is.nan(RDA.reference)]<-0
    RDA.reference<-RDA.reference[,colSums(RDA.reference[1:nrow(Reference),])!=0]

    anna.ref<-rda(scale(Reference,T,T)[,],RDA.reference,scale=F)

    sig.pca<-prcomp(Reference,center = T, scale = T)
    sig<-length(which(eigenvals(sig.pca)>bstick(sig.pca)))

    #sig<-length(which(anova(rda(x~RDA.reference,scale=F,data=mydf),by="axis")[,4]<0.05))
    #sig<-length(which(anova(anna.ref,by="axis")[,4]<0.05))
    anna.test.x<-predict(object=anna.ref,newdata=(Test-colMeans(Reference))/apply(Reference,2,sd),model="CCA",type="wa")
    anna.test<-(anna.test.x - apply(anna.ref$CCA$wa,2,min)) / (apply(anna.ref$CCA$wa,2,max) - apply(anna.ref$CCA$wa,2,min))
    anna.ref.x<-decostand(anna.ref$CCA$wa,method="range")
    sig<-min(sig, ncol(anna.ref.x))
    temp<-rbind(sweep(sweep(anna.ref.x,2,anna.test,"-"),2, anna.ref$CCA$eig/sum(anna.ref$CCA$eig),"*")[,1:sig],rep(0,sig))
    anna.dist<-sort(as.matrix(dist(temp))[,nrow(Reference)+1])
    anna.dist<-anna.dist[-c(which(names(anna.dist)==""))]
  }

  final.dist<-NULL
  if (adaptive==T) {
    l<-NULL
    for (i in 3:if (is.null(k)) length(anna.dist) else k){
      l[i]<-diff(anna.dist[1:i],lag=1)[i-2]
      if (l[i]>(1/(i^ad.factor))){
        break
      }
    }
    final.dist<-anna.dist[1:length(l)]

    #if (k==F){
    #  final.dist<-anna.dist[1:Ckmeans.1d.dp(anna.dist[1:nrow(Reference)])$size[1]]
    #}
    #else {
    #  final.dist<-anna.dist[1:Ckmeans.1d.dp(anna.dist[1:nrow(Reference)])$size[1]]
    #}
    #if (length(final.dist)>35){
    #  final.dist<-anna.dist[1:Ckmeans.1d.dp(final.dist,2)$size[1]]
    #}
  } else {
      l<-NULL
    for (i in 3:if (is.null(k)) length(anna.dist) else k){
      l[i]<-diff(anna.dist[1:i],lag=1)[i-2]
      if (l[i]>(ad.constant*(1/(i^ad.factor)))){
        break
      }
    }
    if (mean(anna.dist[1:k])>mean(anna.dist[1:length(l)])) {
      warning("Some reference sites may be poorly matched ecologically to test site. Examine plot for confirmation. Consider using adaptive nearest-neighbour breaks.")
    }
    final.dist<-anna.dist[1:k]
  }
  
    if (final.dist[1]>mean(final.dist)){
    warning("Test site may not match reference set ecologically. Examine plot for confirmation. Consider additional reference sites, or additional ecological data.")
  }

  output<-NULL
  output$final.dist<-final.dist
  output$method<- if (is.null(RDA.reference)) {"ANNA"} else {"RDA-ANNA"}
  output$adaptive<- if (adaptive) {TRUE} else {FALSE}
  output$k<- k
  output$sig.axis<-sig
  output$all.dist<-anna.dist
  output$ref.scores<-anna.ref
  output$test.scores<-anna.test.x
  output$env.data<-rbind(Reference,Test)
  class(output)<-"match.object"
  return(output)
}

#' @export
print.match.object<-function(match.object){
  cat(paste0("Method: ",match.object$method, " with ", match.object$sig.axis, " significant axes"),"\n\n")
  if (!is.null(match.object$k) & match.object$adaptive==T) {
    cat(paste0("Adaptive threshold with ",match.object$k," upper maximum nearest neigbours."))
  }
  if (is.null(match.object$k) & match.object$adaptive==T) {
    cat(paste0("Adaptive threshold with no upper maximum\n\n"))
  }
  if (!is.null(match.object$k) & match.object$adaptive==F) {
    cat(paste0("Fixed threshold of ",match.object$k, " nearest neighbours\n\n"))
  }
    cat("Selected nearest neighbour reference sites:\n")
    print(match.object$final.dist)
}

#' @export
summary.match.object<-function(match.object){
  print(match.object)
}
