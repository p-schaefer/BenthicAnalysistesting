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
#' @param distance.decay Logical, whether to use adaptive breaks to select nearest neighbour reference sites.
#' @param dd.factor Power factor used in the inverse-power function if adaptive==T
#' @param dd.constant Constant factor used in the inverse-power function if adaptive==T
#' @param scale Scale ordination
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

site.matchUI<-function(Test, Reference, k=NULL, distance.decay=T, dd.factor=2, dd.constant=1, RDA.reference=NULL, scale=T) {
  if (distance.decay==F & (is.null(k) | k==0)){
    stop("Need either distance.decay == TRUE or k>0")
  }

  if (is.null(RDA.reference)){
    anna.ref<-vegan::rda(Reference,scale = scale)
    sig<-length(which(anna.ref$CA$eig>vegan::bstick(anna.ref)))
    var.explained<-as.numeric(anna.ref$CA$eig/sum(anna.ref$CA$eig))[1:sig]
    anna.ref.points<-data.frame(anna.ref$CA$u)[,1:sig]
    anna.ref.points<-data.frame(vegan::decostand(anna.ref.points,method="range"))
    anna.ref.points<-t(apply(anna.ref.points,1,function(x,y=var.explained){x*y}))
    
    anna.test.points<-data.frame(predict(anna.ref,Test, type="wa",scale=scale))[,1:sig]
    anna.test.points<-(anna.test.points - apply(anna.ref$CA$u,2,min)) / (apply(anna.ref$CA$u,2,max) - apply(anna.ref$CA$u,2,min))
    anna.test.points<-t(apply(anna.test.points,1,function(x,y=var.explained){x*y}))
    
    dist<-as.matrix(dist(rbind(anna.test.points,anna.ref.points)))
    dist<-dist[1:nrow(anna.test.points),(nrow(anna.test.points)+1):ncol(dist)]
  }
  
  if (!is.null(RDA.reference)){
    if (any((rownames(RDA.reference)==rownames(Reference))==F)) {
      stop("site name mismatch between RDA.reference and Reference")
    }
    if (ncol(RDA.reference)>(ncol(Reference)-1)){
      stop("Too many metrics for number of environmental variables")
    }
    if (ncol(RDA.reference)>(nrow(Reference-1))){
      stop("Too many metrics selected for number of Reference Sites")
    }
    if (ncol(RDA.reference)>((1/2)*nrow(Reference-1))){
      warning("Number of indicator metrics greater than 0.5* number of reference sites")
    }
    RDA.reference<-scale(RDA.reference,T,T)
    RDA.reference[is.na(RDA.reference)]<-0
    RDA.reference<-RDA.reference[,colSums(RDA.reference[1:nrow(Reference),])!=0]
    
    if (scale){
      Reference.rda<-scale(Reference)
    } else {
      Reference.rda<-Reference
    }
    anna.ref<-vegan::rda(Reference.rda~RDA.reference,scale=F)
    
    sig<-vegan::anova.cca(anna.ref,by="axis",permutations=999)
    sig<-length(which(sig$`Pr(>F)`<=0.05))
    
    var.explained<-as.numeric(anna.ref$CCA$eig/sum(anna.ref$CCA$eig))[1:sig]
    
    anna.ref.points<-data.frame(anna.ref$CCA$wa)[,1:sig]
    anna.ref.points<-data.frame(vegan::decostand(anna.ref.points,method="range"))
    anna.ref.points<-t(apply(anna.ref.points,1,function(x,y=var.explained){x*y}))
    
    anna.test.points<-data.frame(predict(anna.ref,Test,model="CCA", type="wa",scale=scale))[,1:sig]
    anna.test.points<-(anna.test.points - apply(anna.ref$CCA$wa,2,min)) / (apply(anna.ref$CCA$wa,2,max) - apply(anna.ref$CCA$wa,2,min))
    anna.test.points<-t(apply(anna.test.points,1,function(x,y=var.explained){x*y}))
    
    dist<-as.matrix(dist(rbind(anna.test.points,anna.ref.points)))
    dist<-dist[1:nrow(anna.test.points),(nrow(anna.test.points)+1):ncol(dist)]
  }
  

  dist.tf<-dist
  dd.number<-data.frame(sites=rownames(Test))
  dd.number$dd.number<-NA
  for (n in 1:nrow(dist)){
    site<-dist[n,]
    site<-site[order(site)]
    if (distance.decay==T){
      ref.TF<-c(T,T,T,diff(site[-c(1:3)])<=dd.constant/(1:(length(site)-4))^dd.factor)
      names(ref.TF)[1:3]<-names(site)[1:3]
      ref.TF[which.min(ref.TF):length(ref.TF)]<-F
      if (!is.null(k)) {
        ref.TF[(K+1):length(ref.TF)]<-F
      }
    } else {
      ref.TF<-site
      if (!is.null(k)){
        ref.TF[1:k]<-T
        ref.TF[(k+1):length(ref.TF)]<-F
      }
    }
    dist.tf[n,]<-ref.TF[colnames(dist.tf)]
    dd.number$dd.number[n]<-which.min(ref.TF)
  }
  dist.tf<-as.logical(dist.tf)
  

  output<-NULL
  output$distance.matrix<-dist
  output$TF.matrix<-dist.tf
  
  output$method<- if (is.null(RDA.reference)) {"ANNA"} else {"RDA-ANNA"}
  output$distance.decay<- distance.decay
  output$k<- k
  output$sig.axis<-sig
  output$var.explained<-var.explained
  output$dd.number<-dd.number
  output$ordination<-anna.ref
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
  print(match.object$dist.tf)
}

#' @export
summary.match.object<-function(match.object){
  print(match.object)
}
