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
  if (distance.decay==F & is.null(k)){
    stop("Need either distance.decay == TRUE or k>0")
  }
  
  if(!is.null(k)){
    if (identical(k==0,logical(0))){
      k<-NULL
    }
    if (k==0){
      k<-NULL
    }
  }
  
  if (scale){
    Reference.rda<-Reference
    Reference.rda[,!sapply(Reference.rda,is.factor)]<-data.frame(scale(Reference.rda[,!sapply(Reference.rda,is.factor)]))
    scales<-scale(Reference[,!sapply(Reference,is.factor)])
      
    Test.rda<-Test
    Test.rda[,!sapply(Test.rda,is.factor)]<-data.frame(t(apply(Test.rda[,!sapply(Test.rda,is.factor)],1,function(x,center=attr(scales,"scaled:center"),scale=attr(scales,"scaled:scale")) (x-center)/scale )))
  } else {
    Reference.rda<-Reference
    Test.rda<-Test
  }
  
  if (any(sapply(Reference.rda,is.factor))){
    factors<-rbind(Reference.rda,Test.rda)
    fact.matrix<-data.frame(model.matrix(~.,data.frame(factors[,sapply(factors,is.factor)])))[,-c(1)]
    Reference.rda<-Reference.rda[,!sapply(Reference.rda,is.factor)]
    Test.rda<-Test.rda[,!sapply(Test.rda,is.factor)]
    
    Reference.rda<-cbind(Reference.rda,fact.matrix[1:nrow(Reference.rda),])
    Test.rda<-cbind(Test.rda,fact.matrix[(nrow(Reference.rda)+1):nrow(fact.matrix),])
  }

  if (is.null(RDA.reference)){
    anna.ref<-vegan::rda(Reference.rda,scale = F)
    sig<-length(which(anna.ref$CA$eig>vegan::bstick(anna.ref)))
    var.explained<-as.numeric(anna.ref$CA$eig/sum(anna.ref$CA$eig))[1:sig]
    
    anna.ref.points<-data.frame(anna.ref$CA$u)[,1:sig]
    anna.ref.points<-apply(anna.ref.points,2,scales::rescale,to=c(-1,1))
    anna.ref.points<-t(apply(anna.ref.points,1,function(x,y=var.explained){x*y}))
    
    anna.test.points<-data.frame(predict(anna.ref,Test.rda, type="wa",scale=F))[,1:sig]
    anna.test.points<-data.frame(sapply(1:sig,function(i) scales::rescale(anna.test.points[,i],to=c(-1,1),from=range(anna.ref$CA$u[,i]))))
    colnames(anna.test.points)<-colnames(anna.ref.points)
    rownames(anna.test.points)<-rownames(Test.rda)
    #anna.test.points<-data.frame(t(apply(anna.test.points,1,function(x,min.y=apply(anna.ref$CA$u[,1:sig],2,min),max.y=apply(anna.ref$CA$u[,1:sig],2,max)) (x-min.y)/(max.y-min.y))))
    #anna.test.points<-apply(anna.test.points,2,scales::rescale,to=c(-1,1), from=range(anna.ref$CA$u))
    anna.test.points<-t(apply(anna.test.points,1,function(x,y=var.explained){x*y}))
    
    env.scores<-anna.ref$CA$v[,1:sig]
    env.scores<-apply(env.scores,2,scales::rescale,to=c(-1,1))
    env.scores<-t(apply(env.scores,1,function(x,y=var.explained){x*y}))
  }
  
  if (!is.null(RDA.reference)){
    if (any((rownames(RDA.reference)==rownames(Reference))==F)) {
      stop("site name mismatch between RDA.reference and Reference")
    }
    if (ncol(RDA.reference)>(ncol(Reference)-1)){
      stop("Too many metrics for number of environmental variables")
    }
    if (ncol(RDA.reference)>(nrow(Reference)-1)){
      stop("Too many metrics selected for number of Reference Sites")
    }
    if (ncol(RDA.reference)>((1/2)*nrow(Reference)-1)){
      warning("Number of indicator metrics greater than 0.5* number of reference sites")
    }
    RDA.reference<-scale(RDA.reference,T,T)
    RDA.reference[is.na(RDA.reference)]<-0
    RDA.reference<-RDA.reference[,colSums(RDA.reference[1:nrow(Reference),])!=0]
    
    anna.ref<-vegan::rda(Reference.rda~RDA.reference,scale=F)
    
    #sig<-anova(anna.ref,by="axis",permutations=999)
    #sig<-length(which(sig$`Pr(>F)`<=0.05))
    sig<-length(anna.ref$CCA$eig)
    
    var.explained<-as.numeric(anna.ref$CCA$eig/sum(anna.ref$CCA$eig))[1:sig]
    
    anna.ref.points<-data.frame(anna.ref$CCA$wa)[,1:sig]
    anna.ref.points<-apply(anna.ref.points,2,scales::rescale,to=c(-1,1))
    anna.ref.points<-t(apply(anna.ref.points,1,function(x,y=var.explained){x*y}))
    
    anna.test.points<-data.frame(predict(anna.ref,Test.rda,model="CCA", type="wa",scale=scale))[,1:sig]
    anna.test.points<-data.frame(sapply(1:sig,function(i) scales::rescale(anna.test.points[,i],to=c(-1,1),from=range(anna.ref$CCA$wa[,i]))))
    colnames(anna.test.points)<-colnames(anna.ref.points)
    rownames(anna.test.points)<-rownames(Test.rda)
    #anna.test.points<-data.frame(t(apply(anna.test.points,1,function(x,min.y=apply(anna.ref$CCA$wa[,1:sig],2,min),max.y=apply(anna.ref$CCA$wa[,1:sig],2,max)) (x-min.y)/(max.y-min.y))))
    #anna.test.points<-apply(anna.test.points,2,scales::rescale,to=c(-1,1))
    anna.test.points<-t(apply(anna.test.points,1,function(x,y=var.explained){x*y}))
    
    env.scores<-anna.ref$CCA$v[,1:sig]
    env.scores<-apply(env.scores,2,scales::rescale,to=c(-1,1))
    env.scores<-t(apply(env.scores,1,function(x,y=var.explained){x*y}))
    
    met.scores<-anna.ref$CCA$biplot[,1:sig]
    met.scores<-apply(met.scores,2,scales::rescale,to=c(-1,1))
    met.scores<-t(apply(met.scores,1,function(x,y=var.explained){x*y}))
    
  }
  
  dist.all<-as.matrix(dist(rbind(anna.test.points,anna.ref.points)))
  dist<-dist.all[1:nrow(anna.test.points),(nrow(anna.test.points)+1):ncol(dist.all)]
  
  dist.tf<-dist
  dd.number<-data.frame(sites=rownames(Test))
  dd.number$dd.number<-NA
  for (n in 1:nrow(dist)){
    site<-dist[n,]
    site<-site[order(site)]
    if (distance.decay==T){
      ref.TF<-c(T,T,T,T,diff(site[-c(1:3)])<=dd.constant/(1:(length(site)-4))^dd.factor)
      names(ref.TF)[1:4]<-names(site)[1:4]
      ref.TF[which.min(ref.TF):length(ref.TF)]<-F
      if (!is.null(k)) {
        ref.TF[(k+1):length(ref.TF)]<-F
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
  
  dist.tf<-data.frame(apply(dist.tf,2,as.logical))
  rownames(dist.tf)<-rownames(Test)
  colnames(dist.tf)<-rownames(Reference)
  ordination.scores<-data.frame(rbind(anna.ref.points,anna.test.points))
  ordination.scores$Class<-c(rep("Reference",nrow(Reference)),rep("Test",nrow(Test)))

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
  output$ordination.scores<-ordination.scores
  output$env.ordination.scores<-data.frame(env.scores)
  if (!is.null(RDA.reference)){
    output$met.ordination.scores<-met.scores
  }
  
  class(output)<-"match.object"
  return(output)
}

#' @export
print.match.object<-function(match.object){
  cat(paste0("Method: ",match.object$method, " with ", match.object$sig.axis, " significant axes"),"\n\n")
  if (!is.null(match.object$k) & match.object$distance.decay==T) {
    cat(paste0("Adaptive threshold with ",match.object$k," upper maximum nearest neigbours."))
  }
  if (is.null(match.object$k) & match.object$distance.decay==T) {
    cat(paste0("Adaptive threshold with no upper maximum\n\n"))
  }
  if (!is.null(match.object$k) & match.object$distance.decay==F) {
    cat(paste0("Fixed threshold of ",match.object$k, " nearest neighbours\n\n"))
  }
  #cat("Selected nearest neighbour reference sites:\n")
  #print(match.object$dist.tf)
}

#' @export
summary.match.object<-function(match.object){
  print(match.object)
}
