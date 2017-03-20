##########################################################################################################


#' Mahalanobis distance plot
#'
#' Distance plot of Mahalanobis distance of test site as well as reference set, F-distribution, critical values and Jacknife confidence interval.
#' @param tsa.object An object from \code{tsa.test}
#' @keywords Test Site Analysis, Distance Plot
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
#' Evaluate Results
#' plot(tsa.results)


plot.tsa.object<-function(tsa.object){
  par(mar=c(5.1,4.1,4.1,2.1))
  tsa.dist<-tsa.object$mahalanobis.distance
  nInd<-as.numeric(tsa.object$general.results["Number of Metrics",])
  nRef<-as.numeric(tsa.object$general.results["Number of Reference Sites",])
  tsa.lambda<-as.numeric(tsa.object$tsa.results["TSA Lambda",])
  test.site<-tsa.object$general.results["Test Site",]

  d1<-density(tsa.dist[1:(length(tsa.dist)-1)])
  d2<-density(((nInd*(nRef-1))*rf(1000000, df1=nInd, df2=(nRef-nInd), ncp=tsa.lambda))/((nRef-nInd)*nRef))
  plot(d1,main=paste0(test.site),yaxt="n",xlab="Mahalanobis Distance",ylab="",xlim=c(-1,(max(tsa.dist)+3)))
  polygon(d1,col="grey80")
  lines(d2,lty=2,cex=2,col="grey70")
  abline(v=((nInd*(nRef-1))*qf(0.95, df1=nInd, df2=(nRef-nInd), ncp=tsa.lambda, log=FALSE)/((nRef-nInd)*nRef)), lty=2, col='red')
  abline(v=((nInd*(nRef-1))*qf(0.05, df1=nInd, df2=(nRef-nInd), ncp=tsa.lambda, log=FALSE)/((nRef-nInd)*nRef)), lty=2, col='orange')
  points(tsa.dist[length(tsa.dist)],0, pch="*",col='black',cex=2,lwd=2)
  if (any(names(tsa.object)=="jacknife")) {
    segments(x0=tsa.object$jacknife[2,],y0=0.05,x1=tsa.object$jacknife[3,],y1=0.05,col="black",lwd=2)
    text(tsa.object$jacknife[2,],0.05,labels=paste0("Jacknife Consistency ",substr(tsa.object$jacknife[1,],1,3),"%"),pos=3, offset=0.5,cex=0.70,col='black')
  }

  text(tsa.dist[length(tsa.dist)],0, labels="test-site",pos=3, offset=0.5,cex=1,col='black')
}
######################################################################################################

#' z-score Boxplots
#'
#' Boxplots of z-scores used in Test Site Analysis. Metrics used in calculation of TSA are circled, metrics that
#' contribute significantly to an impairment score are marked with an * along the x axis.
#' @param tsa.object An object from tsa.test()
#' @param ... Arguments passed to boxplot()
#' @keywords Test Site Analysis, Boxplot
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
#' Evaluate Results
#' boxplot(tsa.results)

boxplot.tsa.object <- function(tsa.object,...) {
  def.par <- par(no.readonly = TRUE)
  tsa.stand<-tsa.object$z.scores
  nInd<-ncol(tsa.stand)
  nRef<-nrow(tsa.stand)-1

  part.tsa<-if (!is.null(tsa.object$partial.tsa)) {tsa.object$partial.tsa} else {NULL}
  all.met<-colnames(tsa.stand)
  sel.met<-unlist(strsplit(substr(tsa.object$general.results["Selected Indicator Metrics",],1,(nchar(tsa.object$general.results["Selected Indicator Metrics",])-2)),split=", "))

  cols<-colorRampPalette(brewer.pal(12, "Paired"))(nInd)
  text<-paste(seq(1:ncol(tsa.stand)),colnames(tsa.stand),sep=".")
  b1<-ceiling(length(text)/3)
  b2<-ceiling(length(text)*2/3)
  
  l<-rbind(c(1,1,1),c(1,1,1),c(2,3,4))
  layout(l)

  #suppressWarnings(split.screen(c(2,1)))
  #split.screen(c(1, 3), screen = 2)
  #screen(1)
  par(mar = c(1.9,0.8,1.2,0.8))
  boxplot(tsa.stand[1:nRef,],col=cols,outline=F,yaxt="n",ylim=c(min(tsa.stand)*1.3,max(tsa.stand)*1.1),names=seq(1:nInd),cex.axis=1.2,main="",...)
  title(main=paste0(rownames(tsa.stand)[(nRef+1)]," Boxplot"),cex=1.5)
  points(seq(1:nInd),tsa.stand[(nRef+1),],col="red",pch=19,cex=1)

  points(which(colnames(tsa.stand)%in%sel.met),tsa.stand[nrow(tsa.stand),sel.met],col="black",pch="O",cex=1.75)#This line circles points that are used in analysis
  if (any(part.tsa$p<0.05)) {
    points(which(colnames(tsa.stand)%in%rownames(part.tsa)[part.tsa$p<0.05]),rep((min(tsa.stand)*1.2),length(rownames(part.tsa)[part.tsa$p<0.05])),col="red",pch="*",cex=2)
  }
  #screen(3)
  #par(mar = c(0,0,0,0))
  plot(1, type="n", axes=F, xlab="", ylab="")
  legend("center",text[1:b1],cex=0.9,fill=cols[1:b1],bty="n",x.intersp=0.7,y.intersp=0.7)
  #screen(4)
  #par(mar = c(0,0,0,0))
  plot(1, type="n", axes=F, xlab="", ylab="")
  legend("center",text[(b1+1):b2],cex=0.9,fill=cols[(b1+1):b2],bty="n",x.intersp=0.7,y.intersp=0.7)
  #screen(5)
  #par(mar = c(0,0,0,0))
  plot(1, type="n", axes=F, xlab="", ylab="")
  legend("center",text[(b2+1):length(text)],cex=0.9,fill=cols[(b2+1):length(text)],bty="n",x.intersp=0.7,y.intersp=0.7)

  #close.screen(all=T)
  par(def.par)
}

###############################################################################

#' PCoA of Mahalanobis Distance
#'
#' PCoA of Mahalanobis Distance.
#' @param tsa.object An object from tsa.test()
#' @param vectors A logical argument indicating whether vectors for indicator metrics should be plotted
#' @param supplemental An optional data frame for plotting supplimental environmental vectors
#' @keywords Test Site Analysis, PCoA
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
#' Evaluate Results
#' pcoa.tsa(tsa.results)

pcoa.tsa<-function(tsa.object,vectors=T,supplemental=NULL){
  par(mar=c(5.1,4.1,4.1,2.1))
  supp<-data.frame(supplemental)
  mets<-tsa.object$raw.data
  mets<-mets[,unlist(strsplit(substr(tsa.object$general.results["Selected Indicator Metrics",],1,(nchar(tsa.object$general.results["Selected Indicator Metrics",])-2)),split=", "))]
  if (any(rownames(supp)!=rownames(mets))) {
    stop("Missmatch in rownames with supplemental vectors")
  }
  nInd<-ncol(mets)
  nRef<-nrow(mets)-1
  refsites<-c(rep(1,nRef),0)

  plot1<-capscale(D2.dist(mets,(cov(mets[1:nRef,])),inverted=F)~1,add=F,sqrt.dist=F)
  fig<-ordiplot(plot1,type="n",main=paste(rownames(mets[max(nrow(mets)),])," PCoA Plot",sep=""),
                xlab=paste("MDS ",substr((eigenvals(plot1)[1]/sum(eigenvals(plot1)))*100,1,4),"%"),
                ylab=paste("MDS ",substr((eigenvals(plot1)[2]/sum(eigenvals(plot1)))*100,1,4),"%"))
  points(fig,what="sites",cex=0.8,select=refsites==1,col="black",pch=19)
  points(fig,what="sites",cex=0.8,select=refsites==0,col="red",pch=19)
  suppressWarnings(ordiellipse(plot1,refsites,kind="sd",conf=0.99,draw="line",col="grey20",lty=5,show.groups=1))
  text(fig,what="sites",select=refsites==1,col="black",cex=0.8,pos=3)
  text(fig,what="sites",select=refsites==0,col="red",cex=0.9,pos=3)
  if (vectors==T) {
    plot(envfit(plot1,mets[,colSums(mets)>0],display="sites",na.rm=F,permutations=0),cex=0.8,col="orange")
  }
  if (!is.null(supplemental)) {
    plot(envfit(plot1,supp[,colSums(supp)>0],display="sites",na.rm=F,permutations=0),cex=0.8,col="red")
  }
}

#################################################################################################

#' Scatterplot of ANNA and RDA-ANNA site matching functions
#'
#' Scatterplot of ANNA and RDA-ANNA site matching functions.
#' @param site.match.object An object from \code{\link[BenthicAnalysis]{"site.match"}}
#' @param axis A two character vector indicating which axis to plot
#' @keywords ANNA, RDA-ANNA
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
#' Evaluate Results
#' sitematch.plot(nn.sites)

sitematch.plot<-function(match.object,axis=c(1,2)) {
  par(mar=c(5.1,4.1,4.1,2.1))
  final.dist<-match.object$final.dist
  anna.dist<-match.object$all.dist
  anna.ref<-match.object$ref.scores
  anna.test.x<-match.object$test.scores
  method<-match.object$method

  if (method=="RDA-ANNA"){
    full.data<-rbind(anna.ref$CCA$wa[,eval(axis)],anna.test.x[eval(axis)])
    plot(anna.ref$CCA$wa[,1],anna.ref$CCA$wa[,2],
         xlim=c((min(full.data[,1])*1.15),(max(full.data[,1])*1.15)),
         ylim=c((min(full.data[,2])*1.15),(max(full.data[,2])*1.15)),
         xlab=paste0("RDA ",paste0(axis[1])," (", paste0(substr(anna.ref$CCA$eig/sum(anna.ref$CCA$eig)*100,1,4)[axis[1]]),"%)"),
         ylab=paste0("RDA ",paste0(axis[2])," (", paste0(substr(anna.ref$CCA$eig/sum(anna.ref$CCA$eig)*100,1,4)[axis[2]]),"%)"),
         main=paste0("Nearest Neighbour Ordination for ",rownames(anna.test.x) ))
    points(anna.ref$CCA$wa[names(final.dist),axis[1]],anna.ref$CCA$wa[names(final.dist),axis[2]],pch=19)
    text(x=anna.ref$CCA$wa[names(final.dist),axis[1]],y=anna.ref$CCA$wa[names(final.dist),axis[2]],
         labels=names(final.dist),
         pos=2,offset=0.5,
         cex=0.8,col="grey40")
    points(anna.test.x[,axis[1]],anna.test.x[,axis[2]],pch=19,col="red")
    text(x=anna.test.x[,axis[1]],y=anna.test.x[,axis[2]],labels=rownames(anna.test.x),pos=2,offset=0.5,
         cex=0.8,col="red")
  }

  if (method=="ANNA"){
    full.data<-rbind(anna.ref$x[,eval(axis)],anna.test.x[eval(axis)])
    plot(anna.ref$x[,1],anna.ref$x[,2],
         xlim=c((min(full.data[,1])*1.15),(max(full.data[,1])*1.15)),
         ylim=c((min(full.data[,2])*1.15),(max(full.data[,2])*1.15)),
         xlab=paste0("RDA ",paste0(axis[1])," (", paste0(substr(eigenvals(anna.ref)/sum(eigenvals(anna.ref))*100,1,4)[axis[1]]),"%)"),
         ylab=paste0("RDA ",paste0(axis[2])," (", paste0(substr(eigenvals(anna.ref)/sum(eigenvals(anna.ref))*100,1,4)[axis[2]]),"%)"),
         main=paste0("Nearest Neighbour Ordination for ",rownames(anna.test.x) ))
    points(anna.ref$x[names(final.dist),axis[1]],anna.ref$x[names(final.dist),axis[2]],pch=19)
    text(x=anna.ref$x[names(final.dist),axis[1]],y=anna.ref$x[names(final.dist),axis[2]],
         labels=names(final.dist),
         pos=2,offset=0.5,
         cex=0.8,col="grey40")
    points(anna.test.x[,axis[1]],anna.test.x[,axis[2]],pch=19,col="red")
    text(x=anna.test.x[,axis[1]],y=anna.test.x[,axis[2]],labels=rownames(anna.test.x),pos=2,offset=0.5,
         cex=0.8,col="red")
  }
}


##########################################################################################################

#' Plot of nearest-neighbour reference site distances to test site
#'
#' Plot of nearest-neighbour reference site distances to test site
#' @param site.match.object An object from \code{\link[BenthicAnalysis]{"site.match"}}
#' @keywords ANNA, RDA-ANNA
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
#' Evaluate Results
#' plot(nn.sites)

plot.match.object<-function(match.object){
  par(mar=c(5.1,4.1,4.1,2.1))
  anna.dist<-match.object$all.dist
  final.dist<-match.object$final.dist
  k<-match.object$k
  adaptive<-match.object$adaptive
  test.site<-rownames(match.object$test.scores)

  plot(anna.dist,xlab="Rank",ylab="Distance", main=paste0("Nearest-Neighbour Distance Plot for ",test.site))
  if (adaptive) {
    abline(v=length(final.dist),lty=2,col="grey40")
  } else {
    abline(v=k,lty=2,col="grey40")
  }
}
