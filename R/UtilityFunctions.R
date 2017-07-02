
#' @export
grep.paste<-function (data) {#paste0 the vector seperated by |
  substr(paste0(data,collapse="",sep="|"),1,(nchar(paste0(data,collapse="",sep="|"))-1))
}
#' @export
is.nan.data.frame <- function(x) {
  do.call(cbind, lapply(x, is.nan))  
}

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

#' @export
richness.ftrait<- function (x, y) {
  data(trait.feeding, envir = environment())
  output <- NULL
  for (i in 1:nrow(x)) {
    trait <- trait.feeding[which(trait.feeding$TAXON %in% 
                                   toupper(unlist(lapply(colnames(x)[which(x[i, ] > 0)], function(x) substr(x, start = (gregexpr(pattern = ";",x)[[1]][1] + 1), stop = nchar(x)))))), ]
    #trait <- trait[trait$TRAITVAL == y, ]
    
    output[i]<-nrow(trait[trait$TRAITVAL == y, ])
    
  }
  return(output)
}

#' @export
richness.htrait<- function (x, y) {
  data(trait.habit, envir = environment())
  output <- NULL
  for (i in 1:nrow(x)) {
    trait <- trait.habit[which(trait.habit$TAXON %in% 
                                   toupper(unlist(lapply(colnames(x)[which(x[i, ] > 0)], function(x) substr(x, start = (gregexpr(pattern = ";",x)[[1]][1] + 1), stop = nchar(x)))))), ]
    #trait <- trait[trait$TRAITVAL == y, ]
    
    output[i]<-nrow(trait[trait$TRAITVAL == y, ])
    
  }
  return(output)
}

#' @export
pcout.ba <-function(x,makeplot=FALSE,
           explvar=0.99,crit.M1=1/3,crit.c1=2.5,crit.M2=1/4,crit.c2=0.99,
           cs=0.25,outbound=0.25, ...){
    
    #################################################################
    p=ncol(x)
    n=nrow(x)
    
    x.mad=apply(x,2,mad)
    if (any(x.mad==0))
      stop("More than 50% equal values in one or more variables!")
    
    #################################################################
    # PHASE 1:
    
    # Step 1: robustly sphere the data:
    x.sc <- scale(x,apply(x,2,median),x.mad)
    
    # Step 2: PC decomposition; compute p*, robustly sphere:
    ##x.wcov <- cov(x.sc)
    ##x.eig <- eigen(x.wcov)
    ##p <- ncol(x)
    ##p1 <- (1:p)[(cumsum(x.eig$val)/sum(x.eig$val)>explvar)][1]
    
    # or faster:
    x.svd <- svd(scale(x.sc,TRUE,FALSE))
    a <- x.svd$d^2/(n-1)
    p1 <- (1:p)[(cumsum(a)/sum(a)>explvar)][1]
    
    x.pc <- x.sc%*%x.svd$v[,1:p1]
    xpc.sc <- scale(x.pc,apply(x.pc,2,median),apply(x.pc,2,mad))
    
    # Step 3: compute robust kurtosis weights, transform to distances:
    wp <- abs(apply(xpc.sc^4,2,mean)-3)
    
    xpcw.sc <- xpc.sc%*%diag(wp/sum(wp))
    xpc.norm <- sqrt(apply(xpcw.sc^2,1,sum))
    x.dist1 <- xpc.norm*sqrt(qchisq(0.5,p1))/median(xpc.norm)
    
    # Step 4: determine weights according to translated biweight:
    M1 <- quantile(x.dist1,crit.M1)
    const1 <- median(x.dist1)+crit.c1*mad(x.dist1)
    w1 <- (1-((x.dist1-M1)/(const1-M1))^2)^2
    w1[x.dist1<M1] <- 1
    w1[x.dist1>const1] <- 0
    
    #################################################################
    # PHASE 2:
    
    # Step 5: compute Euclidean norms of PCs and their distances:
    xpc.norm <- sqrt(apply(xpc.sc^2,1,sum))
    x.dist2 <- xpc.norm*sqrt(qchisq(0.5,p1))/median(xpc.norm)
    
    # Step 6: determine weight according to translated biweight:
    M2 <- sqrt(qchisq(crit.M2,p1))
    const2 <- sqrt(qchisq(crit.c2,p1))
    w2 <- (1-((x.dist2-M2)/(const2-M2))^2)^2
    w2[x.dist2<M2] <- 1
    w2[x.dist2>const2] <- 0
    
    #################################################################
    # Combine PHASE1 and PHASE 2: compute final weights:
    wfinal <- (w1+cs)*(w2+cs)/((1+cs)^2)
    wfinal01 <- round(wfinal+0.5-outbound)
    
    #################################################################
    # Generate plot:
    if (makeplot){
      op <- par(mfrow=c(3,2), mar=c(4,4,2,2))
      on.exit(par(op))
      
      # location outliers
      plot(x.dist1,xlab="Index",ylab="Distance (location)", ...)
      abline(h=const1)
      abline(h=M1,lty=2)
      plot(w1,xlab="Index",ylab="Weight (location)",ylim=c(0,1), ...)
      abline(h=0)
      abline(h=1,lty=2)
      
      # scatter outliers
      plot(x.dist2,xlab="Index",ylab="Distance (scatter)", ...)
      abline(h=const2)
      abline(h=M2,lty=2)
      plot(w2,xlab="Index",ylab="Weight (scatter)",ylim=c(0,1), ...)
      abline(h=0)
      abline(h=1,lty=2)
      
      # combined weights
      plot(wfinal,xlab="Index",ylab="Weight (combined)",ylim=c(0,1), ...)
      abline(h=cs)
      plot(wfinal01,xlab="Index",ylab="Final 0/1 weight",ylim=c(0,1), ...)
    }
    
    list(wfinal01=wfinal01,wfinal=wfinal,wloc=w1,wscat=w2,
         x.dist1=x.dist1,x.dist2=x.dist2,M1=M1,const1=const1,M2=M2,const2=const2)
  }

#' @export
D2.dist<- function (data, cov, inverted = FALSE) {
  if (!inherits(data, c("data.frame", "matrix"))) 
    stop("data must be a data.frame or matrix!")
  stopifnot(is.matrix(cov))
  if (ncol(data) != ncol(cov)) 
    stop("incompatible dimensions!")
  x <- as.matrix(data)
  n <- nrow(x)
  D2 <- matrix(0, n, n)
  dimnames(D2) <- list(rownames(data), rownames(data))
  if (!inverted) {
    for (i in 1:n) {
      for (j in 1:n) {
        if (i > j) 
          D2[i, j] <- crossprod((x[i, ] - x[j, ]), solve(cov, 
                                                         (x[i, ] - x[j, ])))
      }
    }
  }
  else {
    for (i in 1:n) {
      for (j in 1:n) {
        if (i > j) 
          D2[i, j] <- crossprod((x[i, ] - x[j, ]), crossprod(cov, 
                                                             (x[i, ] - x[j, ])))
      }
    }
  }
  return(as.dist(D2))
}
