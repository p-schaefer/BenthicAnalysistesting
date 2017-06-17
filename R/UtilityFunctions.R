
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
