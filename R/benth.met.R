#' Benthic Metric Calculation
#'
#' Calculation of a variety of metrics for determing impairment of Benthic Macroinvertebrate Communities.
#' @param x data.frame with sampling events in rows and taxa in columns. No row or column names should be defined
#' @param taxa.sep Character that separates taxa names
#' @param HBI Custom sensitivity values for HBI calculation. Must follow format of data(HBI1,envir = environment())
#' @return $Summary.Metrics - Calculated indicator metrics
#' @return $Raw.Data - Raw taxon data
#' @return $Taxa.List - Concatenated taxon names
#' @return $Site.List - Concatenated site names
#' @keywords Benthic Metrics
#' @export
#' @examples
#' data(YKBioData,envir = environment())
#' benth.met(YKBioData,2,2)

benth.metUI<-function(x,taxa.sep=";",HBI=NULL) {
  if (is.null(HBI)) {
    data(HBI1,envir = environment())
    CEFI<-HBI[,c(3,6,7)]
    CEFI<-na.omit(CEFI)
    HBI<-HBI[,c(3,5)]
    HBI<-na.omit(HBI)
  } else {
    HBI<-data.frame(HBI)
  }
  
  if (taxa.sep!=";"){
    colnames(x)<-gsub(taxa.sep,";",rownames(x))
  }
  rownames(x)<-gsub(" ",".",rownames(x))

  #x[,1:site.fields]<-apply(x[,1:site.fields],2,as.character)
  
  #if (!any(colnames(x) %in% c("V1","X1"))) {
  #  colnames(x)<-gsub(pattern=".", replace=";" ,colnames(x),fixed=T)
  #  x<-rbind(as.character(colnames(x)),x)
  #}
  

  #if (site.fields>1){
  #  site.names<-apply(as.matrix(x[(tax.fields+1):nrow(x),1:site.fields]),1,FUN=paste0,collapse="",sep="-")# get site names
  #  site.names<-substr(site.names,start=1,stop=nchar(site.names)-1)
  #  site.names<-gsub(" ","",site.names)
  #} else if (site.fields==1){
  #  site.names<-x[(tax.fields+1):nrow(x),1]
  #  site.names<-gsub(" ","",site.names)
  #}
  
  #taxa.names<-apply(as.matrix(x[1:tax.fields,(site.fields+1):ncol(x)]),2,FUN=paste0,collapse="",sep=";")# get taxa names
  #taxa.names<-substr(taxa.names,start=1,stop=nchar(taxa.names)-1)
  
  taxa.names<-colnames(x)
  site.names<-rownames(x)
  
  taxa<-x
  #if (nrow(x)-site.fields==1){
  #  taxa<-t(x.frame(apply(taxa,2,as.numeric)))
  #} else {
  #  taxa<-data.frame(apply(taxa,2,as.numeric))
  #}
  #colnames(taxa)<-taxa.names
  #rownames(taxa)<-site.names

  taxa.pa<-vegan::decostand(taxa,method="pa")
  taxa.rel<-sweep(taxa,rowSums(taxa),MARGIN=1,FUN="/")
  taxa.intol<-taxa[,taxa.names[grep(grep.paste(HBI[which(HBI[,2]<5),1]),taxa.names)]]
  n.taxa<-ncol(taxa)

  #summ<-data.frame(matrix(nrow=nrow(taxa),ncol=39))
  summ<-NULL
  #colnames(summ)<-c("Richness","Simpson","Shannon",
  #                  "Percent Dominance","Percent Oligochaeta",
  #                  "Percent Chironomidae","Percent Isopoda","Percent Amphipoda",

  #                  "Percent Coleoptera", "Coleo as Elmidae", "Trich as Hydropsychidae","Ephem as Baetidae",
  #                 "Intolerants Richness","Percent Intolerants",
#
  #                  "Percent EPT","EPT Richness",
  #                  "Ephem Richness","Percent Ephem",
  #                  "Plec Richness","Percent Plec",
  #                  "Trich Richness","Percent Trich",
  #                  "EPT per EPT and Chir","Percent Non Chir Dip","Percent CIGH","HBI","CEFI",
  #                  "Percent Predator", "Percent Scraper", "Percent Shredder", "Percent Filter","Percent Gatherer","Scraper:Shredder+Collector",
  #                  "Percent Clinger","Percent Burrower","Percent Sprawler","Burrower:Sprawler+Clinger")
  
  abund<-rowSums(taxa)

  summ$Richness<-richness.calc(taxa)
  summ$Simpson<-vegan::diversity(taxa,index="simpson")
  summ$Shannon<-vegan::diversity(taxa,index="shannon")
  summ$'Percent.Dominance'<-(apply(taxa, 1, max))/abund
  summ$'Percent.Oligochaeta'<-adapt.sum(taxa[,grep(grep.paste(c("Oligochaetae","Oligochaeta","Oligochaete","Lumbriculida","Haplotaxida","Clitellata","Naidina")),colnames(taxa))])/abund
  summ$'Percent.Chironomidae'<-adapt.sum(taxa[,grep(grep.paste(c("Chironomidae","Chironominae","Tanypodinae","Diamesinae","Orthocladiinae","Podonominae","Prodiamesinae","Telmatogetoninae")),colnames(taxa))])/abund
  summ$'Percent.Isopoda'<-adapt.sum(taxa[,grep("Isopoda",colnames(taxa))])/abund
  summ$'Percent.Amphipoda'<-adapt.sum(taxa[,grep("Amphipoda",colnames(taxa))])/abund
  summ$'Percent.Coleoptera'<-adapt.sum(taxa[,grep("Coleoptera",colnames(taxa))])/abund
  summ$'Coleo.as.Elmidae'<-adapt.sum(taxa[,grep("Elmidae",colnames(taxa))])/adapt.sum(taxa[,grep("Coleoptera",colnames(taxa))])
  
  summ$'Trich.as.Hydropsychidae'<-adapt.sum(taxa[,grep("Hydropsychidae",colnames(taxa))])/adapt.sum(taxa[,grep("Trichoptera",colnames(taxa))])
  summ$'Ephem.as.Baetidae'<-adapt.sum(taxa[,grep("Baetidae",colnames(taxa))])/adapt.sum(taxa[,grep("Ephemeroptera",colnames(taxa))])

  summ$'Percent.EPT'<-adapt.sum(taxa[,grep(paste0("Ephemeroptera|Plecoptera|Trichoptera"),colnames(taxa))])/abund
  summ$'Percent.mEPT'<-(adapt.sum(taxa[,grep(paste0("Ephemeroptera|Plecoptera|Trichoptera"),colnames(taxa))])-adapt.sum(taxa[,grep(paste0("Baetidae|Hydropsychidae"),colnames(taxa))]))/abund
  summ$'Percent.ICHAEBO'<-(summ$'Percent.Oligochaeta'+summ$'Percent.Chironomidae'+summ$'Percent.Isopoda'+summ$'Percent.Amphipoda')+(adapt.sum(taxa[,grep(paste0("Baetidae|Hydropsychidae"),colnames(taxa))])/abund)
  summ$'EPT.Richness'<-richness.calc(taxa[,grep("Ephemeroptera|Plecoptera|Trichoptera",colnames(taxa))])
  summ$'Ephem.Richness'<-richness.calc(taxa[,grep("Ephemeroptera",colnames(taxa))])
  summ$'Percent.Ephem'<-adapt.sum(taxa[,grep("Ephemeroptera",colnames(taxa))])/abund
  summ$'Plec.Richness'<-richness.calc(taxa[,grep("Plecoptera",colnames(taxa))])
  summ$'Percent.Plec'<-adapt.sum(taxa[,grep("Plecoptera",colnames(taxa))])/abund
  summ$'Trich.Richness'<-richness.calc(taxa[,grep("Trichoptera",colnames(taxa))])
  summ$'Percent.Trich'<-adapt.sum(taxa[,grep("Trichoptera",colnames(taxa))])/abund
  
  summ$'EPT.per.EPT.and.Chir'<-adapt.sum(taxa[,grep(paste0("Ephemeroptera|Plecoptera|Trichoptera"),colnames(taxa))])/
    (adapt.sum(taxa[,grep(paste0("Ephemeroptera|Plecoptera|Trichoptera"),colnames(taxa))])+
              adapt.sum(taxa[,grep(grep.paste(c("Chironomidae","Chironominae","Tanypodinae","Diamesinae","Orthocladiinae","Podonominae","Prodiamesinae","Telmatogetoninae")),colnames(taxa))]))
  
  summ$'Percent.Non.Chir.Dip'<-1-(adapt.sum(taxa[,grep(grep.paste(c("Chironomidae","Chironominae","Tanypodinae","Diamesinae","Orthocladiinae","Podonominae","Prodiamesinae","Telmatogetoninae")),colnames(taxa))])/adapt.sum(taxa[,grep(paste0("Diptera"),colnames(taxa))]))
  summ$'Percent.CIGH'<-adapt.sum(taxa[,grep("Corixidae|Hirudinea|Isopoda|Gastropoda",colnames(taxa))])/abund
  
  summ$'Intolerants.Richness'<-apply(taxa.intol, 1, function(x) length(which(x>0)))
  summ$'Percent.Intolerants'<-adapt.sum(taxa.intol)/abund
  
  
  taxa.hbi<-taxa[,taxa.names[grep(grep.paste(HBI[,1]),taxa.names)]]
  t1<-lapply(colnames(taxa.hbi), function(x) substr(x, start=(gregexpr(pattern =';',x)[[1]][1]+1),stop=nchar(x))) #find matches between taxa.hbi and HBI in HBI
  t3<-match(t1,HBI[,1])
  if (any(is.na(t3))) {
    t2<-lapply(taxa.names[grep(grep.paste(t1[which(is.na(t3))]),taxa.names)], function(x) substr(x, stop=(gregexpr(pattern =';',x)[[1]][1]-1),start=1))
    t3[is.na(t3)]<-match(t2,HBI[,1])
  }
  summ$'HBI'<-apply(taxa.hbi,1,function(x) sum(x*HBI[t3, 2])/sum(x))
  

  taxa.rel.cefi<-taxa.rel
  taxa.rel.cefi[taxa.rel.cefi<=0.05]<-0
  taxa.rel.cefi<-taxa.rel.cefi[,taxa.names[grep(grep.paste(CEFI[,1]),taxa.names)]]
  t1<-lapply(colnames(taxa.rel.cefi), function(x) substr(x, start=(gregexpr(pattern =';',x)[[1]][1]+1),stop=nchar(x))) #find matches between taxa.hbi and HBI in HBI
  t3<-match(t1,CEFI[,1])
  if (any(is.na(t3))) {
    t2<-lapply(taxa.names[grep(grep.paste(t1[which(is.na(t3))]),taxa.names)], function(x) substr(x, stop=(gregexpr(pattern =';',x)[[1]][1]-1),start=1))
    t3[is.na(t3)]<-match(t2,HBI[,1])
  }
  summ$'CEFI'<-apply(taxa.rel.cefi,1,function(x) sum(x*CEFI[t3, 2]*CEFI[t3, 3])/sum(x*CEFI[t3, 3]))
  
  summ$'Percent.Predator'<-number.ftrait(taxa,"PREDATOR")/abund
  summ$'Percent.Scraper'<-(number.ftrait(taxa,"SCRAPER")+number.ftrait(taxa,"SCRAPER/GRAZER"))/abund
  summ$'Percent.Shredder'<-number.ftrait(taxa,"SHREDDER")/abund
  summ$'Percent.Filter'<-number.ftrait(taxa,"COLLECTOR-FILTERER")/abund
  summ$'Percent.Gatherer'<-number.ftrait(taxa,"COLLECTOR-GATHERER")/abund
  #summ[,33]<-(number.ftrait(taxa,"SCRAPER")+number.ftrait(taxa,"SCRAPER/GRAZER"))/(number.ftrait(taxa,"SHREDDER")+number.ftrait(taxa,"COLLECTOR"))
  summ$'Scraper.to.Shredder.Collector'<-log(summ$'Percent.Scraper'/(summ$'Percent.Shredder'+summ$'Percent.Gatherer'))
  
  
  summ$'Percent.Clinger'<-(number.htrait(taxa,"CLINGER"))/abund
  summ$'Percent.Burrower'<-(number.htrait(taxa,"BURROWER"))/abund
  summ$'Percent.Sprawler'<-(number.htrait(taxa,"SPRAWLER"))/abund
  summ$'Burrower.to.Sprawler.Clinger'<-log(summ$'Percent.Burrower'/(summ$'Percent.Clinger'+summ$'Percent.Sprawler'))
  
  summ<-as.data.frame(summ)
  rownames(summ)<-rownames(taxa)

  summ[is.nan.data.frame(summ)]<-NA
  summ[is.inf.data.frame(summ)]<-NA
  summ<-summ[,(apply(summ,2,my.sum)!=0 & !is.na(apply(summ,2,my.sum)))]

  output<-NULL
  output$Summary.Metrics<-summ
  output$Raw.Data<-taxa
  output$Taxa.List<-as.character(taxa.names)
  output$Site.List<-as.character(site.names)
  class(output)<-"benth.metric"
  return(output)
}

#' @export
print.benth.metric<-function(benth.metric) {
  print(benth.metric$Summary.Metrics)
}

#' @export
summary.benth.metric<-function(benth.metric) {
  summary(benth.metric$Summary.Metrics)
}
