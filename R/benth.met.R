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

benth.metUI<-function(x,taxa.sep=";",HBI=NULL,CEFI=NULL,f.trait=NULL,h.trait=NULL) {
  
  #MAKE SURE TAXA NAMES KEEP TRAILING ";"
  

  if (any(grepl("1",colnames(x)))) {
    stop("Duplicate taxa names are not permitted")
  }
  
  if (taxa.sep!=";"){
    colnames(x)<-gsub(taxa.sep,";",colnames(x),fixed = T)
  }
  rownames(x)<-gsub(" ","_",rownames(x),fixed = T)

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
  
  attributes<-benth.attributes(x=x,taxa.sep=taxa.sep,HBI=HBI,CEFI=CEFI,f.trait=f.trait,h.trait=h.trait)
  
  taxa.names<-toupper(colnames(x))
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
  taxa.intol<-taxa[,which(attributes$HBI<5)]
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
  
  summ$'HBI'<-apply(taxa,1,function(x) sum((x*attributes$HBI)/sum(x,na.rm=T),na.rm=T))
  
  taxa.rel.cefi<-taxa.rel
  taxa.rel.cefi[taxa.rel.cefi<=0.05]<-0
  summ$'CEFI'<-apply(taxa.rel.cefi,1,function(x) sum((x*attributes$CEFI.V*attributes$CEFI.W)/sum(x*attributes$CEFI.W,na.rm=T),na.rm=T))
  summ$'CEFI'[summ$'CEFI'==0]<-NA
  summ$'Predator.Percent'<-apply(taxa,1,function(x) sum(x[which(attributes$Feeding=="PREDATOR")]))/abund
  summ$'Predator.Richness'<-apply(taxa,1,function(x) length(which(x[which(attributes$Feeding=="PREDATOR")]>0)))
  summ$'ScraperGrazer.Percent'<-apply(taxa,1,function(x) sum(x[which(attributes$Feeding=="SCRAPER/GRAZER")]))/abund
  summ$'ScraperGrazer.Richness'<-apply(taxa,1,function(x) length(which(x[which(attributes$Feeding=="SCRAPER/GRAZER")]>0)))
  summ$'Shredder.Percent'<-apply(taxa,1,function(x) sum(x[which(attributes$Feeding=="SHREDDER")]))/abund
  summ$'Shredder.Richness'<-apply(taxa,1,function(x) length(which(x[which(attributes$Feeding=="SHREDDER")]>0)))
  summ$'Filterer.Percent'<-apply(taxa,1,function(x) sum(x[which(attributes$Feeding=="COLLECTOR-FILTERER")]))/abund
  summ$'Filterer.Richness'<-apply(taxa,1,function(x) length(which(x[which(attributes$Feeding=="COLLECTOR-FILTERER")]>0)))
  summ$'Gatherer.Percent'<-apply(taxa,1,function(x) sum(x[which(attributes$Feeding=="COLLECTOR-GATHERER")]))/abund
  summ$'Gatherer.Richness'<-apply(taxa,1,function(x) length(which(x[which(attributes$Feeding=="COLLECTOR-GATHERER")]>0)))
  summ$'ScraperGrazer.to.Shredder.Collector'<-log(summ$'ScraperGrazer.Percent'/(summ$'Shredder.Percent'+summ$'Gatherer.Percent'))
  summ$'Clinger.Percent'<-apply(taxa,1,function(x) sum(x[which(attributes$Habitat=="CLINGER")]))/abund
  summ$'Clinger.Richness'<-apply(taxa,1,function(x) length(which(x[which(attributes$Habitat=="CLINGER")]>0)))
  summ$'Burrower.Percent'<-apply(taxa,1,function(x) sum(x[which(attributes$Habitat=="BURROWER")]))/abund
  summ$'Burrower.Richness'<-apply(taxa,1,function(x) length(which(x[which(attributes$Habitat=="BURROWER")]>0)))
  summ$'Sprawler.Percent'<-apply(taxa,1,function(x) sum(x[which(attributes$Habitat=="SPRAWLER")]))/abund
  summ$'Sprawler.Richness'<-apply(taxa,1,function(x) length(which(x[which(attributes$Habitat=="SPRAWLER")]>0)))
  summ$'Burrower.to.Sprawler.Clinger'<-log(summ$'Burrower.Percent'/(summ$'Clinger.Percent'+summ$'Sprawler.Percent'))

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
