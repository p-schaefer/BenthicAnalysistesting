#' @export
benth.attributes<- function(x,taxa.sep=";",HBI=NULL,CEFI=NULL,f.trait=NULL,h.trait=NULL) {
  if (is.null(HBI)) {
    #HBI_RAW<-BenthicAnalysistesting::HBI_RAW
    data(HBI_RAW,envir = environment())
  } else {
    HBI_RAW<-data.frame(HBI)
  }
  if (is.null(CEFI)) {
    #CEFI<-BenthicAnalysistesting::CEFI
    data(CEFI,envir = environment())
    CEFI$Taxa<-toupper(CEFI$Taxa)
  } else {
    CEFI<-data.frame(CEFI)
  }
  if (is.null(f.trait)) {
    data(trait.feeding,envir = environment())
    #trait.feeding<-BenthicAnalysistesting::trait.feeding
    trait.feeding$TAXON<-gsub(".",";",trait.feeding$TAXON,fixed=T)
  } else {
    trait.feeding<-data.frame(trait.feeding)
  }
  if (is.null(h.trait)) {
    #trait.habit<-BenthicAnalysistesting::trait.habit
    data(trait.habit,envir = environment())
    trait.habit$TAXON<-gsub(".",";",trait.habit$TAXON,fixed=T)
  } else {
    trait.habit<-data.frame(trait.habit)
  }
  
  if (any(grepl("1",colnames(x)))) {
    stop("Duplicate taxa names are not permitted")
  }
  
  if (taxa.sep!=";"){
    colnames(x)<-gsub(taxa.sep,";",colnames(x),fixed = T)
  }
  rownames(x)<-gsub(" ","_",rownames(x),fixed = T)

  taxa.names<-toupper(colnames(x))
  taxa.levels<-max(lengths(regmatches(colnames(x), gregexpr(";", colnames(x),fixed=T))))+1
  taxa.heirarchy<-matrix(scan(textConnection(taxa.names), sep=";", what=""),ncol=taxa.levels,byrow=T)
  taxa.heirarchy<-cbind(taxa.heirarchy,NA,NA,NA,NA,NA,NA)
  colnames(taxa.heirarchy)[(taxa.levels+1):(taxa.levels+6)]<-c("HBI","FBI","CEFI.V","CEFI.W","Feeding","Habitat")
  
  for (i in 1:nrow(taxa.heirarchy)){
    for (n in taxa.levels:1){
      if (!is.na(taxa.heirarchy[i,taxa.levels+1])|!is.na(taxa.heirarchy[i,taxa.levels+2])){
        next #lower taxonomic level already recognized
      }
      
      q.taxa<-taxa.heirarchy[i,n]
      
      if (q.taxa==""){
        next #an unspecified taxon
      }
      
      HBI_index<-grep(q.taxa,HBI_RAW$SCIENTIFIC_NAME)#Index HBI database for taxa
      
      if (length(HBI_index)>1){
        HBI_index_test<-which(q.taxa==HBI_RAW$SCIENTIFIC_NAME)
        if (!identical(HBI_index_test,integer(0))){
          HBI_index<-HBI_index_test
        }
      }
      
      if (length(HBI_index)>1 & n>1) { #account for multiple matches by including higher level taxon
        q.taxa<-paste0(taxa.heirarchy[i,c(n-1,n)],sep=";",collapse="")
        q.taxa<-substr(q.taxa,start=1,stop=(nchar(q.taxa)-1))
        HBI_index<-which(q.taxa==HBI_RAW$SCIENTIFIC_NAME)
        if (length(HBI_index)>1) { #if still multiple matches, skip to next
          next
        }
      }
      
      if (length(HBI_index)>1){
        if (all(HBI_RAW$FBI[HBI_index]==HBI_RAW$FBI[HBI_index][1])) {
          taxa.heirarchy[i,"FBI"]<-HBI_RAW$FBI[HBI_index][1]
          if (all(HBI_RAW$HBI[HBI_index]==HBI_RAW$HBI[HBI_index][1])) {
            taxa.heirarchy[i,"HBI"]<-HBI_RAW$HBI[HBI_index][1]
          }
          next
        } else {
          next
        }
      }

      if (!identical(HBI_index,integer(0)) & length(HBI_index)==1){
        taxa.heirarchy[i,"HBI"]<-as.numeric(HBI_RAW$HBI[HBI_index])
        taxa.heirarchy[i,"FBI"]<-as.numeric(HBI_RAW$FBI[HBI_index])
        next
      }
    }
  }
  
  for (i in 1:nrow(taxa.heirarchy)){
    for (n in taxa.levels:1){
      if (!is.na(taxa.heirarchy[i,taxa.levels+3])|!is.na(taxa.heirarchy[i,taxa.levels+4])){
        next
      }
      q.taxa<-taxa.heirarchy[i,n]
      if (q.taxa==""){
        next
      }
      CEFI_index<-grep(q.taxa,CEFI$Taxa)
      if (length(CEFI_index)>1){
        CEFI_index_test<-which(q.taxa==CEFI$Taxa)
        if (!identical(CEFI_index_test,integer(0))){
          CEFI_index<-CEFI_index_test
        }
      }

      if (!identical(CEFI_index,integer(0)) & length(CEFI_index)==1){
        taxa.heirarchy[i,"CEFI.V"]<-as.numeric(CEFI$CEFI.V[CEFI_index])
        taxa.heirarchy[i,"CEFI.W"]<-as.numeric(CEFI$CEFI.W[CEFI_index])
        next
      }
    }
  }
  
  for (i in 1:nrow(taxa.heirarchy)){
    for (n in taxa.levels:1){
      if (!is.na(taxa.heirarchy[i,taxa.levels+5])){
        next
      }
      q.taxa<-taxa.heirarchy[i,n]
      if (q.taxa==""){
        next
      }
      feeding_index<-grep(q.taxa,trait.feeding$TAXON)
      if (length(feeding_index)>1){
        feeding_index_test<-which(q.taxa==trait.feeding$TAXON)
        if (!identical(feeding_index_test,integer(0))){
          feeding_index<-feeding_index_test
        }
      }

      if (length(feeding_index)>1 & n>1) { #account for multiple matches by including higher level taxon
        q.taxa<-paste0(taxa.heirarchy[i,c(n-1,n)],sep=";",collapse="")
        q.taxa<-substr(q.taxa,start=1,stop=(nchar(q.taxa)-1))
        feeding_index<-which(q.taxa==trait.feeding$TAXON)
        if (length(feeding_index)>1) { #if still multiple matches, skip to next
          next
        }
      }

      if (!identical(feeding_index,integer(0)) & length(feeding_index)==1){
        taxa.heirarchy[i,"Feeding"]<-as.character(trait.feeding$TRAITVAL[feeding_index])
        next
      }
    }
  }
  
  for (i in 1:nrow(taxa.heirarchy)){
    for (n in taxa.levels:1){
      if (!is.na(taxa.heirarchy[i,taxa.levels+6])){
        next
      }
      q.taxa<-taxa.heirarchy[i,n]
      if (q.taxa==""){
        next
      }
      habitat_index<-grep(q.taxa,trait.habit$TAXON)
      if (length(habitat_index)>1){
        habitat_index_test<-which(q.taxa==trait.habit$TAXON)
        if (!identical(habitat_index_test,integer(0))){
          habitat_index<-habitat_index_test
        }
      }

      if (length(habitat_index)>1 & n>1) { #account for multiple matches by including higher level taxon
        q.taxa<-paste0(taxa.heirarchy[i,c(n-1,n)],sep=";",collapse="")
        q.taxa<-substr(q.taxa,start=1,stop=(nchar(q.taxa)-1))
        habitat_index<-which(q.taxa==trait.habit$TAXON)
        if (length(habitat_index)>1) { #if still multiple matches, skip to next
          next
        }
      }

      if (!identical(habitat_index,integer(0)) & length(habitat_index)==1){
        taxa.heirarchy[i,"Habitat"]<-as.character(trait.habit$TRAITVAL[habitat_index])
        next
      }
    }
  }
  taxa.heirarchy<-data.frame(taxa.heirarchy)
  taxa.heirarchy[,c(1:taxa.levels,(taxa.levels+5),(taxa.levels+6))]<-apply(taxa.heirarchy[,c(1:taxa.levels,(taxa.levels+5),(taxa.levels+6))],2,as.character)
  taxa.heirarchy[,(taxa.levels+1):(taxa.levels+4)]<-apply(taxa.heirarchy[,(taxa.levels+1):(taxa.levels+4)],2,as.numeric)
  return(taxa.heirarchy)
}