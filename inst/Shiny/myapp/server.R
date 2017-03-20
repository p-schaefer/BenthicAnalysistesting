#left off working on combining add.met and met.sel in batch mode


shinyServer(function(input, output, session) {
  #########################################################
  #DATA INPUT
  ########################################################
  
  #########################################################
  #Input biological Data
  ########################################################
  
  #bio.data<-reactiveValues()
  bio.data<- reactive({
    inbioFile <- input$inbioFile
    if (is.null(inbioFile)){
      return(NULL)
    } else {
      if (input$metdata==F) {
        d<-benth.met(x=read.csv(inbioFile$datapath, header=F,strip.white=TRUE), tax.fields=input$taxa.names, site.fields=input$site.names, HBI = NULL)
        d$Raw.Data1<-cbind(d$Site.List,d$Raw.Data)
        colnames(d$Raw.Data1)[1]<-"Sites"
        d$Summary.Metrics1<-cbind(d$Site.List,d$Summary.Metrics)
        colnames(d$Summary.Metrics1)[1]<-"Sites"
        d$untransformed.metrics<-d$Summary.Metrics1
        d$transformations<-matrix(ncol=2,nrow=ncol(d$Summary.Metrics))
        d$transformations[1:ncol(d$Summary.Metrics),2]<-rep("None",ncol(d$Summary.Metrics))
        d$transformations[1:ncol(d$Summary.Metrics),1]<-colnames(d$Summary.Metrics)
        colnames(d$transformations)<-c("Metric","Trans")
        d
        
      } else {
        x<-read.csv(inbioFile$datapath, header=T,strip.white=TRUE)
        site.fields<-input$site.names
        if (site.fields>1){
          site.names<-apply(as.matrix(x[,1:site.fields]),1,FUN=paste0,collapse="",sep="-")# get site names
          site.names<-substr(site.names,start=1,stop=nchar(site.names)-1)
          taxa<-data.frame(x[,(site.fields+1):ncol(x)])
          rownames(taxa)<-site.names
          d<-NULL
          d$Summary.Metrics<-taxa
          d$Raw.Data<-NULL
          d$Taxa.List<-NULL
          d$Site.List<-site.names
          d$Raw.Data1<-NULL
          d$Summary.Metrics1<-cbind(d$Site.List,d$Summary.Metrics)
          colnames(d$Summary.Metrics1)[1]<-"Sites"
          d$untransformed.metrics<-d$Summary.Metrics
          d$transformations<-matrix(ncol=2,nrow=ncol(d$Summary.Metrics))
          d$transformations[1:ncol(d$Summary.Metrics),2]<-rep("None",ncol(d$Summary.Metrics))
          d$transformations[1:ncol(d$Summary.Metrics),1]<-colnames(d$Summary.Metrics)
          colnames(d$transformations)<-c("Metric","Trans")
          class(d)<-"benth.metric"
          d
          
        } else if (site.fields==1){
          site.names<-x[1:nrow(x),1]
          taxa<-data.frame(x)
          rownames(taxa)<-site.names
          d<-NULL
          d$Summary.Metrics<-taxa
          d$Raw.Data<-NULL
          d$Taxa.List<-NULL
          d$Site.List<-site.names
          d$Raw.Data1<-NULL
          d$Summary.Metrics1<-d$Summary.Metrics
          colnames(d$Summary.Metrics1)[1]<-"Sites"
          d$untransformed.metrics<-d$Summary.Metrics
          d$transformations<-matrix(ncol=2,nrow=ncol(d$Summary.Metrics))
          d$transformations[1:ncol(d$Summary.Metrics),2]<-rep("None",ncol(d$Summary.Metrics))
          d$transformations[1:ncol(d$Summary.Metrics),1]<-colnames(d$Summary.Metrics)
          colnames(d$transformations)<-c("Metric","Trans")
          class(d)<-"benth.metric"
          d
          
        }
        #if (nrow(x)-site.fields==1){
        #  taxa<-t(data.frame(apply(taxa,2,as.numeric)))
        #} else {
        #  taxa<-data.frame(apply(taxa,2,as.numeric))
        #}
      }
      #observeEvent(input$apply.trans,{
      #  if (bio.data.t$modified=="YES"){
      #    d<-reactiveValuesToList(bio.data.t)
      #    #class(d)<-"benth.metric"
      #    d
      #  }
      #})
      #output$testing<-renderTable({d$Summary.Metrics})
      }
  })
  

  output$bio.data.view <- renderDataTable({
    validate(
      need(!is.null(bio.data()),"")
    )
    bio.data()$Raw.Data1
  })
  
  output$metric.data.view <- renderDataTable({
    bio.data()$Summary.Metrics1
  })
  
  output$metric.summary.view <- renderPrint({
    summary(bio.data()$Summary.Metrics)
  })
  
  observeEvent(input$downloadmetricData,{
    if (!is.null(bio.data.t)){
      dir<-choose.dir()
      write.csv(bio.data.t$untransformed.metrics,file=paste0(dir,"/Metrics.csv"))
    }
    
  })
  
  observeEvent(input$downloadtransmetricData,{
    if (!is.null(bio.data.t)){
      dir<-choose.dir()
      write.csv(bio.data.t$Summary.Metrics,file=paste0(dir,"/Transformed Metrics.csv"))
    }
  })
  
  
  output$downloadbioData <- downloadHandler(
    filename = function() { paste(input$inbioFile,'.csv', sep='') },
    content = function(file) {
      write.csv(bio.data()$Raw.Data, file)
    }
  )
  
  output$downloadmetricData <- downloadHandler(
    filename = function() { paste(input$inbioFile,'.csv', sep='') },
    content = function(file) {
      write.csv(bio.data()$Summary.Metrics, file)
    }
  )
  
  #######################################################################
  #Transform Biological Data
  #######################################################################
  
  output$sel.met.for.trans<-renderUI({
    if (is.null(bio.data())){
      helpText("Input data required")
    }
    selectInput('met.for.trans',label=h4("Select") ,choices=colnames(bio.data()$Summary.Metrics), multiple=F, selectize=T)
  })
  
  met.for.trans<-reactive({
    met.for.trans<-input$met.for.trans
    if (is.null(bio.data())){
      return(NULL)
    } else{
      met.for.trans
    }
  })
  
  trans.metric<-reactive({
    
    validate(
      need(if (input$trans=="Log10" & any(bio.data()$Summary.Metrics[,met.for.trans()]==0)) {FALSE} else {TRUE}, "Metric contains 0's, try log10(x+1)"),
      need(if (input$trans=="Inverse" & any(bio.data()$Summary.Metrics[,met.for.trans()]==0)) {FALSE} else {TRUE}, "Metric contains 0's"),
      need(if (input$trans=="Arcsine Sqare Root" & (bio.data()$Summary.Metrics[,met.for.trans()]<0 || (bio.data()$Summary.Metrics[,met.for.trans()]>1))) {FALSE} else {TRUE}, "Transofmration only available for values between 0-1"),
      need(input$trans!="Delete","")
    )

    if (input$trans=="None"){
      t.metric<-try(bio.data()$untransformed.metrics[,met.for.trans()],silent = T)
    }
    if (input$trans=="Log10"){
      t.metric<-try(log(bio.data()$untransformed.metrics[,met.for.trans()]),silent = T)
    }
    if (input$trans=="Log10+1" ){
      t.metric<-try(log(bio.data()$untransformed.metrics[,met.for.trans()]+1),silent = T)
    }
    if (input$trans=="Square Root" ){
      t.metric<-try(sqrt(bio.data()$untransformed.metrics[,met.for.trans()]),silent = T)
    }
    if (input$trans=="Inverse" ){
      t.metric<-try(1/(bio.data()$untransformed.metrics[,met.for.trans()]),silent = T)
    }
    if (input$trans=="Arcsine Sqare Root"){
      t.metric<-try(asin(sqrt(bio.data()$untransformed.metrics[,met.for.trans()])),silent = T)
    }
    if (input$trans=="Logit"){
      t.metric<-try(car::logit(bio.data()$untransformed.metrics[,met.for.trans()]),silent = T)
    }
    
    validate(
      need(!(is(t.metric,"try-error")),"")
    )
    if (!is(t.metric,"try-error")) {
      t.metric
    } else {
      NULL
    }
  })
  
  output$met.trans.plot1<-renderPlot({
    validate(
      need(met.for.trans() != "", "Please select a metric"),
      need(input$trans!="Delete", "Metric to be deleted")
    )
    hist(as.numeric(trans.metric()),prob=F,col="grey", main=paste0(input$trans," ",met.for.trans()), xlab="")
    #lines(density(as.numeric(bio.data()$Summary.Metrics[,met.for.trans()])),lwd=2,col="blue")
  })
  
  output$met.trans.plot2<-renderPlot({
    validate(
      need(met.for.trans() != "", "Please select a metric"),
      need(input$trans!="Delete", "Metric to be deleted")
    )
    qqnorm(as.numeric(trans.metric()), main=paste0(input$trans," ",met.for.trans()), xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
    qqline(as.numeric(trans.metric()), datax = FALSE, distribution = qnorm, probs = c(0.25, 0.75), qtype = 7)
  })
  
  

  bio.data.t<-reactiveValues()
  bio.data.t$Summary.Metrics=NULL
  bio.data.t$transformations=NULL
  bio.data.t$modified=NULL
  bio.data.t$Raw.Data<-NULL
  bio.data.t$Taxa.List<-NULL
  bio.data.t$Site.List<-NULL
  bio.data.t$Raw.Data1<-NULL
  bio.data.t$untransformed.metrics<-NULL
  
  observe({
    bio.data.t$Summary.Metrics<-bio.data()$Summary.Metrics
    bio.data.t$transformations<-bio.data()$transformations
    bio.data.t$Summary.Metrics1<-bio.data()$Summary.Metrics1
    bio.data.t$Raw.Data<-bio.data()$Raw.Data
    bio.data.t$Taxa.List<-bio.data()$Taxa.List
    bio.data.t$Site.List<-bio.data()$Site.List
    bio.data.t$Raw.Data1<-bio.data()$Raw.Data1
    bio.data.t$untransformed.metrics<-bio.data()$untransformed.metrics
  })

  observeEvent(input$apply.trans,{
    #validate(
    #  need(!has_warning(trans.metric()),"")
    #)
    if(is.null(bio.data.t$modified)){
      bio.data.t$Summary.Metrics<-bio.data()$Summary.Metrics
      bio.data.t$transformations<-bio.data()$transformations
      bio.data.t$Summary.Metrics1<-bio.data()$Summary.Metrics1
      bio.data.t$Raw.Data<-bio.data()$Raw.Data
      bio.data.t$Taxa.List<-bio.data()$Taxa.List
      bio.data.t$Site.List<-bio.data()$Site.List
      bio.data.t$Raw.Data1<-bio.data()$Raw.Data1
      bio.data.t$untransformed.metrics<-bio.data()$untransformed.metrics
      #class(bio.data.t)<-"benth.metric"
      
      if (input$trans=="Delete"){
        bio.data.t$transformations[which(bio.data.t$transformations[,1]%in%met.for.trans()),2]<-paste0(input$trans)
        bio.data.t$Summary.Metrics<-bio.data.t$Summary.Metrics[,-c(which(colnames(bio.data.t$Summary.Metrics)%in%met.for.trans()))]
        bio.data.t$Summary.Metrics1<-bio.data.t$Summary.Metrics1[,-c(which(colnames(bio.data.t$Summary.Metrics1)%in%met.for.trans()))]
      } else {
        bio.data.t$transformations[which(bio.data.t$transformations[,1]%in%met.for.trans()),2]<-paste0(input$trans)
        bio.data.t$Summary.Metrics[,which(colnames(bio.data.t$Summary.Metrics)%in%met.for.trans())]<-trans.metric()
        bio.data.t$Summary.Metrics1[,which(colnames(bio.data.t$Summary.Metrics1)%in%met.for.trans())]<-trans.metric()
        colnames(bio.data.t$Summary.Metrics)[which(colnames(bio.data.t$Summary.Metrics)%in%met.for.trans())]<-paste0(input$trans," ",met.for.trans())
        colnames(bio.data.t$Summary.Metrics1)[which(colnames(bio.data.t$Summary.Metrics1)%in%met.for.trans())]<-paste0(input$trans," ",met.for.trans())
        
      }
      
    } else {
      if (input$trans=="Delete"){
        bio.data.t$transformations[which(bio.data.t$transformations[,1]%in%met.for.trans()),2]<-paste0(input$trans)
        bio.data.t$Summary.Metrics<-bio.data.t$Summary.Metrics[,-c(which(colnames(bio.data.t$Summary.Metrics)%in%met.for.trans()))]
        bio.data.t$Summary.Metrics1<-bio.data.t$Summary.Metrics1[,-c(which(colnames(bio.data.t$Summary.Metrics1)%in%met.for.trans()))]
      } else {
        bio.data.t$transformations[which(bio.data.t$transformations[,1]%in%met.for.trans()),2]<-paste0(input$trans)
        bio.data.t$Summary.Metrics[,which(colnames(bio.data.t$Summary.Metrics)%in%met.for.trans())]<-trans.metric()
        bio.data.t$Summary.Metrics1[,which(colnames(bio.data.t$Summary.Metrics1)%in%met.for.trans())]<-trans.metric()
        colnames(bio.data.t$Summary.Metrics)[which(colnames(bio.data.t$Summary.Metrics)%in%met.for.trans())]<-paste0(input$trans," ",met.for.trans())
        colnames(bio.data.t$Summary.Metrics1)[which(colnames(bio.data.t$Summary.Metrics1)%in%met.for.trans())]<-paste0(input$trans," ",met.for.trans())
      }
    }
    bio.data.t$modified<-"YES"
  })
  
  output$met.trans.table<-renderTable({
    bio.data.t$transformations
  })
  
  output$trans.summary.stats<-renderPrint({
    validate(
      need(met.for.trans() != "", "Please select a metric"),
      need(input$trans!="Delete", "Metric to be deleted")
    )
    summary(trans.metric())
  })
  
  output$transformed.data<-renderDataTable({
    bio.data.t$Summary.Metrics1
  })
  
  
  #########################################################
  #Input Environmental Data
  ########################################################
  
  env.data<-reactive({
    if (is.null(bio.data.t)){
      return(NULL)
    }
    inenvFile <- input$inenvFile
    if (is.null(inenvFile)){
      return(NULL)
    } else {
      d<-read.csv(inenvFile$datapath, header=T,strip.white=TRUE)

      if (input$site.names>1){
        site.names<-apply(as.matrix(d[,1:input$site.names]),1,FUN=paste0,collapse="",sep="-")# get site names
        site.names<-substr(site.names,start=1,stop=nchar(site.names)-1)
      } else if (input$site.names==1){
        site.names<-d[,1]
      }
      
      if (any(bio.data.t$Site.List%in%site.names==F)) {
        validate(
          need(if(any(bio.data.t$Site.List%in%site.names==F)){FALSE}else{TRUE},"Site name mismatch between biological data and environmental data" )
        )
        stop("Site name mismatch between biological data and environmental data")
      } else {
        d<-d[,-c(1:input$site.names)]
        d<-cbind(site.names,d)
        colnames(d)[1]<-"Sites"
      }
      rownames(d)<-site.names
      d
    }
  })
  
  output$env.data.view <- renderDataTable({
    if (is.null(bio.data.t)){
      return(NULL)
    }
    env.data()
  })
  
  output$env.summary.view <- renderPrint({
    if (is.null(bio.data.t)){
      return(NULL)
    }
    summary(env.data())
  })
  
  output$downloadenvData <- downloadHandler(
    filename = function() { paste(input$inenvFile,'.csv', sep='') },
    content = function(file) {
      write.csv(env.data, file)
    }
  )
  
  #########################################################
  #User matched Reference sites
  ########################################################
  
  user.site.match<-reactive({
    inuser.site.matchFile <- input$inrefmatchFile
    if (is.null(inuser.site.matchFile)){
      return(NULL)
    }
    x<-data.frame(read.csv(inuser.site.matchFile$datapath, header=T,strip.white=TRUE))
    
    if (any(as.vector(x[,1])%in%bio.data.t$Site.List==F)) {
      stop("Site mismatch between biological data and user site matched data")
    } 
    if (any(apply(x,1,function(x) any(duplicated(x[-c(1, which(x==""))]))))) {
      stop("Duplicate Reference sample names detected for one or more test samples")
    } 
    if (any(apply(x,1,function(x) any(duplicated(x[x!=""]))))) {
      stop("A sample cannot be classified as both test and reference")
    } 
    
    x
    
  })
  
  output$usersitematch.table<-renderTable({
    if (is.null(user.site.match())){
      return(NULL)
    }
    user.site.match()
  })
  
  output$usersitematchavail<-reactive({
    if (is.null(user.site.match())){
      return(0)
    } else {
      return(1)
    }
  })
  outputOptions(output, 'usersitematchavail', suspendWhenHidden=FALSE)
  
  
  #########################################################
  #Identify Reference Sites
  ########################################################
  refID.data<-reactive({
    if (is.null(bio.data.t)){
      return(NULL)
    }
    inrefIDFile <- input$inrefIDFile
    if (is.null(inrefIDFile)){
      return(NULL)
    } else {
      d<-read.csv(inrefIDFile$datapath, header=T,strip.white=TRUE)
      site.fields<-input$site.names
      
      if (site.fields>1){
        site.names<-apply(as.matrix(d[,1:site.fields]),1,FUN=paste0,collapse="",sep="-")# get site names
        site.names<-substr(site.names,start=1,stop=nchar(site.names)-1)
        x<-data.frame(cbind(site.names,d[,ncol(d)]))
        colnames(x)<-c("Sites","Reference")
        
      } else if (site.fields==1){
        x<-data.frame(d)
        colnames(x)<-c("Sites","Reference")
      }
      
      if (any(bio.data.t$Site.List%in%x[,1]==F)) {
        validate(
          need(if(any(bio.data.t$Site.List%in%x[,1]==F)){FALSE}else{TRUE},"Site name mismatch between biological data and environmental data")
        )
        stop("Site name mismatch between biological data and environmental data")
      } else {
      }
      x
    }
  })
  
  output$choose_columns <- renderUI({
    if (is.null(bio.data.t)){
      return(NULL)
    }
    
    if (!is.null(user.site.match())) {
      refID <- user.site.match()[,1]
      colnames <- rownames(bio.data.t$Summary.Metrics)
      
      b1<-ceiling(length(colnames)*1/4)
      b2<-ceiling(length(colnames)*2/4)
      b3<-ceiling(length(colnames)*3/4)
      
      c1<-colnames[1:b1]
      c2<-colnames[(b1+1):b2]
      c3<-colnames[(b2+1):b3]
      c4<-colnames[(b3+1):length(colnames)]
      
      splitLayout(checkboxGroupInput("c1", "",choices  = c1,selected = c1[!c1%in%refID]),
                  checkboxGroupInput("c2", "",choices  = c2,selected = c2[!c2%in%refID]),
                  checkboxGroupInput("c3", "",choices  = c3,selected = c3[!c3%in%refID]),
                  checkboxGroupInput("c4", "",choices  = c4,selected = c4[!c4%in%refID]))
    } else {
      if(!is.null(refID.data())){
        refID <- refID.data()[which(refID.data()[,2]==0),1]
        colnames <- rownames(bio.data.t$Summary.Metrics)
        
        b1<-ceiling(length(colnames)*1/4)
        b2<-ceiling(length(colnames)*2/4)
        b3<-ceiling(length(colnames)*3/4)
        
        c1<-colnames[1:b1]
        #s1<-refID[1:b1]
        c2<-colnames[(b1+1):b2]
        #s2<-refID[(b1+1):b2]
        c3<-colnames[(b2+1):b3]
        #s3<-refID[(b2+1):b3]
        c4<-colnames[(b3+1):length(colnames)]
        #s4<-refID[(b3+1):length(colnames)]
        
        splitLayout(checkboxGroupInput("c1", "",choices  = c1,selected = c1[!c1%in%refID]),
                    checkboxGroupInput("c2", "",choices  = c2,selected = c2[!c2%in%refID]),
                    checkboxGroupInput("c3", "",choices  = c3,selected = c3[!c3%in%refID]),
                    checkboxGroupInput("c4", "",choices  = c4,selected = c4[!c4%in%refID]))
      } else {
        colnames <- rownames(bio.data.t$Summary.Metrics)
        b1<-ceiling(length(colnames)*1/4)
        b2<-ceiling(length(colnames)*2/4)
        b3<-ceiling(length(colnames)*3/4)
        
        c1<-colnames[1:b1]
        c2<-colnames[(b1+1):b2]
        c3<-colnames[(b2+1):b3]
        c4<-colnames[(b3+1):length(colnames)]
        
        splitLayout(checkboxGroupInput("c1", "", choices  = c1),
                    checkboxGroupInput("c2", "", choices  = c2),
                    checkboxGroupInput("c3", "", choices  = c3),
                    checkboxGroupInput("c4", "", choices  = c4))
      }
    }
  })
  
  sel.ref<-reactive ({
    if (is.null(bio.data.t)){
      return(NULL)
    }
    if (!is.null(input$c1)|!is.null(input$c2)|!is.null(input$c3)|!is.null(input$c4)){
      c(input$c1,input$c2,input$c3,input$c4)
    }
    #if (!is.null(user.site.match())){
    #  rownames(bio.data.t$Summary.Metrics)[!rownames(bio.data.t$Summary.Metrics)%in%user.site.match()[,1]]
    #}
  })

  output$usersitematchwasmodified<-reactive({
    refsites<-rownames(bio.data.t$Summary.Metrics)[!rownames(bio.data.t$Summary.Metrics)%in%user.site.match()[,1]]
    if ((any(sel.ref()%in%refsites==F) | any(refsites%in%sel.ref()==F)) & !is.null(user.site.match())){
      return(0)
    } else {
      return(1)
    }
  })
  outputOptions(output, 'usersitematchwasmodified', suspendWhenHidden=FALSE)

  output$selrefID <- renderPrint({
    if (is.null(bio.data.t)){
      return(NULL)
    }
    sel.ref()
  })
  
  output$seltestID <- renderPrint({
    if (is.null(bio.data.t)){
      return(NULL)
    }
    rownames(bio.data.t$Summary.Metrics)[which(!rownames(bio.data.t$Summary.Metrics)%in%sel.ref())]
  })
  
  test.site.choices<-reactive({
    if (is.null(bio.data.t)){
      return(NULL)
    }
    rownames(bio.data.t$Summary.Metrics)[which(!rownames(bio.data.t$Summary.Metrics)%in%sel.ref())]
  })
  
  
  #########################################################
  #INDIVIDUAL SITE ANALYSIS
  ########################################################
  
  #########################################################
  #Test Site Selection
  ########################################################
  
  output$sel.test.site<-renderUI({
    if (is.null(bio.data.t)){
      helpText("Input data required")
    }
    selectInput("test.site", label = h4("Select"), 
                  choices = test.site.choices(), 
                  selected = 1)
    
  })
  
  test.site<-reactive({
    test.site<-input$test.site
    if (is.null(bio.data.t)){
      return(NULL)
    } else{
      test.site
    }
  })
  
  #########################################################
  #Reference site matching
  ########################################################
  
  use.user.site.match<-reactive({input$user.ref.sitematch})
  k<-reactive({input$k.sel})
  adaptive<-reactive({input$adaptive})
  nn.method<-reactive({input$nn.method})
  
  output$nnmethodselected<-reactive({
    if (input$nn.method=="RDA-ANNA"){
      return(0)
    } else {
      return(1)
    }
  })
  outputOptions(output, 'nnmethodselected', suspendWhenHidden=FALSE)
  
  nn.sites<-reactive({
    if (is.null(test.site())){
      return(NULL)
    }
    sel.mets<-sel.mets()
    if (is.null(sel.mets) & nn.method()=="RDA-ANNA"){
      validate(
        need(if(is.null(sel.mets) & nn.method()=="RDA-ANNA"){FALSE}else{TRUE},"Must select indicator metrics before RDA-ANNA is possible")
      )
      stop("Must select indicator metrics before RDA-ANNA is possible")
      return(NULL)
    }
    if (!use.user.site.match()) {
      if (!is.null(env.data())) {
        
        nn.sites<-try(site.match(Test=env.data()[which(env.data()[,"Sites"]%in%test.site()),-c(1)],
                                 Reference=env.data()[which(env.data()[,"Sites"]%in%sel.ref()),-c(1)],
                                 k= if (k()!=0 & k()<nrow(env.data()[which(env.data()[,"Sites"]%in%sel.ref()),-c(1)])) k() else NULL,
                                 adaptive=adaptive(),
                                 RDA.reference = if (nn.method()=="RDA-ANNA") {bio.data.t$Summary.Metrics[which(rownames(bio.data.t$Summary.Metrics)%in%sel.ref()),sel.mets]} else {NULL} ),silent=T)
        validate(
          need(!(is(nn.sites,"try-error")),paste0(attr(nn.sites,"condition")$message))
        )
        nn.sites
      } else {
        return(NULL)
      }

    } else {
      if (!is.null(env.data())) {
        nn.sites<-try(site.match(Test=env.data()[which(env.data()[,"Sites"]%in%test.site()),-c(1)],
                             Reference=env.data()[which(env.data()[,"Sites"]%in%sel.ref()),-c(1)],
                             k= if (k()!=0 & k()<nrow(env.data()[which(env.data()[,"Sites"]%in%sel.ref()),-c(1)])) k() else NULL,
                             adaptive=adaptive(), RDA.reference = if (nn.method()=="RDA-ANNA") {bio.data.t$Summary.Metrics[which(rownames(bio.data.t$Summary.Metrics)%in%sel.ref()),sel.mets]} else {NULL}),silent=T)
        validate(
          need(!(is(nn.sites,"try-error")),paste0(attr(nn.sites,"condition")$message))
        )

        nn.sites$final.dist<-NULL
        nn.sites$final.dist<-nn.sites$all.dist[paste(unlist(user.site.match()[which(user.site.match()[,1]%in%test.site()),2:length(which(user.site.match()[which(user.site.match()[,1]%in%test.site()),]!=""))]))]
        nn.sites
        
      } else {
        nn.sites<-NULL
        nn.sites$final.dist<-paste(unlist(user.site.match()[which(user.site.match()[,1]%in%test.site()),2:length(which(user.site.match()[which(user.site.match()[,1]%in%test.site()),]!=""))]))
        names(nn.sites$final.dist)<-nn.sites$final.dist
        nn.sites$method<-NULL
        nn.sites$method<-"User Selected"
        class(nn.sites)<-"match.object"
        nn.sites
        
      }
    }
  })
  
  
  output$site.match.axis<-renderUI({
    if (is.null(nn.sites())){
      return(NULL)
    } 
    if (nn.sites()$method=="User Selected"){
      return(NULL)
    }
    splitLayout(
      radioButtons("site.axis1", "Axis 1",choices  = 1:nn.sites()$sig.axis,selected = 1),
      radioButtons("site.axis2", "Axis 2",choices  = 1:nn.sites()$sig.axis,selected = 2))
  })
  
  nnord.axis<-reactive ({
    if (is.null(nn.sites())){
      return(NULL)
    } 
    if (nn.sites()$method=="User Selected"){
      return(NULL)
    }
    as.numeric(c(input$site.axis1,input$site.axis2))
  })
  
  
  output$nn.ord<-renderPlot({
    if (is.null(nn.sites())){
      return(NULL)
    } 
    
    method<-nn.sites()$method
    
    if (method=="User Selected"){
      return(NULL)
    }
    if (method=="RDA-ANNA"){
      validate(
        need((input$site.axis1 != "" & 
                input$site.axis2 != "" &
                nn.sites() != ""), "Please select an axis")
      )
      
      final.dist<-nn.sites()$final.dist
      anna.dist<-nn.sites()$all.dist
      anna.ref<-nn.sites()$ref.scores
      anna.test.x<-nn.sites()$test.scores
      
      full.data<-rbind(anna.ref$CCA$wa[,c(as.numeric(input$site.axis1),as.numeric(input$site.axis2))],anna.test.x[c(as.numeric(input$site.axis1),as.numeric(input$site.axis2))])
      plot(x=anna.ref$CCA$wa[,as.numeric(input$site.axis1)],y=anna.ref$CCA$wa[,as.numeric(input$site.axis2)],
           xlim=c((min(full.data[,1])*1.15),(max(full.data[,1])*1.15)),
           ylim=c((min(full.data[,2])*1.15),(max(full.data[,2])*1.15)),
           xlab=paste0("RDA ",paste0(as.numeric(input$site.axis1))," (", paste0(substr(anna.ref$CCA$eig/sum(anna.ref$CCA$eig)*100,1,4)[as.numeric(input$site.axis1)]),"%)"),
           ylab=paste0("RDA ",paste0(as.numeric(input$site.axis2))," (", paste0(substr(anna.ref$CCA$eig/sum(anna.ref$CCA$eig)*100,1,4)[as.numeric(input$site.axis2)]),"%)"),
           main=paste0("Nearest Neighbour Ordination for ",rownames(anna.test.x) ))
      points(anna.ref$CCA$wa[names(final.dist),as.numeric(input$site.axis1)],anna.ref$CCA$wa[names(final.dist),as.numeric(input$site.axis2)],pch=19)
      text(x=anna.ref$CCA$wa[names(final.dist),as.numeric(input$site.axis1)],y=anna.ref$CCA$wa[names(final.dist),as.numeric(input$site.axis2)],
           labels=names(final.dist),
           pos=2,offset=0.5,
           cex=0.8,col="grey40")
      points(anna.test.x[,as.numeric(input$site.axis1)],anna.test.x[,as.numeric(input$site.axis2)],pch=19,col="red")
      text(x=anna.test.x[,as.numeric(input$site.axis1)],y=anna.test.x[,as.numeric(input$site.axis2)],labels=rownames(anna.test.x),pos=2,offset=0.5,
           cex=0.8,col="red")
    }
    if (method=="ANNA"){
      validate(
        need((input$site.axis1 != "" & 
                input$site.axis2 != "" &
                nn.sites() != ""), "Please select an axis")
      )
      
      final.dist<-nn.sites()$final.dist
      anna.dist<-nn.sites()$all.dist
      anna.ref<-nn.sites()$ref.scores
      anna.test.x<-nn.sites()$test.scores
      
      full.data<-rbind(anna.ref$x[,c(as.numeric(input$site.axis1),as.numeric(input$site.axis2))],anna.test.x[c(as.numeric(input$site.axis1),as.numeric(input$site.axis2))])
      plot(x=anna.ref$x[,as.numeric(input$site.axis1)],y=anna.ref$x[,as.numeric(input$site.axis2)],
           xlim=c((min(full.data[,1])*1.15),(max(full.data[,1])*1.15)),
           ylim=c((min(full.data[,2])*1.15),(max(full.data[,2])*1.15)),
           xlab=paste0("RDA ",paste0(as.numeric(input$site.axis1))," (", paste0(substr(eigenvals(anna.ref)/sum(eigenvals(anna.ref))*100,1,4)[as.numeric(input$site.axis1)]),"%)"),
           ylab=paste0("RDA ",paste0(as.numeric(input$site.axis2))," (", paste0(substr(eigenvals(anna.ref)/sum(eigenvals(anna.ref))*100,1,4)[as.numeric(input$site.axis2)]),"%)"),
           main=paste0("Nearest Neighbour Ordination for ",rownames(anna.test.x) ))
      points(anna.ref$x[names(final.dist),as.numeric(input$site.axis1)],anna.ref$x[names(final.dist),as.numeric(input$site.axis2)],pch=19)
      text(x=anna.ref$x[names(final.dist),as.numeric(input$site.axis1)],y=anna.ref$x[names(final.dist),as.numeric(input$site.axis2)],
           labels=names(final.dist),
           pos=2,offset=0.5,
           cex=0.8,col="grey40")
      points(anna.test.x[,as.numeric(input$site.axis1)],anna.test.x[,as.numeric(input$site.axis2)],pch=19,col="red")
      text(x=anna.test.x[,as.numeric(input$site.axis1)],y=anna.test.x[,as.numeric(input$site.axis2)],labels=rownames(anna.test.x),pos=2,offset=0.5,
           cex=0.8,col="red")
    }
})
  
  output$nn.dist<-renderPlot({
    if (is.null(nn.sites())|use.user.site.match()==T){
      return(NULL)
    }
    if (nn.sites()$method!="User Selected"){
      anna.dist<-nn.sites()$all.dist
      final.dist<-nn.sites()$final.dist
      adaptive<-nn.sites()$adaptive
      k<-nn.sites()$k
      test.site<-env.data()[which(env.data()[,"Sites"]%in%test.site()),-c(1)]
      
      plot(anna.dist,xlab="Rank",ylab="Distance", main=paste0("Nearest-Neighbour Distance Plot for ",test.site()))
      if (adaptive) {
        abline(v=length(final.dist),lty=2,col="grey40")
      } else {
        abline(v=k,lty=2,col="grey40")
      }
    }
  })
  
  output$nn.table<-renderPrint({
    if (is.null(nn.sites())){
      return(NULL)
    }
    if (nn.sites()$method!="User Selected"){
      print(nn.sites()$ref.scores)
    }
  })
  
  output$nn.table2<-renderPrint({
    if (is.null(nn.sites())){
      return(NULL)
    }
    nn.sites()$final.dist
    })
  
  #########################################################
  #Identify Indicator Metrics
  ########################################################

  observe({
    if (nn.method()!="RDA-ANNA" & is.null(bio.data.t)){
      bio.data.t$Summary.Metrics<-cbind(bio.data.t$Summary.Metrics, 
                                        add.met(Test=bio.data.t$Raw.Data[test.site(),],Reference=bio.data.t$Raw.Data[names(nn.sites()$final.dist),]),original=F)
    }
  })
  
  output$choose_columns1<-renderUI({
    if (is.null(bio.data.t)){
      return(NULL)
    }
    if (input$metdata==F){
      #colnames<-c(colnames(bio.data.t$Summary.Metrics),c("O:E","Bray-Curtis","CA1","CA2"))
      if (nn.method()=="RDA-ANNA"){
        colnames<-colnames(bio.data.t$Summary.Metrics)
      } else {
        colnames<-c(colnames(bio.data.t$Summary.Metrics),c("O:E","Bray-Curtis","CA1","CA2"))
      }
    } else {
      colnames<-colnames(bio.data.t$Summary.Metrics)
    }
    
    checkboxGroupInput("b1", "", choices  = colnames, selected=NULL)#colnames[colnames%in%sel.mets()])
  })
  
  observe({
    if (input$selectallmet>0){
      if (input$metdata==F){
        if (nn.method()=="RDA-ANNA"){
          colnames<-colnames(bio.data.t$Summary.Metrics)
        } else {
          colnames<-c(colnames(bio.data.t$Summary.Metrics),c("O:E","Bray-Curtis","CA1","CA2"))
        }
      } else {
        colnames<-colnames(bio.data.t$Summary.Metrics)
      }
      
      updateCheckboxGroupInput(session=session, inputId="b1", choices=colnames, selected=colnames)
    }
  })
  
  observe({
    if (input$selectnonemet>0){
      if (input$metdata==F){
        if (nn.method()=="RDA-ANNA"){
          colnames<-colnames(bio.data.t$Summary.Metrics)
        } else {
          colnames<-c(colnames(bio.data.t$Summary.Metrics),c("O:E","Bray-Curtis","CA1","CA2"))
        }
      } else {
        colnames<-colnames(bio.data.t$Summary.Metrics)
      }
      
      updateCheckboxGroupInput(session=session, inputId="b1", choices=colnames, selected=NULL)
    }
  })
  
  
  sel.mets<-reactive({
    input$b1
    #c(input$b1,input$b2,input$b3,input$b4)
  })
  
  output$indicator.pairs.plot<-renderPlot({
    validate(
      need((length(input$b1)>2&input$b1!=""),"Select at least 2 indicator metrics")
    )

    pairs(bio.data.t$Summary.Metrics[names(nn.sites()$final.dist),sel.mets()], 
           diag.panel=panel.hist,
           upper.panel=panel.smooth,
           lower.panel=panel.cor) 
  })
  
  #########################################################
  #Test Site Analysis
  ########################################################
  distance<-reactive({input$distance})
  outlier.rem<-reactive({input$outlier.rem})
  m.select<-reactive({input$mselect})
  outbound<-reactive({as.numeric(input$outbound.input)})
  
  tsa.results<-reactive({
    if (is.null(nn.sites())){
      return(NULL)
    }
    #if (input$metdata==T & input$mselect==T) {
    #  validate(
    #    need(if(input$metdata==T & input$mselect==T){FALSE}else{TRUE},"Automated metric selection not available when input data are indicator metrics")
    #  )
    #  stop("Automated metric selection not available when input data are indicator metrics")
    #}
    
    if (input$distance==T & is.null(env.data())){
      validate(
        need(if(input$distance==T & is.null(env.data())){FALSE}else{TRUE},"Weighing by ecological distance requires ecological data")
      )
      stop("Weighing by ecological distance requires ecological data")
    }
    bio.data.t1<-NULL
    bio.data.t1$Summary.Metrics<-bio.data.t$Summary.Metrics[c(names(nn.sites()$final.dist),test.site()),]
    if(nn.method()!="RDA-ANNA" & input$metdata==F){
      bio.data.t1$Summary.Metrics<-cbind(bio.data.t1$Summary.Metrics, add.met(Test=bio.data.t$Raw.Data[test.site(),],
                                              Reference=bio.data.t$Raw.Data[names(nn.sites()$final.dist),]),original=F)
    }
    
    tsa.results<-try(tsa.test(Test=bio.data.t1$Summary.Metrics[test.site(),sel.mets()],
                              Reference=bio.data.t1$Summary.Metrics[names(nn.sites()$final.dist),sel.mets()],
                              distance= if (distance()) nn.sites()$final.dist else NULL,
                              outlier.rem= outlier.rem(),
                              m.select= m.select(),
                              na.cutoff=0.7,outbound=outbound()),silent=T)
    validate(
      need(!(is(tsa.results,"try-error")),paste0(attr(tsa.results,"condition")$message))
    )
    
    tsa.results
    
    #if (input$metdata==F){
    #  tsa.results<-try(tsa.test(Test=bio.data.t$Summary.Metrics[test.site(),sel.mets()],
    #                        Reference=Reference=bio.data.t$Summary.Metrics[names(nn.sites()$final.dist),sel.mets()],
    #                        distance= if (distance()) nn.sites()$final.dist else NULL,
    #                        outlier.rem= outlier.rem(),
    #                        m.select= m.select(),
    #                        na.cutoff=0.7,outbound=outbound()),silent=T)
    #  validate(
    #    need(!(is(tsa.results,"try-error")),paste0(attr(tsa.results,"condition")$message))
    #  )
      
    #  tsa.results
    #} else {
    #  taxa.data<-rbind(bio.data.t$Summary.Metrics[names(nn.sites()$final.dist),],bio.data.t$Summary.Metrics[test.site(),])
    #  tsa.results<-try(tsa.test(Test=taxa.data[(1+length(nn.sites()$final.dist)),sel.mets()],
    #                        Reference=taxa.data[1:length(nn.sites()$final.dist),sel.mets()],
    #                        distance= if (distance()) nn.sites()$final.dist else NULL,
    #                        outlier.rem= outlier.rem(),
    #                        m.select= m.select(),
    #                        na.cutoff=0.7,outbound=outbound()),silent=T)
    #  validate(
    #    need(!(is(tsa.results,"try-error")),paste0(attr(tsa.results,"condition")$message))
    #  )
    #  tsa.results
    #}
  })
  
  output$display.ref.sites<-renderText({
    matrix(tsa.results()$ref.sites)
  })
  
  output$display.outlier.ref.sites<-renderText({
    validate(
      need(outlier.rem(),"")
    )
    matrix(tsa.results()$outlier.ref.sites)
  })
  
  output$tsa.distplot<-renderPlot({
    if (is.null(tsa.results())){
      return(NULL)
    }
    
    tsa.dist<-tsa.results()$mahalanobis.distance
    nInd<-as.numeric(tsa.results()$general.results["Number of Metrics",])
    nRef<-as.numeric(tsa.results()$general.results["Number of Reference Sites",])
    tsa.lambda<-as.numeric(tsa.results()$tsa.results["TSA Lambda",])
    test.site<-tsa.results()$general.results["Test Site",]
    
    d1<-density(tsa.dist[1:(length(tsa.dist)-1)])
    d2<-density(((nInd*(nRef-1))*rf(1000000, df1=nInd, df2=(nRef-nInd), ncp=tsa.lambda))/((nRef-nInd)*nRef))
    plot(d1,main=paste0(test.site),yaxt="n",xlab="Mahalanobis Distance",ylab="",xlim=c(-1,(max(tsa.dist)+3)))
    polygon(d1,col="grey80")
    lines(d2,lty=2,cex=2,col="grey70")
    abline(v=((nInd*(nRef-1))*qf(0.95, df1=nInd, df2=(nRef-nInd), ncp=tsa.lambda, log=FALSE)/((nRef-nInd)*nRef)), lty=2, col='red')
    abline(v=((nInd*(nRef-1))*qf(0.05, df1=nInd, df2=(nRef-nInd), ncp=tsa.lambda, log=FALSE)/((nRef-nInd)*nRef)), lty=2, col='orange')
    points(tsa.dist[length(tsa.dist)],0, pch="*",col='black',cex=2,lwd=2)
    if (any(names(tsa.results())=="jacknife")) {
      segments(x0=tsa.results()$jacknife[2,],y0=0.01,x1=tsa.results()$jacknife[3,],y1=0.01,col="black",lwd=2)
      text(tsa.results()$jacknife[2,],0.01,labels=paste0("Jacknife Consistency ",substr(tsa.results()$jacknife[1,],1,3),"%"),pos=3, offset=0.5,cex=0.85,col='black')
    }
    
    text(tsa.dist[length(tsa.dist)],0, labels="test-site",pos=3, offset=0.5,cex=1,col='black')
  })
  
  output$tsa.boxplot<-renderPlot({
    if (is.null(tsa.results())){
      return(NULL)
    }
    bio.data.t1<-NULL
    bio.data.t1$Summary.Metrics<-bio.data.t$Summary.Metrics[c(names(nn.sites()$final.dist),test.site()),]
    if(nn.method()!="RDA-ANNA" & input$metdata==F){
      bio.data.t1$Summary.Metrics<-cbind(bio.data.t1$Summary.Metrics, add.met(Test=bio.data.t$Raw.Data[test.site(),],
                                                                              Reference=bio.data.t$Raw.Data[names(nn.sites()$final.dist),]),original=F)
    }
    
    tsa.stand<-tsa.zscore(Test=bio.data.t1$Summary.Metrics[test.site(),],Reference=bio.data.t1$Summary.Metrics[names(nn.sites()$final.dist),])
    
    #if (input$metdata==F){
    #  tsa.stand<-tsa.zscore(Test=add.met(Test=bio.data.t$Raw.Data[test.site(),],Reference=bio.data.t$Raw.Data[names(nn.sites()$final.dist),])[(1+length(nn.sites()$final.dist)),],
    #                        Reference=add.met(Test=bio.data.t$Raw.Data[test.site(),],Reference=bio.data.t$Raw.Data[names(nn.sites()$final.dist),])[1:length(nn.sites()$final.dist),])
    #} else {
    #  tsa.stand<-tsa.zscore(Test=bio.data.t$Summary.Metrics[test.site(),],Reference=bio.data.t$Summary.Metrics[names(nn.sites()$final.dist),])
    #}

    nInd<-ncol(tsa.stand)
    nRef<-nrow(tsa.stand)-1
    
    part.tsa<-if (!is.null(tsa.results()$partial.tsa)) {tsa.results()$partial.tsa} else {NULL}
    all.met<-colnames(tsa.stand)
    sel.met<-unlist(strsplit(substr(tsa.results()$general.results["Selected Indicator Metrics",],1,(nchar(tsa.results()$general.results["Selected Indicator Metrics",])-2)),split=", "))
    
    cols<-colorRampPalette(brewer.pal(12, "Paired"))(nInd)
    text<-paste(seq(1:ncol(tsa.stand)),colnames(tsa.stand),sep=".")
    b1<-ceiling(length(text)/3)
    b2<-ceiling(length(text)*2/3)
    
    suppressWarnings(split.screen(c(2,1)))
    split.screen(c(1, 3), screen = 2)
    screen(1)
    par(mar = c(1.9,0.8,1.2,0.8))
    boxplot(tsa.stand[1:nRef,],col=cols,outline=F,yaxt="n",ylim=c(min(tsa.stand)*1.3,max(tsa.stand)*1.1),names=seq(1:nInd),cex.axis=0.6,main="")
    title(main=paste0(rownames(tsa.stand)[(nRef+1)]," Boxplot"),cex.main=0.75)
    points(seq(1:nInd),tsa.stand[(nRef+1),],col="red",pch=19,cex=1)
    
    points(which(colnames(tsa.stand)%in%sel.met),tsa.stand[nrow(tsa.stand),sel.met],col="black",pch="O",cex=1.5)#This line circles points that are used in analysis
    if (any(part.tsa$p<0.05)) {
      points(which(colnames(tsa.stand)%in%rownames(part.tsa)[part.tsa$p<0.05]),rep((min(tsa.stand)*1.2),length(rownames(part.tsa)[part.tsa$p<0.05])),col="red",pch="*",cex=2)
    }
    screen(3)
    par(mar = c(0,0,0,0))
    plot(1, type="n", axes=F, xlab="", ylab="")
    legend("center",text[1:b1],cex=0.9,fill=cols[1:b1],bty="n",x.intersp=0.8,y.intersp=0.8)
    screen(4)
    par(mar = c(0,0,0,0))
    legend("center",text[(b1+1):b2],cex=0.9,fill=cols[(b1+1):b2],bty="n",x.intersp=0.8,y.intersp=0.8)
    screen(5)
    par(mar = c(0,0,0,0))
    plot(1, type="n", axes=F, xlab="", ylab="")
    legend("center",text[(b2+1):length(text)],cex=0.9,fill=cols[(b2+1):length(text)],bty="n",x.intersp=0.8,y.intersp=0.8)
    
    close.screen(all=T)
    
  })
  
  output$tsa.pcoa<-renderPlot({
    #supp<-data.frame(supplemental)
    if (is.null(tsa.results())){
      return(NULL)
    }
    supplemental<-NULL
    supp<-NULL
    vectors=T
    mets<-tsa.results()$raw.data
    mets<-mets[,unlist(strsplit(substr(tsa.results()$general.results["Selected Indicator Metrics",],1,(nchar(tsa.results()$general.results["Selected Indicator Metrics",])-2)),split=", "))]
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
    suppressWarnings(ordiellipse(plot1,refsites,kind="sd",conf=0.95,draw="line",col="grey20",lty=5,show.groups=1))
    text(fig,what="sites",select=refsites==1,col="black",cex=0.8,pos=3)
    text(fig,what="sites",select=refsites==0,col="red",cex=0.9,pos=3)
    if (vectors==T) {
      plot(envfit(plot1,mets[,colSums(mets)>0],display="sites",na.rm=F,permutations=0),cex=0.8,col="orange")
    }
    if (!is.null(supplemental)) {
      plot(envfit(plot1,supp[,colSums(supp)>0],display="sites",na.rm=F,permutations=0),cex=0.8,col="red")
    }
  })

  output$tsa.results<-renderPrint({
    if (is.null(tsa.results())){
      return(NULL)
    }
    tsa.results()$tsa.results
  })
  output$ptsa.results<-renderPrint({
    if (is.null(tsa.results())){
      return(NULL)
    }
    tsa.results()$partial.tsa
  })
  
  output$tsa.jack<-renderPrint({
    if (is.null(tsa.results())){
      return(NULL)
    }
    tsa.results()$jacknife
  })
  
  output$tsa.ca<-renderPlot({
    if (is.null(tsa.results())){
      return(NULL)
    }
    if (input$metdata){
      return(NULL)
    }
    Reference<-bio.data.t$Raw.Data[names(tsa.results()$mahalanobis.distance)[1:(length(tsa.results()$mahalanobis.distance)-1)],]
    nRef<-nrow(Reference)
    Test<-bio.data.t$Raw.Data[names(tsa.results()$mahalanobis.distance)[(length(tsa.results()$mahalanobis.distance))],]
    raw.data<-rbind(Reference,Test)
    pRef<-colSums(decostand(Reference,"pa"))/nrow(Reference)
    
    ca.ord<-cca(log(raw.data[,names(which(pRef>=0.1))]+1))
    #ca1<-c(ca.ord$CA$u[,1],predict(ca.ord,log(Test[,names(which(pRef>=0.1))]+1),type="wa")[1])
    #names(ca1)[nRef+1]<-rownames(Test)
    #ca2<-c(ca.ord$CA$u[,2],predict(ca.ord,log(Test[,names(which(pRef>=0.1))]+1),type="wa")[2])
    #names(ca2)[nRef+1]<-rownames(Test)
    ca1<-ca.ord$CA$u[,1]
    ca2<-ca.ord$CA$u[,2]
    
    plot(ca.ord,type="n",main=paste(rownames(ca1)[(nRef+1)]," CA Plot",sep=""),
         xlab=paste("CA1 ",substr((eigenvals(ca.ord)[1]/sum(eigenvals(ca.ord)))*100,1,4),"%"),
         ylab=paste("CA2 ",substr((eigenvals(ca.ord)[2]/sum(eigenvals(ca.ord)))*100,1,4),"%"))
    text(x=ca.ord$CA$v[,1],y=ca.ord$CA$v[,2],labels=rownames(ca.ord$CA$v),col="grey50",cex=0.7)
    points(x=ca1[1:nRef],y=ca2[1:nRef],cex=0.8,col="black",pch=19)
    points(x=ca1[(nRef+1)],y=ca2[(nRef+1)],cex=0.8,col="red",pch=19)
    text(x=ca1[1:nRef],y=ca2[1:nRef],labels=names(ca1)[1:nRef],col="black",cex=0.8,pos=3)
    text(x=ca1[(nRef+1)],y=ca2[(nRef+1)],labels=names(ca1)[(nRef+1)],col="red",cex=0.9,pos=3)
    suppressWarnings(ordiellipse(ca.ord,c(rep(1,nrow(Reference)),0),kind="sd",conf=0.95,draw="line",col="grey20",lty=5,show.groups=1))
  })
  
  output$print.sel.met<-renderPrint({
    if (is.null(tsa.results())){
      return(NULL)
    }
    if (outlier.rem()){
      Reference<-tsa.results()$raw.data.with.outliers[names(nn.sites()$final.dist),]
      Test<-tsa.results()$raw.data.with.outliers[test.site(),]
    } else {
      Reference<-tsa.results()$raw.data[names(nn.sites()$final.dist),]
      Test<-tsa.results()$raw.data[test.site(),]
    }
    
    sel.met<-try(metric.select(Test=Test,Reference=Reference,outlier.rem = outlier.rem(), rank=F,outbound= if(outlier.rem()){outbound()}else {0}),silent=T)
    
    validate(
      need(!(is(sel.met,"try-error")),paste0(attr(sel.met,"condition")$message))
    )
    
    sel.met
  })
  
  #########################################################
  #BATCH ANALYSIS
  ########################################################
  output$batch.nn.method<-renderUI({
    if (is.null(env.data()) & is.null(user.site.match())){
      return("Enter Data First")
    } else{
      if (!is.null(env.data()) & !is.null(user.site.match())){
        opts<-c("ANNA","RDA-ANNA","User Selected")
        sel<-"ANNA"
      } else {
        if (!is.null(user.site.match())){
          opts<-c("User Selected")
          sel<-"User Selected"
        } 
        if (!is.null(env.data())){
          opts<-c("ANNA","RDA-ANNA")
          sel<-"ANNA"
        }}
      radioButtons("nnmethod",choices=opts,label="Site matching method:",selected=sel,inline=T)
      
    }
  })  
  
  output$envdataavail<-reactive({
    if (!is.null(env.data()) & !is.null(user.site.match())) {
      if (!is.null(input$nnmethod)){
        if (input$nnmethod=="User Selected"){
          return(1)
        }
      }
    } 
  })
  
  output$sitematchdataavail<-reactive({
    if (!is.null(user.site.match())) {
      return(1)
    } 
  })
  output$envdataavail1<-reactive({
    if (!is.null(env.data())) {
      return(1)
    } 
  })
  
  
  output$onlyenv<- reactive({
    if (is.null(user.site.match()) & !is.null(env.data())) {
      return(1)
    } 
  })
  output$onlyuser<- reactive({
    if (!is.null(user.site.match()) & is.null(env.data())) {
      return(1)
    } 
  })
  output$envanduser<- reactive({
    if (!is.null(user.site.match()) & !is.null(env.data())) {
      return(1)
    } 
  })
  output$dataavail<-reactive({
    if (!is.null(env.data()) | !is.null(user.site.match())){
      return(1)
    }
  })
  
  output$ecodistwithuserrefsites<-reactive({
    nnmethod<-input$nnmethod
    ab.distance<-input$ab.distance
    if (is.null(nnmethod) | is.null(ab.distance)){
      return(0)
    } else if (nnmethod=="User Selected" & ab.distance==T){
      return(1)
    } 
  })
  
  outputOptions(output, 'ecodistwithuserrefsites', suspendWhenHidden=FALSE)
  
  outputOptions(output, 'envdataavail', suspendWhenHidden=FALSE)
  outputOptions(output, 'sitematchdataavail', suspendWhenHidden=FALSE)
  outputOptions(output, 'onlyenv', suspendWhenHidden=FALSE)
  outputOptions(output, 'onlyuser', suspendWhenHidden=FALSE)
  outputOptions(output, 'envanduser', suspendWhenHidden=FALSE)
  outputOptions(output, 'dataavail', suspendWhenHidden=FALSE)
  outputOptions(output, 'envdataavail1', suspendWhenHidden=FALSE)
  
  
  output$ab.choose_columns1<-renderUI({
    if (is.null(bio.data.t)){
      return(NULL)
    }
    if (!is.null(input$nnmethod)){
      if (input$metdata==T){
        colnames<-colnames(bio.data.t$Summary.Metrics)[-c(1)]
      } else {
        if (!is.null(env.data())){
          if (input$nnmethod=="ANNA"){
        colnames<-c(colnames(bio.data.t$Summary.Metrics),c("O:E","Bray-Curtis","CA1","CA2"))
          }
          if (input$nnmethod=="RDA-ANNA"){
        colnames<-colnames(bio.data.t$Summary.Metrics)
          }
        }
        if (!is.null(user.site.match())){
          if (input$nnmethod=="User Selected"){
            colnames<-c(colnames(bio.data.t$Summary.Metrics),c("O:E","Bray-Curtis","CA1","CA2"))
          }
        }
      }
      checkboxGroupInput("ab.b1", "", choices  = colnames, selected=NULL)
    }
  })
  
  
  
  observe({
    if (input$ab.selectallmet>0){
      if (input$metdata==T){
        colnames<-colnames(bio.data.t$Summary.Metrics)[-c(1)]
      } else {
        if (!is.null(env.data())){
          if (input$nnmethod=="ANNA"){
            colnames<-c(colnames(bio.data.t$Summary.Metrics),c("O:E","Bray-Curtis","CA1","CA2"))
          }
          if (input$nnmethod=="RDA-ANNA"){
            colnames<-colnames(bio.data.t$Summary.Metrics)
          }
        }
        if (!is.null(user.site.match())){
          if (input$nnmethod=="User Selected"){
            colnames<-c(colnames(bio.data.t$Summary.Metrics),c("O:E","Bray-Curtis","CA1","CA2"))
          }
        }
      }
      updateCheckboxGroupInput(session=session, inputId="ab.b1", choices=colnames, selected=colnames)
    }
  })
  
  observe({
    if (input$ab.selectnonemet>0){
      if (input$metdata==T){
        colnames<-colnames(bio.data.t$Summary.Metrics)[-c(1)]
      } else {
        if (!is.null(env.data())){
          if (input$nnmethod=="ANNA"){
            colnames<-c(colnames(bio.data.t$Summary.Metrics),c("O:E","Bray-Curtis","CA1","CA2"))
          }
          if (input$nnmethod=="RDA-ANNA"){
            colnames<-colnames(bio.data.t$Summary.Metrics)
          }
        }
        if (!is.null(user.site.match())){
          if (input$nnmethod=="User Selected"){
            colnames<-c(colnames(bio.data.t$Summary.Metrics),c("O:E","Bray-Curtis","CA1","CA2"))
          }
        }
      }
      updateCheckboxGroupInput(session=session, inputId="ab.b1", choices=colnames, selected=NULL)
    }
  })
  
  ab.sel.mets<-reactive({
    input$ab.b1
  })
  
  sel.dir<-reactive({
    if (input$ab.dir>0){
      choose.dir()
    } else {
      "No Directory Selected"
    }
    })
  
  output$show.sel.dir<-renderText({sel.dir()})
  
  ab.adaptive<-reactive({
    ab.adaptive1<-input$ab.adaptive1
    ab.adaptive2<-input$ab.adaptive2
    if (ab.adaptive1==F | ab.adaptive2==F){
      return(FALSE)
    } else {
      return(TRUE)
    }
  })
  
  ab.k<-reactive({
    ab.k.sel1<-input$ab.k.sel1
    ab.k.sel2<-input$ab.k.sel2
    if (ab.k.sel1!=0){
      return(ab.k.sel1)
    }
    if (ab.k.sel2!=0){
      return(ab.k.sel2)
    }
    if (is.null(ab.k.sel1) & is.null(ab.k.sel2)){
      return(0)
    }
    if (ab.k.sel1==0 & ab.k.sel2==0){
      return(0)
    }
  })
  
  results<-reactive({
    if (input$ab.go>0){
      if (sel.dir()=="No Directory Selected"){
        return(NULL)
      } else {
        withProgress(message="Working", value=0, {
          results<-data.frame(matrix(nrow=length(test.site.choices()),ncol=16))
          rownames(results)<-test.site.choices()
          colnames(results)<-c("Impairment Rank","Interval Test","Equivalence Test","Randomization Test","Test Site D2","Lower Critical","Upper Critical",
                               "TSA Lambda","TSA F Value","Number of Metrics","Number of Reference Sites","Nearest-Neighbour Method","Jacknife Consistency",
                               "Indicator Metrics","Significant Metrics","Reference Sites")
          new.dirs<-isolate(c(input$ab.nnscatter.plot,input$ab.nndist.plot,input$ab.tsadist.plot,input$ab.tsabox.plot,input$ab.tsascatter.plot,
                      input$ab.cascatter.plot,input$ab.multi.plot))
          new.dir.names<-c("NN Ordination", "NN Distance", "TSA Distance", "TSA Boxplot", "TSA Ordination", "CA Ordination", "Multiplot")
          
          ab.sel.mets<-isolate(ab.sel.mets())
          nnmethod<-isolate(input$nnmethod)
          ab.adaptive<-isolate(ab.adaptive())
          ab.k<-isolate(ab.k())
          ab.outlier.rem<-isolate(input$ab.outlier.rem)
          ab.m.select<-isolate(input$ab.m.select)
          ab.distance<-isolate(input$ab.distance)

          if (any(new.dirs==T)){
            for (i in new.dir.names[new.dirs]){
              dir.create(paste0(sel.dir(),"/",i),showWarnings = F)
            }
          }
          
          for (i in test.site.choices()){
            n<-which(test.site.choices()%in%i)
            incProgress(1/length(test.site.choices()), detail = paste("In progress: ", i))
            

            if (nnmethod=="RDA-ANNA"){
              nn.sites<-try(site.match(Test=env.data()[which(env.data()[,"Sites"]%in%i),-c(1)],
                                   Reference=env.data()[which(env.data()[,"Sites"]%in%sel.ref()),-c(1)],
                                   k= if (ab.k!=0 & ab.k<nrow(env.data()[which(env.data()[,"Sites"]%in%sel.ref()),-c(1)])) ab.k else NULL,
                                   adaptive=ab.adaptive,
                                   RDA.reference=isolate(bio.data.t$Summary.Metrics[which(rownames(bio.data.t$Summary.Metrics)%in%sel.ref()),ab.sel.mets])),
                            silent=T)
              if(is(nn.sites,"try-error")){
                results[i,1]<-nn.sites[1]
                next
              } else {
                nn.sites$method1<-NULL
                nn.sites$method1<-"RDA-ANNA"
              }
            }
            if (nnmethod=="ANNA"){
              
              nn.sites<-try(site.match(Test=env.data()[which(env.data()[,"Sites"]%in%i),-c(1)],
                                   Reference=env.data()[which(env.data()[,"Sites"]%in%sel.ref()),-c(1)],
                                   k= if (ab.k!=0 & ab.k<nrow(env.data()[which(env.data()[,"Sites"]%in%sel.ref()),-c(1)])) ab.k else NULL,
                                   adaptive=ab.adaptive),
                            silent=T)
              if(is(nn.sites,"try-error")){
                results[i,1]<-nn.sites[1]
                next
              } else {
                nn.sites$method1<-NULL
                nn.sites$method1<-"ANNA"
                
              }
            }
            if (nnmethod=="User Selected"){
              if (!is.null(env.data())) {
                nn.sites<-try(site.match(Test=env.data()[which(env.data()[,"Sites"]%in%i),-c(1)],
                                     Reference=env.data()[which(env.data()[,"Sites"]%in%sel.ref()),-c(1)],
                                     k= if (ab.k!=0 & ab.k<nrow(env.data()[which(env.data()[,"Sites"]%in%sel.ref()),-c(1)])) ab.k else NULL,
                                     adaptive=ab.adaptive,
                                     RDA.reference = if (isolate(input$nnmethod.user=="RDA-ANNA")) {bio.data.t$Summary.Metrics[which(rownames(bio.data.t$Summary.Metrics)%in%sel.ref()),ab.sel.mets]} else {NULL}),
                              silent=T)
                if(is(nn.sites,"try-error")){
                  results[i,1]<-nn.sites[1]
                  next
                } else {
                  nn.sites$method1<-NULL
                  nn.sites$method1<-"User Selected with distance"
                }
              } else {
                nn.sites<-NULL
                nn.sites$final.dist<-paste(unlist(user.site.match()[which(user.site.match()[,1]%in%i),2:length(which(user.site.match()[which(user.site.match()[,1]%in%i),]!=""))]))
                names(nn.sites$final.dist)<-nn.sites$final.dist
                nn.sites$method<-NULL
                nn.sites$method<-"User Selected"
                nn.sites$method1<-NULL
                nn.sites$method1<-"User Selected"
                nn.sites
                
              }
            }
            
            bio.data.t1<-NULL
            bio.data.t1$Summary.Metrics<-bio.data.t$Summary.Metrics[c(names(nn.sites$final.dist),i),]
            if(nnmethod!="RDA-ANNA"){
              a.met<-add.met(Test=bio.data.t$Raw.Data[i,],Reference=bio.data.t$Raw.Data[names(nn.sites$final.dist),],original=F)
              bio.data.t1$Summary.Metrics<-cbind(bio.data.t1$Summary.Metrics, a.met)
            }

            if (input$metdata==F){
              tsa.results<-try(tsa.test(bio.data.t1$Summary.Metrics[i,ab.sel.mets],
                                    Reference=bio.data.t1$Summary.Metrics[names(nn.sites$final.dist),ab.sel.mets],
                                    distance= if (ab.distance) nn.sites$final.dist else NULL,
                                    outlier.rem= ab.outlier.rem,
                                    m.select= ab.m.select,
                                    na.cutoff=0.7),
                               silent=T)
              
              #tsa.results<-try(tsa.test(Test=add.met(Test=bio.data.t$Raw.Data[i,],Reference=bio.data.t$Raw.Data[names(nn.sites$final.dist),])[(1+length(nn.sites$final.dist)),ab.sel.mets],
              #                      Reference=add.met(Test=bio.data.t$Raw.Data[i,],Reference=bio.data.t$Raw.Data[names(nn.sites$final.dist),])[1:length(nn.sites$final.dist),ab.sel.mets],
              #                      distance= if (ab.distance) nn.sites$final.dist else NULL,
              #                      outlier.rem= ab.outlier.rem,
              #                      m.select= ab.m.select,
              #                      na.cutoff=0.7),
              #                 silent=T)
              if(is(tsa.results,"try-error")){
                results[i,1]<-tsa.results[1]
                next
              } else {
              }
            } else {
              tsa.results<-try(tsa.test(bio.data.t1$Summary.Metrics[i,ab.sel.mets],
                                        Reference=bio.data.t1$Summary.Metrics[names(nn.sites$final.dist),ab.sel.mets],
                                        distance= if (ab.distance) nn.sites$final.dist else NULL,
                                        outlier.rem= ab.outlier.rem,
                                        m.select= ab.m.select,
                                        na.cutoff=0.7),
                               silent=T)
              
              #taxa.data<-rbind(bio.data.t1$Summary.Metrics.Data[names(nn.sites$final.dist),],bio.data.t1$Summary.Metrics[i,])
              #tsa.results<-try(tsa.test(Test=taxa.data[(1+length(nn.sites$final.dist)),ab.sel.mets],
              #                      Reference=taxa.data[1:length(nn.sites$final.dist),ab.sel.mets],
              #                      distance= if (ab.distance) nn.sites$final.dist else NULL,
              #                      outlier.rem= ab.outlier.rem,
              #                      m.select= ab.m.select,
              #                      na.cutoff=0.7),
              #                 silent=T)
              if(is(tsa.results,"try-error")){
                results[i,1]<-tsa.results[1]
                next
              } else{
              }
            }
            results[i,]<-c(tsa.results$tsa.results[1,],tsa.results$tsa.results[2,],tsa.results$tsa.results[3,],tsa.results$tsa.results[4,],
                           tsa.results$tsa.results[5,],tsa.results$tsa.results[6,],tsa.results$tsa.results[7,],tsa.results$tsa.results[8,],
                           tsa.results$tsa.results[9,],tsa.results$general.results[5,],tsa.results$general.results[6,],nn.sites$method1,
                           tsa.results$jacknife[1,],tsa.results$general.results[3,],tsa.results$general.results[4,],tsa.results$general.results[2,])
            
            if (any(new.dirs==T)){
              if (isolate(input$ab.nnscatter.plot)){
                if (!is.null(env.data())){
                  jpeg(filename=paste0(sel.dir(),"/NN Ordination/",i,"-nnord.jpeg"),width=640,height=480)
                  sitematch.plot(nn.sites)
                  dev.off()
                }
              }
              if (isolate(input$ab.nndist.plot)){
                if (nnmethod!="User Selected") {
                  jpeg(filename=paste0(sel.dir(),"/NN Distance/",i,"-nndist.jpeg"),width=640,height=480)
                  plot(nn.sites)
                  dev.off()
                }
              }
              if (isolate(input$ab.tsadist.plot)){
                jpeg(filename=paste0(sel.dir(),"/TSA Distance/",i,"-tsadist.jpeg"),width=640,height=480)
                plot(tsa.results)
                dev.off()
              }
              if (isolate(input$ab.tsabox.plot)){
                if (input$metdata==T){
                  jpeg(filename=paste0(sel.dir(),"/TSA Boxplot/",i,"-tsabox.jpeg"),width=640,height=480)
                  boxplot(tsa.results)
                  dev.off()
                } 
                if (input$metdata==F){
                  tsa.stand<-tsa.zscore(Test=bio.data.t1$Summary.Metrics[i,],
                                        Reference=bio.data.t1$Summary.Metrics[rownames(tsa.results$raw.data[-c(nrow(tsa.results$raw.data)),]),])
                  nInd<-ncol(tsa.stand)
                  nRef<-nrow(tsa.stand)-1
                  part.tsa<-if (!is.null(tsa.results$partial.tsa)) {tsa.results$partial.tsa} else {NULL}
                  all.met<-colnames(tsa.stand)
                  sel.met<-unlist(strsplit(substr(tsa.results$general.results["Selected Indicator Metrics",],1,(nchar(tsa.results$general.results["Selected Indicator Metrics",])-2)),split=", "))
                  
                  cols<-colorRampPalette(brewer.pal(12, "Paired"))(nInd)
                  text<-paste(seq(1:ncol(tsa.stand)),colnames(tsa.stand),sep=".")
                  b1<-ceiling(length(text)/3)
                  b2<-ceiling(length(text)*2/3)

                  l<-rbind(c(1,1,1),c(1,1,1),c(2,3,4))
                  jpeg(filename=paste0(sel.dir(),"/TSA Boxplot/",i,"-tsabox.jpeg"),width=640,height=480)
                  
                  layout(l)
                  
                  par(mar = c(1.9,0.8,1.2,0.8))
                  boxplot(tsa.stand[1:nRef,],col=cols,outline=F,yaxt="n",ylim=c(min(tsa.stand)*1.3,max(tsa.stand)*1.1),names=seq(1:nInd),cex.axis=1.2,main="")
                  title(main=paste0(rownames(tsa.stand)[(nRef+1)]," Boxplot"),cex=1.5)
                  points(seq(1:nInd),tsa.stand[(nRef+1),],col="red",pch=19,cex=1)
                  
                  points(which(colnames(tsa.stand)%in%sel.met),tsa.stand[nrow(tsa.stand),sel.met],col="black",pch="O",cex=1.75)#This line circles points that are used in analysis
                  if (any(part.tsa$p<0.05)) {
                    points(which(colnames(tsa.stand)%in%rownames(part.tsa)[part.tsa$p<0.05]),rep((min(tsa.stand)*1.2),length(rownames(part.tsa)[part.tsa$p<0.05])),col="red",pch="*",cex=2)
                  }
                  par(mar = c(0,0,0,0))
                  plot(1, type="n", axes=F, xlab="", ylab="")
                  legend("center",text[1:b1],cex=1.25,fill=cols[1:b1],bty="n",x.intersp=1,y.intersp=1)
                  par(mar = c(0,0,0,0))
                  plot(1, type="n", axes=F, xlab="", ylab="")
                  legend("center",text[(b1+1):b2],cex=1.25,fill=cols[(b1+1):b2],bty="n",x.intersp=1,y.intersp=1)
                  par(mar = c(0,0,0,0))
                  plot(1, type="n", axes=F, xlab="", ylab="")
                  legend("center",text[(b2+1):length(text)],cex=1.25,fill=cols[(b2+1):length(text)],bty="n",x.intersp=1,y.intersp=1)
                  dev.off()
                }
              }
              if (isolate(input$ab.tsascatter.plot)){
                jpeg(filename=paste0(sel.dir(),"/TSA Ordination/",i,"-tsaord.jpeg"),width=640,height=480)
                pcoa.tsa(tsa.results)
                dev.off()
              }
              if (isolate(input$ab.cascatter.plot)){
                jpeg(filename=paste0(sel.dir(),"/CA Ordination/",i,"-caord.jpeg"),width=640,height=480)
                Reference<-bio.data.t$Raw.Data[names(tsa.results$mahalanobis.distance)[1:(length(tsa.results$mahalanobis.distance)-1)],]
                nRef<-nrow(Reference)
                Test<-bio.data.t$Raw.Data[names(tsa.results$mahalanobis.distance)[(length(tsa.results$mahalanobis.distance))],]
                raw.data<-rbind(Reference,Test)
                pRef<-colSums(decostand(Reference,"pa"))/nrow(Reference)
                
                ca.ord<-cca(log(raw.data[,names(which(pRef>=0.1))]+1))
                ca1<-ca.ord$CA$u[,1]
                ca2<-ca.ord$CA$u[,2]
                
                plot(ca.ord,type="n",main=paste(rownames(ca1)[(nRef+1)]," CA Plot",sep=""),
                     xlab=paste("CA1 ",substr((eigenvals(ca.ord)[1]/sum(eigenvals(ca.ord)))*100,1,4),"%"),
                     ylab=paste("CA2 ",substr((eigenvals(ca.ord)[2]/sum(eigenvals(ca.ord)))*100,1,4),"%"))
                text(x=ca.ord$CA$v[,1],y=ca.ord$CA$v[,2],labels=rownames(ca.ord$CA$v),col="grey50",cex=0.7)
                points(x=ca1[1:nRef],y=ca2[1:nRef],cex=0.8,col="black",pch=19)
                points(x=ca1[(nRef+1)],y=ca2[(nRef+1)],cex=0.8,col="red",pch=19)
                text(x=ca1[1:nRef],y=ca2[1:nRef],labels=names(ca1)[1:nRef],col="black",cex=0.8,pos=3)
                text(x=ca1[(nRef+1)],y=ca2[(nRef+1)],labels=names(ca1)[(nRef+1)],col="red",cex=0.9,pos=3)
                suppressWarnings(ordiellipse(ca.ord,c(rep(1,nrow(Reference)),0),kind="sd",conf=0.95,draw="line",col="grey20",lty=5,show.groups=1))
                dev.off()
              }
              if (isolate(input$ab.multi.plot)){
                if (input$metdata==T){
                  l2<-rbind(c(1,2,2,2),c(1,2,2,2),c(3,2,2,2),c(4,4,5,5),c(4,4,5,5),c(6,6,7,7),c(6,6,7,7))
                } else {
                  l2<-rbind(c(1,2,2,2),c(1,2,2,2),c(1,3,4,5),c(6,6,7,7),c(6,6,7,7),c(8,8,9,9),c(8,8,9,9))
                }
                pdf(file=paste0(sel.dir(),"/Multiplot/",i,"-multi.pdf"),width=8,height=11)
                
                layout(l2)
                txt<-c("Test Sample:",i,paste0("Status - ", tsa.results$tsa.results[1,]),"","Reference samples:",rownames(tsa.results$raw.data[-c(nrow(tsa.results$raw.data)),]))
                par(mar = c(0,0,0,0))
                textplot(txt,halign="left", valign="top")

                if (input$metdata==T){
                  boxplot(tsa.results)
                } else {
                  tsa.stand<-tsa.zscore(Test=bio.data.t1$Summary.Metrics[i,],
                                        Reference=bio.data.t1$Summary.Metrics[rownames(tsa.results$raw.data[-c(nrow(tsa.results$raw.data)),]),])
                  nInd<-ncol(tsa.stand)
                  nRef<-nrow(tsa.stand)-1
                  part.tsa<-if (!is.null(tsa.results$partial.tsa)) {tsa.results$partial.tsa} else {NULL}
                  all.met<-colnames(tsa.stand)
                  sel.met<-unlist(strsplit(substr(tsa.results$general.results["Selected Indicator Metrics",],1,(nchar(tsa.results$general.results["Selected Indicator Metrics",])-2)),split=", "))
                  
                  cols<-colorRampPalette(brewer.pal(12, "Paired"))(nInd)
                  text<-paste(seq(1:ncol(tsa.stand)),colnames(tsa.stand),sep=".")
                  b1<-ceiling(length(text)/3)
                  b2<-ceiling(length(text)*2/3)
                  
                  par(mar = c(1.9,0.8,1.2,0.8))
                  boxplot(tsa.stand[1:nRef,],col=cols,outline=F,yaxt="n",ylim=c(min(tsa.stand)*1.3,max(tsa.stand)*1.1),names=seq(1:nInd),cex.axis=1.2,main="")
                  title(main=paste0(rownames(tsa.stand)[(nRef+1)]," Boxplot"),cex=1.5)
                  points(seq(1:nInd),tsa.stand[(nRef+1),],col="red",pch=19,cex=1)
                  
                  points(which(colnames(tsa.stand)%in%sel.met),tsa.stand[nrow(tsa.stand),sel.met],col="black",pch="O",cex=1.75)#This line circles points that are used in analysis
                  if (any(part.tsa$p<0.05)) {
                    points(which(colnames(tsa.stand)%in%rownames(part.tsa)[part.tsa$p<0.05]),rep((min(tsa.stand)*1.2),length(rownames(part.tsa)[part.tsa$p<0.05])),col="red",pch="*",cex=2)
                  }
                  par(mar = c(0,0,0,0))
                  plot(1, type="n", axes=F, xlab="", ylab="")
                  legend("center",text[1:b1],cex=1,fill=cols[1:b1],bty="n",x.intersp=0.85,y.intersp=0.85)
                  par(mar = c(0,0,0,0))
                  plot(1, type="n", axes=F, xlab="", ylab="")
                  legend("center",text[(b1+1):b2],cex=1,fill=cols[(b1+1):b2],bty="n",x.intersp=0.85,y.intersp=0.85)
                  par(mar = c(0,0,0,0))
                  plot(1, type="n", axes=F, xlab="", ylab="")
                  legend("center",text[(b2+1):length(text)],cex=1,fill=cols[(b2+1):length(text)],bty="n",x.intersp=0.85,y.intersp=0.85)
                }
                
                if (isolate(input$multiplot1.sel=="None")){
                  plot(1, type="n", axes=F, xlab="", ylab="")
                }
                if (isolate(input$multiplot1.sel=="Nearest-Neighbour Ordination")){
                  sitematch.plot(nn.sites)
                }
                if (isolate(input$multiplot1.sel=="Nearest Neighbour Distance")){
                  plot(nn.sites)
                }
                if (isolate(input$multiplot1.sel=="TSA Distance")){
                  plot(tsa.results)
                }
                if (isolate(input$multiplot1.sel=="TSA Ordination")){
                  pcoa.tsa(tsa.results)
                }
                if (isolate(input$multiplot1.sel=="CA Ordination")){
                  Reference<-bio.data.t$Raw.Data[names(tsa.results$mahalanobis.distance)[1:(length(tsa.results$mahalanobis.distance)-1)],]
                  nRef<-nrow(Reference)
                  Test<-bio.data.t$Raw.Data[names(tsa.results$mahalanobis.distance)[(length(tsa.results$mahalanobis.distance))],]
                  raw.data<-rbind(Reference,Test)
                  pRef<-colSums(decostand(Reference,"pa"))/nrow(Reference)
                  
                  ca.ord<-cca(log(raw.data[,names(which(pRef>=0.1))]+1))
                  ca1<-ca.ord$CA$u[,1]
                  ca2<-ca.ord$CA$u[,2]
                  
                  plot(ca.ord,type="n",main=paste(rownames(ca1)[(nRef+1)]," CA Plot",sep=""),
                       xlab=paste("CA1 ",substr((eigenvals(ca.ord)[1]/sum(eigenvals(ca.ord)))*100,1,4),"%"),
                       ylab=paste("CA2 ",substr((eigenvals(ca.ord)[2]/sum(eigenvals(ca.ord)))*100,1,4),"%"))
                  text(x=ca.ord$CA$v[,1],y=ca.ord$CA$v[,2],labels=rownames(ca.ord$CA$v),col="grey50",cex=0.7)
                  points(x=ca1[1:nRef],y=ca2[1:nRef],cex=0.8,col="black",pch=19)
                  points(x=ca1[(nRef+1)],y=ca2[(nRef+1)],cex=0.8,col="red",pch=19)
                  text(x=ca1[1:nRef],y=ca2[1:nRef],labels=names(ca1)[1:nRef],col="black",cex=0.8,pos=3)
                  text(x=ca1[(nRef+1)],y=ca2[(nRef+1)],labels=names(ca1)[(nRef+1)],col="red",cex=0.9,pos=3)
                  suppressWarnings(ordiellipse(ca.ord,c(rep(1,nrow(Reference)),0),kind="sd",conf=0.95,draw="line",col="grey20",lty=5,show.groups=1))
                }
                if (isolate(input$multiplot2.sel=="None")){
                  plot(1, type="n", axes=F, xlab="", ylab="")
                }
                if (isolate(input$multiplot2.sel=="Nearest-Neighbour Ordination")){
                  sitematch.plot(nn.sites)
                }
                if (isolate(input$multiplot2.sel=="Nearest Neighbour Distance")){
                  plot(nn.sites)
                }
                if (isolate(input$multiplot2.sel=="TSA Distance")){
                  plot(tsa.results)
                }
                if (isolate(input$multiplot2.sel=="TSA Ordination")){
                  pcoa.tsa(tsa.results)
                }
                if (isolate(input$multiplot2.sel=="CA Ordination")){
                  Reference<-bio.data.t$Raw.Data[names(tsa.results$mahalanobis.distance)[1:(length(tsa.results$mahalanobis.distance)-1)],]
                  nRef<-nrow(Reference)
                  Test<-bio.data.t$Raw.Data[names(tsa.results$mahalanobis.distance)[(length(tsa.results$mahalanobis.distance))],]
                  raw.data<-rbind(Reference,Test)
                  pRef<-colSums(decostand(Reference,"pa"))/nrow(Reference)
                  
                  ca.ord<-cca(log(raw.data[,names(which(pRef>=0.1))]+1))
                  ca1<-ca.ord$CA$u[,1]
                  ca2<-ca.ord$CA$u[,2]
                  
                  plot(ca.ord,type="n",main=paste(rownames(ca1)[(nRef+1)]," CA Plot",sep=""),
                       xlab=paste("CA1 ",substr((eigenvals(ca.ord)[1]/sum(eigenvals(ca.ord)))*100,1,4),"%"),
                       ylab=paste("CA2 ",substr((eigenvals(ca.ord)[2]/sum(eigenvals(ca.ord)))*100,1,4),"%"))
                  text(x=ca.ord$CA$v[,1],y=ca.ord$CA$v[,2],labels=rownames(ca.ord$CA$v),col="grey50",cex=0.7)
                  points(x=ca1[1:nRef],y=ca2[1:nRef],cex=0.8,col="black",pch=19)
                  points(x=ca1[(nRef+1)],y=ca2[(nRef+1)],cex=0.8,col="red",pch=19)
                  text(x=ca1[1:nRef],y=ca2[1:nRef],labels=names(ca1)[1:nRef],col="black",cex=0.8,pos=3)
                  text(x=ca1[(nRef+1)],y=ca2[(nRef+1)],labels=names(ca1)[(nRef+1)],col="red",cex=0.9,pos=3)
                  suppressWarnings(ordiellipse(ca.ord,c(rep(1,nrow(Reference)),0),kind="sd",conf=0.95,draw="line",col="grey20",lty=5,show.groups=1))
                }
                if (isolate(input$multiplot3.sel=="None")){
                  plot(1, type="n", axes=F, xlab="", ylab="")
                }
                if (isolate(input$multiplot3.sel=="Nearest-Neighbour Ordination")){
                  sitematch.plot(nn.sites)
                }
                if (isolate(input$multiplot3.sel=="Nearest Neighbour Distance")){
                  plot(nn.sites)
                }
                if (isolate(input$multiplot3.sel=="TSA Distance")){
                  plot(tsa.results)
                }
                if (isolate(input$multiplot3.sel=="TSA Ordination")){
                  pcoa.tsa(tsa.results)
                }
                if (isolate(input$multiplot3.sel=="CA Ordination")){
                  Reference<-bio.data.t$Raw.Data[names(tsa.results$mahalanobis.distance)[1:(length(tsa.results$mahalanobis.distance)-1)],]
                  nRef<-nrow(Reference)
                  Test<-bio.data.t$Raw.Data[names(tsa.results$mahalanobis.distance)[(length(tsa.results$mahalanobis.distance))],]
                  raw.data<-rbind(Reference,Test)
                  pRef<-colSums(decostand(Reference,"pa"))/nrow(Reference)
                  
                  ca.ord<-cca(log(raw.data[,names(which(pRef>=0.1))]+1))
                  ca1<-ca.ord$CA$u[,1]
                  ca2<-ca.ord$CA$u[,2]
                  
                  plot(ca.ord,type="n",main=paste(rownames(ca1)[(nRef+1)]," CA Plot",sep=""),
                       xlab=paste("CA1 ",substr((eigenvals(ca.ord)[1]/sum(eigenvals(ca.ord)))*100,1,4),"%"),
                       ylab=paste("CA2 ",substr((eigenvals(ca.ord)[2]/sum(eigenvals(ca.ord)))*100,1,4),"%"))
                  text(x=ca.ord$CA$v[,1],y=ca.ord$CA$v[,2],labels=rownames(ca.ord$CA$v),col="grey50",cex=0.7)
                  points(x=ca1[1:nRef],y=ca2[1:nRef],cex=0.8,col="black",pch=19)
                  points(x=ca1[(nRef+1)],y=ca2[(nRef+1)],cex=0.8,col="red",pch=19)
                  text(x=ca1[1:nRef],y=ca2[1:nRef],labels=names(ca1)[1:nRef],col="black",cex=0.8,pos=3)
                  text(x=ca1[(nRef+1)],y=ca2[(nRef+1)],labels=names(ca1)[(nRef+1)],col="red",cex=0.9,pos=3)
                  suppressWarnings(ordiellipse(ca.ord,c(rep(1,nrow(Reference)),0),kind="sd",conf=0.95,draw="line",col="grey20",lty=5,show.groups=1))
                }
                if (isolate(input$multiplot4.sel=="None")){
                  plot(1, type="n", axes=F, xlab="", ylab="")
                }
                if (isolate(input$multiplot4.sel=="Nearest-Neighbour Ordination")){
                  sitematch.plot(nn.sites)
                }
                if (isolate(input$multiplot4.sel=="Nearest Neighbour Distance")){
                  plot(nn.sites)
                }
                if (isolate(input$multiplot4.sel=="TSA Distance")){
                  plot(tsa.results)
                }
                if (isolate(input$multiplot4.sel=="TSA Ordination")){
                  pcoa.tsa(tsa.results)
                }
                if (isolate(input$multiplot4.sel=="CA Ordination")){
                  Reference<-bio.data.t$Raw.Data[names(tsa.results$mahalanobis.distance)[1:(length(tsa.results$mahalanobis.distance)-1)],]
                  nRef<-nrow(Reference)
                  Test<-bio.data.t$Raw.Data[names(tsa.results$mahalanobis.distance)[(length(tsa.results$mahalanobis.distance))],]
                  raw.data<-rbind(Reference,Test)
                  pRef<-colSums(decostand(Reference,"pa"))/nrow(Reference)
                  
                  ca.ord<-cca(log(raw.data[,names(which(pRef>=0.1))]+1))
                  ca1<-ca.ord$CA$u[,1]
                  ca2<-ca.ord$CA$u[,2]
                  
                  plot(ca.ord,type="n",main=paste(rownames(ca1)[(nRef+1)]," CA Plot",sep=""),
                       xlab=paste("CA1 ",substr((eigenvals(ca.ord)[1]/sum(eigenvals(ca.ord)))*100,1,4),"%"),
                       ylab=paste("CA2 ",substr((eigenvals(ca.ord)[2]/sum(eigenvals(ca.ord)))*100,1,4),"%"))
                  text(x=ca.ord$CA$v[,1],y=ca.ord$CA$v[,2],labels=rownames(ca.ord$CA$v),col="grey50",cex=0.7)
                  points(x=ca1[1:nRef],y=ca2[1:nRef],cex=0.8,col="black",pch=19)
                  points(x=ca1[(nRef+1)],y=ca2[(nRef+1)],cex=0.8,col="red",pch=19)
                  text(x=ca1[1:nRef],y=ca2[1:nRef],labels=names(ca1)[1:nRef],col="black",cex=0.8,pos=3)
                  text(x=ca1[(nRef+1)],y=ca2[(nRef+1)],labels=names(ca1)[(nRef+1)],col="red",cex=0.9,pos=3)
                  suppressWarnings(ordiellipse(ca.ord,c(rep(1,nrow(Reference)),0),kind="sd",conf=0.95,draw="line",col="grey20",lty=5,show.groups=1))
                }
                
                dev.off()

              }
            }
            }
        })
      }
      write.csv(results,file=paste0(sel.dir(),"/results.csv"))
      results
    } else {
      return(NULL)
    }
  })
  
  observeEvent(input$ab.go,{
    results()
  })
  
  
  output$ab.results<-renderDataTable({
    temp<-cbind(rownames(results()),results()[,1:15])
    colnames(temp)[1]<-"Site"
    temp
  })
  
  output$abdone<-reactive({
    if (!is.null(results())){
      return(1)
    }
  })
  
  output$seldir<-reactive({
    abselmet<-ab.sel.mets()
    if (sel.dir()!="No Directory Selected" & !is.null(abselmet) & length(abselmet)>2){
      return(1)
    } else {
      return(0)
    }
  })
  
  output$inprogress<-reactive({
    input$ab.go
    if (!is.null(results())){
      return(1)
    } 
  })
  
  
  outputOptions(output, 'abdone', suspendWhenHidden=FALSE)
  outputOptions(output, 'seldir', suspendWhenHidden=FALSE)
  outputOptions(output, 'inprogress', suspendWhenHidden=FALSE)
  
  output$multiplot1<-renderUI({
    opts<-c("None",
            if (!is.null(env.data())) {"Nearest-Neighbour Ordination"},
            if (input$nnmethod!="User Selected"){"Nearest Neighbour Distance"},
            if (input$metdata==F) {"CA Ordination"},
            "TSA Distance", "TSA Ordination")
    
    wellPanel(selectInput("multiplot1.sel", label = "Middle Left", choices = opts, selected = "NONE"),
              selectInput("multiplot2.sel", label = "Middle Right", choices = opts, selected = "NONE"),
              selectInput("multiplot3.sel", label = "Lower Left", choices = opts, selected = "NONE"),
              selectInput("multiplot4.sel", label = "Lower Right", choices = opts, selected = "NONE")
    )
  })
  
#NEED TO INCLUDE THIS LAST LINE TO MAKE SURE IT CLOSES IN STANDALONE VERSION
#session$onSessionEnded(function() { 
#     stopApp()
#     q("no") 
#    })

})


