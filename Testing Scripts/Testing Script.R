#TEST

#load datasets
data(YKEnvData,envir = environment()) #Biological dataset
data(YKBioData,envir = environment()) #Environmental dataset

#Calculate indicator metrics from raw biological data
bio.data.test<-benth.met(YKBioData,2,2)

#Extract just the summary metrics
bio.data<-bio.data.test$Summary.Metrics

#standardize row names between datasets
rownames(YKEnvData)<-bio.data.test$Site.List

#Match a test site (#201) to the nearest neighbour reference set
nn.sites<-site.match(YKEnvData[201,-c(1)],YKEnvData[1:118,-c(1)],k=NULL,adaptive=T)

#Calculate additional metrics based on selected Reference sites
taxa.data<-add.met(Test=bio.data.test$Raw.Data[201,],Reference=bio.data.test$Raw.Data[names(nn.sites$final.dist),])


#TSA test of indicator metrics at test site and reference sites selected used site.match()
tsa.results<-tsa.test(Test=taxa.data[nrow(taxa.data),],Reference=taxa.data[names(nn.sites$final.dist),],distance=nn.sites$final.dist, outlier.rem=T, m.select=T)
tsa.results

met.sel<-metric.select(Test=tsa.results$raw.data[nrow(tsa.results$raw.data),],Reference=tsa.results$raw.data[names(nn.sites$final.dist),],outlier.rem = F, rank = F, outbound = 0.1)


#Evaluate Results
boxplot(tsa.results)
plot(tsa.results)
pcoa.tsa(tsa.results)




