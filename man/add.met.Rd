% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add.met.R
\name{add.met}
\alias{add.met}
\title{Reference specific Benthic Metric Calculation}
\usage{
add.met(Test, Reference, original = F, tax.fields = 2)
}
\arguments{
\item{Test.taxa}{Output of \code{\link[BenthicAnalysistesting]{benth.met}}$Raw.Data for the test site}

\item{Reference.taxa}{Output of \code{\link[BenthicAnalysistesting]{benth.met}}$Raw.Data for the Reference site}
}
\value{
Indicator metrics caluclated from \code{\link[BenthicAnalysistesting]{benth.met}} in addition to O:E ratios and Bray-Curtis Distance.
}
\description{
Calculation of a variety of metrics for determing impairment of Benthic Macroinvertebrate Communities, based on available Reference Condition information.
}
\examples{
#load datasets
data(YKEnvData,envir = environment()) #Biological dataset
data(YKBioData,envir = environment()) #Environmental dataset

#Calculate indicator metrics from raw biological data
bio.data.test<-benth.met(YKBioData,2,2)

#standardize row names between datasets
rownames(YKEnvData)<-rownames(bio.data.test$Site.List)

#Match a test site (#201) to the nearest neighbour reference set
nn.sites<-site.match(YKEnvData[201,],YKEnvData[1:118,],k=NULL,adaptive=T)

#Extract the raw taxa data for additional metric calculation
taxa.data<-rbind(bio.data.test$Raw.Data[names(nn.sites$final.dist),],bio.data.test$Raw.Data[rownames(bio.data.test$Raw.Data[201,]),])

#Calculate additional metrics based on nearest neighbour reference sites
additional.metrics<-add.met(Test=bio.data.test$Raw.Data[rownames(bio.data.test$Raw.Data[201,]),],Reference=bio.data.test$Raw.Data[names(nn.sites$final.dist),])
additional.metrics
}
\keyword{Benthic}
\keyword{Metrics}
