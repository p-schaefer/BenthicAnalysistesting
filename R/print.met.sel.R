print.met.sel<-function(met.sel){
  cat("Selected Metrics:\n")
  print(met.sel$Best.Metrics)
  #cat("\n")
  #cat("All indicative metrics:\n")
  #print(met.sel$Indicative.Metrics)
}