# A wrapper for the limma linear model based extraction of calls

limfit_<-function(fdsn,vals, n=1000, p=0.05){
  require(limma)
  require(sva)
  require(parallel)
  m=length(unique(fdsn))
  ncl=ncol(vals)
  if (ncl!=length(fdsn)) {return("Error: ncol of design and matrix differ!")}
  clin=character(0)
  design = model.matrix(~fdsn)
  cn=paste("f",0:(m-1),sep ="")
  colnames(design)=cn
  fit=lmFit(vals,design)
  fit1=eBayes(fit)
  tptb=topTable(fit1, coef=2, number=n) #,adjust.method="holm"
  nr=nrow(tptb[tptb$adj.P.Val<p,])
  clin=rownames(tptb[1:nr,])
  #res=rbind(tptb1[1,])
  return(tptb) 
}  

