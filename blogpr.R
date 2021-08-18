#
#
#
blogpr<-function(p, sh=FALSE){
  require(parallel)
  require(stringi)
  require(limma)
  require(reshape2)
  cpl=colorRampPalette(c("#001030","#0050AA","#10AA10","#FFFF00","#FFA000","#B50000"))
  lf=list.files(path = p)
  lf=lf[lf!="key.csv"]
  key=read.csv(paste(p,"key.csv", sep=""))
  N=nrow(key)
  fls=lapply(lf,function(f){
    finp=paste(p,f,sep = "")
    fn=stri_extract(f,regex="\\S+(?=\\.)")
    fcon=file(finp)
    f2l=readLines(con=fcon, n=2)
    nlns=as.double(stri_extract(f2l[2], regex="[1-9]+"))
    fhead=readLines(con = fcon, n=nlns+2)
    x=read.delim(finp, skip=nlns+2, header=T, check.names = F, stringsAsFactors = FALSE)
    list(fn,x)
  })
  closeAllConnections()
  flns=lapply(fls,function(f){f[[1]]})
  fls=lapply(fls,function(f){f[[2]]})
  names(fls)=flns
  bl=lapply(seq_along(key$file), function(i){
    f=fls[[key$file[i]]] 
    if (key$ch[i]=="G") f=f[f$Block==key$block[i],c("Column","Row","Name","ID","F532 Median","Flags")] else f=f[f$Block==key$block[i],c(c("Column","Row","Name","ID","F635 Median","Flags"))]
    colnames(f)=c("C", "R","Name","ID","V","Fl")
    key$file=stri_extract_first(key$file, regex=".+(?<=_)\\d+")                                    
    return(list(key$pat[i],key$diag[i],key$ch[i],f,key$file[i],key$zone[i])) #1-patient/2-diag/3-channel/4-[1-Col/2-Row/3-Name/4-ID/5-F635 Med/6-Flags]/5-file/6-zone
    })
  
  coor=lapply(bl, function(x){list("Chip"=unlist(x[5]),"Zone"=unlist(x[6]), "Diagnosis"=unlist(x[2]),"Channel"=unlist(x[3]),"Row"=x[4][[1]][,2], "Column"=x[4][[1]][,1])})  #1-chip/2-zone/3-diag/4-channel/5-row/6-column
  if (all(unlist(lapply(coor, function(x){unlist(coor[[1]][5:6])==unlist(x[5:6])})))) print("Ordered, homogeneous Sets") else return("Error: Nonhomogeneous Sets!")
  fl=lapply(bl,function(i){i[4][[1]][,6]})
  ID=lapply(bl,function(i){i[4][[1]][,4]})
  Nm=lapply(bl,function(i){i[4][[1]][,3]})
  if (all(unlist(lapply(ID, function(x){ID[[1]]==x})))) ID=ID[[1]] else return("Error: Different Maps!")
  if (all(!is.na(unlist(Nm)))) {
    if (all(unlist(lapply(Nm, function(x){Nm[[1]]==x})))) Nm=Nm[[1]] else return("Error: Different Names!")}
  pt=lapply(bl,function(i){as.character(i[[1]])})
  dg=lapply(bl,function(i){as.character(i[[2]])})
  cn=paste(pt,dg, sep="_")
  
  
  dgch=t(as.data.frame(lapply(coor, function(x){c(as.character(x$Diagnosis),as.character(x$Channel))})))
  colnames(dgch)=c("D","C")
  Res0=as.data.frame(lapply(bl, function(x){x[[4]]$V}), col.names=paste(pt,dg, sep="_"))
  meanB=list()
  meanB[[1]]=rowMeans(Res0[,dgch[,1]=="E"&dgch[,2]=="R"])
  meanB[[2]]=rowMeans(Res0[,dgch[,1]=="E"&dgch[,2]=="G"])
  names(meanB)=c("R","G")
  chips=unique(unlist(lapply(coor,function(x){x$Chip})))
  R0=Res0[,dgch[,1]!="E"]
  #R0[]=0
  inf0 =list()
  co=1
  for (c in chips) {
    cx=sapply(coor, function(ci){ci$Chip==c&ci$Diagnosis!="E"})  # which are the reads of the stained zones of this chip 
    cnx=sapply(coor, function(ci){ci$Chip==c&ci$Diagnosis=="E"}) # which are the reads of the same chip before staining
    zone=unique(sapply(coor[cx], function(x){x[[2]]}))           # which zones exist for it
    chn=sapply(coor,function(x){x[x$Chip==c]$Channel})           # which channel was this chip read in
    chn=unique(unlist(chn[!sapply(chn,is.null)]))
    
    # channel of stained and background images should be the same, if rescanning in a different channel the alternative channel files have to be removed assuming each chip has been read in only one channel
    if (length(chn)>1) return("Two color chip processing has not been implemented yet!") 
    
    for (z in zone){                                            # For each zone determine empty chip background read 
      zi=sapply(coor, function(i){i$Zone==z})                   # (either a single read or the mean of several reads or, when missing, the mean of other similar)
      i=cnx&zi
      if (any(i)) {
          if (length(i[i])>1) {
            bkg=rowMeans(Res0[,i])
            NegCol=which(i)
            }
        else {
            bkg=Res0[,i]
            NegCol=which(i)
            }
          }
      else {
        bkg=meanB[as.character(chn)]
        NegCol=0
        }
      
      bkg=unlist(bkg)
      #print(list(c,z,chn,head(bkg)))
      #readline()
      ip=cx&zi
      #if (max(which(ip))>ncol(R0)) return("R0 dimension error!")
      j=which(ip)
      if (length(j)>1){
        pt=proc.time()
      }
      #j=j[!is.na(j[1:ncol(R0)])]
      for (t in j) {                     # Subtract the empty chip background
        R0[,co]=Res0[,t]-bkg 
        inf0[[co]]=list(Chip=c,Zone=z,Channel=as.character(chn),Patient_Diag=colnames(Res0)[t],PosCol=t, NegCol=NegCol) # 1.Chip;2.Zone;3.Channel;4.Diagnosis;5.Positive Column; 6.Negative Column
        co=co+1
      }
    }
  }
  colnames(R0)=unlist(lapply(inf0,function(x){x[[4]]}))
  # R0=apply(R0,2,function(x){x=x-min(x)+1})          
  # rownames(R0)=ID
  cl <- makeCluster(7)
  ex <- Filter(function(x) is.function(get(x, .GlobalEnv)), ls(.GlobalEnv))
  clusterExport(cl, ex)
  clusterExport(cl, list("R0", "coor"), envir=environment())
  B0=parApply(cl,R0, 2, function(clm){
      #gprl1(clm, coor[[1]]$Row, coor[[1]]$Column, ID, shw = TRUE)
      IM=data.frame(R=coor[[1]]$Row, C=coor[[1]]$Column, V=clm)
      ID=data.frame(R=coor[[1]]$Row, C=coor[[1]]$Column, V=ID, stringsAsFactors = FALSE)
      ID[,1]=as.double(ID[,1])
      ID[,2]=as.double(ID[,2])
      bgmy(IM,ID)[,3]
  })                                                           
  stopCluster(cl)
  if (sh==TRUE) {
  for (i in 1:ncol(R0)) {
      xS=data.frame(R=coor[[1]]$Row, C=coor[[1]]$Column, V=R0[,i])
      xB=data.frame(R=coor[[1]]$Row, C=coor[[1]]$Column, V=B0[,i])
      xF=xS
      xF[,3]=xF[,3]-xB[,3]
      xF[,3]=xF[,3]-min(xF[,3])+1
      zr=range(R0)
      img=acast(xS, C~R, value.var = "V")
      image(img, main=c("Signal", colnames(R0)[i]), zlim=zr, col=cpl(256))
      img=acast(xB, C~R, value.var = "V")
      image(img, main=c("Background", colnames(R0)[i]), zlim=zr, col=cpl(256))
      img=acast(xF, C~R, value.var = "V")
      image(img, main=c("Filtered", colnames(R0)[i]), zlim=zr, col=cpl(256))
      hist(xF[,3], breaks = 100)
      rd=readline()
    }
  }
    res=sapply(seq_along(R0[1,]),function(i){
    backgroundCorrect.matrix(R0[,i], B0[,i], method="normexp", normexp.method="mle")
  })
  colnames(res)=colnames(R0)
  return(list(coor,ID,res,fl, Nm, inf0))
}