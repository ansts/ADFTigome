require(limma)
require(parallel)
require(reshape2)
require(d3heatmap)
require(umap)
require(plot3D)
require(gplots)require(rgl)
require(plot3D)
require(corrplot)
require(factoextra)
require(FactoMineR)
require(stringdist)
require(qualV)
require(igraph)
require(pbapply)
require(stringi)
require(matrixStats)


require(multimode)
require(multcomp)
require(Biostrings)
library(lsmeans)
require(boot)
require(Rfast)
require(sva)
require(mixtools)
require(clusterCrit)
require(e1071)
require(abind)






cpl=colorRampPalette(c("#0010303F","#0050AA9F","#10AA109F","#FFFF009F","#FFA0009F","#B500009F"), alpha=T)

#
# Collect the data from the .gpr files in a folder /gpr containing also a key.csv file describing the slides
#

alldata=chpr("gpr//", sh=T)
Coor=alldata[[1]] # 1-chip/2-zone/3-diag/4-channel/5-row/6-column
iD=alldata[[2]]   # ID of spots
Res=alldata[[3]]  # Results - locally normalized data (the background is reconstituted using the duplicate diff. and subtracted using normexp)
rownames(Res)=iD  # 
Fl=alldata[[4]]   # Flags data
Nam=alldata[[5]]  # Names data (not used)

coix=as.data.frame(t(sapply(alldata[[6]], function(x){as.character(x)})), stringsAsFactors=FALSE) # 1.Chip;2.Zone;3.Channel;4.Patient ID _ Diagnosis;5.Positive Column; 6.Negative Column
FLdf=as.data.frame(Fl[1:ncol(Res)])
FL=rowSums(FLdf)==0
Cop=Coor[as.double(coix[,5])]

WD=Res[FL,]
WDa=aggregate.data.frame(WD, by=list(rownames(WD)), FUN=mean)
pp=rownames(Dn)
l=length(pp)

WDam=data.frame(log10(WDa[,2:ncol(WDa)]), stringsAsFactors = FALSE, row.names = WDa[,1])
# normalize for amino acid composition dependent binding (e.g. non-specific stickiness of charged amino acids)
D=pepnorm(WDam)
# normalize arrays 
Dn=normalizeCyclicLoess(D, method="affy", span=.15, iterations=3) 

fdn=rep(1,5)
fdn[1:2]=2
fdn[3:4]=3
fdn=factor(fdn)
calls=limfit_(fdn,Dn, p=0.05)

# single sample tests

smp=c("AD1","AD2","FTD1","FTD2","IVIgM")
calls=sapply(1:4, function(i){
  fdn=rep(1,5)
  fdn[i]=2
  if (i %in% c(1,3)) fdn[i+1]=3 else fdn[i-1]=3
  fdn=as.factor(c(fdn))
  print(fdn)
  x=limfit_(fdn, Dn)
  print(x)
  png(file=paste("vulcanostdres",smp[i],".png", collapse = "", sep="_"), width = 15, height=15, units="cm", res=300)
  plot(x$logFC, -log(x$P.Value), col=(x$adj.P.Val<0.1)+1, pch=16, cex=(x$adj.P.Val<0.1)*0.5+0.25, ylab="-log(p)", xlab="Log ratio", ylim=c(0,9))
  dev.off()
  cls=rownames(x)[x$adj.P.Val<0.1]
  if (length(cls)>0) {
    png(file=paste("traces_best", smp[i], ".png", collapse = "", sep="_"), width = 15, height=15, units="cm", res=300)
    for (n in cls){
      col=which(cls==n)  
      plot(v[n,], ty="l", ylim=range(v), xaxt="n", xlab="", col=col, ylab="Log Intensity")
      par(new=T)
    }
    par(new=F)
    xtick<-seq(1, 5, by=1)
    axis(side=1, at=xtick, labels = FALSE)
    text(x=xtick,  par("usr")[3], 
         labels =smp, pos = 1, xpd = TRUE)
    legend("top",cls, text.col=seq_along(cls), cex=.5)
    dev.off()
  }
  return(cls)
})



# Selection based on graph of reactivity profiles
v=scale(Dn)
colnames(v)=c("AD1","AD2","FTD1","FTD2","IVIgM")
vm=melt(v)
colnames(vm)=c("S","P","V")

# Perform ANOVA on the values of the different patient pools for each pair of peptides  
cl=makeCluster(4)
clusterExport(cl, list("pp","vm","l"), envir = environment())

Fs=pbsapply(pp, function(p1){
      sapply(pp, function(p2){
      y=vm[vm$S %in% c(p1,p2),]
      ly=lm(V~P, data=y)
      Fv=anova(ly)$`F value`[1]
    })
}, cl=cl)
stopCluster(cl)
diag(Fs)=0

cl=makeCluster(4)
clusterExport(cl, list("pp","vm","l"), envir = environment())

Ps=pbsapply(pp, function(p1){
  sapply(pp, function(p2){
    y=vm[vm$S %in% c(p1,p2),]
    ly=lm(V~P, data=y)
    anova(ly)$`Pr(>F)`[1]
  })
}, cl=cl)
stopCluster(cl)
pPs=p.adjust(unlist(Ps[lower.tri(Ps)]))
table(pPs)

Osa=sapply(pp, function(p){
  stringdist(p,pp, method = "osa")
})


png(file="density_by_OSA.png", width = 15, height=15, units="cm", res=300)
for (i in 2:7) {
  plot(density(log10(Fs[lower.tri(Fs)&Osa==i])), xlab="log10(F)", ylim=c(0,1),xlim=range(log10(Fs[lower.tri(Fs)])), col=cpl(7)[i], main="")
  par(new=T)
}
legend("topleft", paste("OSA distance = ",2:7,  sep=""), text.col=cpl(7)[2:7], cex=.5)
lines(c(0.5,0.5), c(0,1))
par(new=F)
dev.off()

png(file=paste("traces", ".png", collapse = "", sep="_"), width = 15, height=15, units="cm", res=300)
cq=c("TFQARSM","QIAWIRQ")
for (p in cq){
  col=which(cq==p)  
  plot(v[p,], ty="l", ylim=range(v), xaxt="n", xlab="", col=col, ylab="Log Intensity")
  par(new=T)
}
par(new=F)
xtick<-seq(1, 5, by=1)
axis(side=1, at=xtick, labels = FALSE)
text(x=xtick,  par("usr")[3], 
     labels =colnames(v), pos = 1, xpd = TRUE)
legend("top",legend=cq, text.col=seq_along(cq), cex=.5)
dev.off()

DOsa=lapply(2:7, function(i)log10(Fs[lower.tri(Fs)&Osa==i]))
plot(t(sapply(DOsa, quantile, c(0.5, 0.9))), cex=log10(lengths(DOsa)))

Lpsel=cbind(pp[which(Fs>3.16, arr.ind = T)[,1]],pp[which(Fs>3.16, arr.ind = T)[,2]])

ps=unique(c(Lpsel))
  
FsAdjM=Fs
FsAdjM[FsAdjM<3.16]=0
diag(FsAdjM)=0

GFs=graph_from_adjacency_matrix(FsAdjM,mode="undirected", weighted = T)
Gfs=simplify(GFs)
cGFs=components(GFs)
GFs=induced.subgraph(GFs, vids=V(GFs)[cGFs$membership==1])

write.graph(GFs, format = "graphml", file="GFs.graphml")

mxcq3=max_cliques(GFs, min=3)

png(file="cliqDist.png", width = 15, height=15, units="cm", res=300)
tmxq3=table(lengths(mxcq3))
nms=as.numeric(tmxq3)
barplot(tmxq3, xlab="Clique Сize", ylab = "Counts")
dev.off()

Dnmxcq=t(sapply(mxcq3, function(m){
  colMedians(Dn[m,])
}))

limfit_(fdn, Dnmxcq, p=0.05)

Dndmxcq=t(sapply(mxcq3, function(m){
  colMedians(Dn[m,1:4]-Dn[m,5])
}))

plot(rowMeans(Dndmxcq[,1:2]),rowMeans(Dndmxcq[,3:4]), cex=0.75, pch=16, col=rgb(0,0,0,0.5))
lines(c(-0.6,0.4), c(0,0))
lines(c(0,0), c(-0.6,0.4))

Dndmxcqp0=sapply(mxcq3, function(m){
  x=Dn[m,1:4]-Dn[m,5]
  x=melt(x)
  x$Var2=unlist(stri_extract_all(x$Var2, regex="(?<=_)\\w+")) # Label by diagnosis
  x=lm(value~Var2, data=x)
  x=anova(x)
  return(x$`Pr(>F)`[1])
})

Dndmxcqp=p.adjust(Dndmxcqp0)
sumAD=rowSums(Dndmxcq[,1:2])
sumFTD=rowSums(Dndmxcq[,3:4])
cqsig=which(Dndmxcqp<0.1)
cqsigAD=which(Dndmxcqp<0.1 & (sumAD>sumFTD))
cqsigFTD=which(Dndmxcqp<0.1 & (sumAD<sumFTD))
colcqsig=sapply((Dndmxcqp<0.1), function(y) {
  rgb(y*1,0,0,0.5)
})
plot(rowMeans(Dndmxcq[,1:2]),rowMeans(Dndmxcq[,3:4]), xlim=c(-0.8,0.8), ylim=c(-0.8,0.8), col=colcqsig, pch=16, cex=abs(xy$residuals)+0.3)
col=(Dndmxcqp<0.05)+(Dndmxcqp<0.1)+1
png(file="vulcanomxcq.png", width = 15, height=15, units="cm", res=300)
plot(log((sumAD-min(sumAD)+1)/(sumFTD-min(sumFTD)+1)),-log(Dndmxcqp0), pch=16, cex=(col>1)*0.2+0.3, col=col, ylab="-Log(p)", xlab="Log(AD/FTD)")
dev.off()



callsAD=mxcq3[cqsigAD]
callsFTD=mxcq3[cqsigFTD]
ppcqsigAD=unique(unlist(callsAD))
ppcqsigADhi=ppcqsigAD[(rowMeans(Dn[ppcqsigAD,1:2])-Dn[ppcqsigAD,5])>0]
ppcqsigADlo=ppcqsigAD[(rowMeans(Dn[ppcqsigAD,1:2])-Dn[ppcqsigAD,5])<=0]
ppcqsigFTD=unique(unlist(callsFTD))
ppcqsigFTDhi=ppcqsigFTD[(rowMeans(Dn[ppcqsigFTD,3:4])-Dn[ppcqsigFTD,5])>0]
ppcqsigFTDlo=ppcqsigFTD[(rowMeans(Dn[ppcqsigFTD,3:4])-Dn[ppcqsigFTD,5])<=0]
vcqsigADhi=(names(V(GFs)) %in% ppcqsigADhi)*4
vcqsigADlo=(names(V(GFs)) %in% ppcqsigADlo)*8
vcqsigFTDhi=(names(V(GFs)) %in% ppcqsigFTDhi)*2
vcqsigFTDlo=(names(V(GFs)) %in% ppcqsigFTDlo)*1
vcqsig=vcqsigADhi+vcqsigADlo+vcqsigFTDhi+vcqsigFTDlo
GFs=set.vertex.attribute(GFs,"cqsig",value=vcqsig)
write.graph(GFs, format = "graphml", file="GFs.graphml")
ppofint=aggregate(names(V(GFs)), by=list(vcqsig), FUN=list)
names(ppofint[[2]])=ppofint[[1]]
ppofint=ppofint[[2]]
mxcqofint=lapply(ppofint[2:7],function(l){
  mxcq3[unique(unlist(lapply(l,grep,mxcq3)))]
})

ppofintt=sapply(ppofint[2:5],function(l){
  x=60-length(l)
  c(l,rep(0,x))
})
write.csv(ppofintt, "ppOI.csv", col.names = F, row.names = F)

boxplot(v[unlist(callsAD),], notch=T, xlab='Sample', varwidth=T, ylab="Scaled Log Intensities", main="Calls AD > FTD")
boxplot(v[unlist(callsFTD),], notch=T, xlab='Sample', varwidth=T,ylab="Scaled Log Intensities", main="Calls FTD > AD")
boxplot(v[unlist(ppcqsigADhi),], notch=T, xlab='Sample', varwidth=T,ylab="Scaled Log Intensities",main="AD High")
boxplot(v[unlist(ppcqsigADlo),], notch=T, xlab='Sample', varwidth=T,ylab="Scaled Log Intensities",main="AD Low")
boxplot(v[unlist(ppcqsigFTDlo),], notch=T, xlab='Sample', varwidth=T,ylab="Scaled Log Intensities",main="FTD Low")
boxplot(v[unlist(ppcqsigFTDhi),], notch=T, xlab='Sample', varwidth=T,ylab="Scaled Log Intensities",main="FTD High")

ppil=melt(ppofint[2:5])
GoFI=induced.subgraph(GFs, ppil[,1])
GoFI=simplify(GoFI)
write.graph(GoFI, format = "graphml", file="GoFI.graphml")

# clique graph
cqOI=mxcq3[cqsig]
cqOIAdjMx=sapply(cqOI, function(c1){
          sapply(cqOI, function(c2){
            length(intersect(c1,c2))/length(union(c1,c2))
  })
})
diag(cqOIAdjMx)=0
cqOIF=graph_from_adjacency_matrix(cqOIAdjMx, mode="undirected", weighted = T)
cqOIF=simplify(cqOIF)
cqlab=rep(0,length(cqsig))
cqlab[cqsig %in% cqsigAD]=1
cqlab[cqsig %in% cqsigFTD]=2
cqOIF=set_vertex_attr(cqOIF,"D", value=cqlab)
write.graph(cqOIF, format = "graphml", file="cqOIF.graphml")
cqOIFnodes=read.csv(file="cqOIF.nodesxx.csv")
tcqOIFnmc=table(cqOIFnodes$modularity_class)
cqi=as.numeric(names(tcqOIFnmc[tcqOIFnmc>1]))
cqmod=sapply(cqi, function(i){
  j=cqOIFnodes$modularity_class==i
  unique(unlist(cqOI[j]))
})

names(cqmod)=cqi
sapply(seq_along(cqmod), function(i){
  l=cqmod[[i]]
  l=t(sapply(l,function(p){
    c(paste(">",p, collapse="", sep=""),p)
  }))
  write.csv(l, file=paste("cqmod_",cqi[i],".csv", collapse="", sep=""), quote = F, row.names = FALSE)
} )

