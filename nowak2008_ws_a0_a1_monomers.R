
library(igraph)
library(Matrix)
library(deSolve)
library(pheatmap)
library(foreach)
library(doParallel)
library(RColorBrewer)
llen <- 7

totV <- 0
for(i in 1:llen){
  totV <- totV + 2^i
}

totV <- totV/2
treex <- igraph::make_tree(totV,children = 2,mode = "out")

start<-"0"
nodes<-NULL
nodes[1]<-start
c<-unlist(neighborhood(treex,nodes="1",mode="out"))[-1]
for(i in 2:length(degree(treex))){
  if(length(c)>0){
    p<-unlist(neighborhood(treex,nodes=c[1],mode = "in"))[-1]
    if(length(p)>0){
      start<-nodes[p[1]]
    }
    nodes[c[1]]<-paste(start,"0",sep = "")
    nodes[c[2]]<-paste(start,"1",sep = "")
    c<-unlist(neighborhood(treex,nodes=as.character(i),mode="out"))[-1]
  }
}

treex0<-treex
nodes0<-nodes




treex<-igraph::make_tree(totV,children = 2,mode = "out")
flag<-0

start<-"1"
nodes<-NULL
nodes[1]<-start
c<-unlist(neighborhood(treex,nodes="1",mode="out"))[-1]
for(i in 2:length(degree(treex))){
  if(length(c)>0){
    p<-unlist(neighborhood(treex,nodes=c[1],mode = "in"))[-1]
    if(length(p)>0){
      start<-nodes[p[1]]
    }
    nodes[c[1]]<-paste(start,"0",sep = "")
    nodes[c[2]]<-paste(start,"1",sep = "")
    c<-unlist(neighborhood(treex,nodes=as.character(i),mode="out"))[-1]
  }
}


treex1<-treex
nodes1<-nodes
set_vertex_attr(graph=treex1,name = "label",value = nodes1)
set_vertex_attr(graph=treex0,name = "label",value = nodes0)

L1<-cbind(nodes1[as_edgelist(treex1)[,1]],nodes1[as_edgelist(treex1)[,2]])
L0<-cbind(nodes0[as_edgelist(treex0)[,1]],nodes0[as_edgelist(treex0)[,2]])


#stoichmatrix +2 in the columns for ctivating monomers
N<-matrix(0,nrow=totV*2,ncol=2*nrow(L1)+2)
rownames(N)<-c(nodes1,nodes0)
for(i in 1:nrow(L1)){
  w<-which(rownames(N)==L1[i,1])
  N[w,i]<- -1
  
  w<-which(rownames(N)==L1[i,2])
  N[w,i]<- 1
  
}






cnt<-nrow(L1)

for(i in 1:nrow(L0)){
  w<-which(rownames(N)==L0[i,1])
  N[w,i+cnt]<- -1
  
  w<-which(rownames(N)==L0[i,2])
  N[w,i+cnt]<- 1
  
}

posUno<-which(rownames(N)=="1")
posZero<-which(rownames(N)=="0")


N[posUno,nrow(L0)+nrow(L1)+1]<-1
N[posZero,nrow(L0)+nrow(L1)+2]<-1

#last rate degradations (same for all)
gamma<-1
a<-1
#
ai1<-.5
ai0<-.5
a1<-1/2
a0<-1/2
A<-rep(ai1,times=ncol(N))
A[seq(1,length(A),by=2)]<-ai0
A[nrow(L0)+nrow(L1)+1]<-a1
A[nrow(L0)+nrow(L1)+2]<-a0

uno<-1
zero<-1
degRates<-rep(x=gamma,times=nrow(N))

rates <- c(A,degRates)
 idM<-diag(rep(x=-1,times=nrow(N)))
NnoDeg<-N
N<-cbind(N,idM)
#should be full rank
#i.e. rk = min(dim(N))
rk<-rankMatrix(N)
#N[1,]<-0
#N[totV+1,]<-0
#check dimensions
Ax<- N %*% rates
#should give n strings = totV*2
dim(Ax)

#rates<-10*matrix(1/nchar(c(L1[,1],L0[,1],0.0)),nrow=ncol(N),ncol=1)
Pars<-list(stoich=N,rates=rates,total=totV,monomers=c(uno,zero))

sist<-function(Time,State,Pars)
{
  #State[c(1,64)]<-1
  with(as.list(c(State,Pars)),{
    #State[c(1,total+1)]<-1
    w<-which(stoich<0,arr.ind=TRUE)    
    v <- State[w[,1]]*rates[w[,2]]
    posUno<-which(rownames(stoich)=="1")
    posZero<-which(rownames(stoich)=="0")
    oioi<-c(monomers[1] * rates[which(stoich[posUno,]>0)],monomers[2] * rates[which(stoich[posZero,]>0)])
    v<-c(v[1:(which(stoich[posUno,]>0)-1)],oioi,v[(which(stoich[posUno,]>0)):length(v)])
    #v<-matrix(c(v,deg),nrow=ncol(stoich),ncol=1)
    #stoich[1,]<-0
    #stoich[total+1,]<-0
    x <- stoich %*% v
    x[x<0]<-0
    #x[c(1,total+1)]<-1
    list(x)
  })
}

yini  <- matrix(0,nrow=totV*2,ncol=1)
yini[1,1]<-0
yini[totV+1,1]<-0
names(yini)<-rownames(N)

times <- seq(0, 30,.25)


outSuperSymm   <- ode(yini, times, sist, Pars,method = "bdf")
#####

ParsAsymm<-Pars
p<-sample(1:ncol(NnoDeg),0.5*ncol(NnoDeg))

#####
multi<-2
#####

ParsAsymm$rates[p]<-ParsAsymm$rates[p]*multi
outAsymm1   <- ode(yini, times, sist, ParsAsymm,method = "bdf")
outmin<-outSuperSymm[,-1]
outminAs1<-outAsymm1[,-1]
cx<-rownames(N)
#cx<-cx[-c(1,64)]
png("SuperSymm_vs_Asymm_boxplot.png", width = 5, height = 5, units = 'in', res = 300)

par(mfrow=c(1,2))
boxplot(log10(outmin[nrow(outmin),])~nchar(cx),xlab = "String length",ylab="Abundance")
boxplot(log10(outminAs1[nrow(outminAs1),])~nchar(cx),xlab = "String length",ylab="Abundance")
dev.off()
CM<-colorRampPalette(RColorBrewer::brewer.pal(8,"YlGnBu"))

pdf("heatmap_SuperSymm.pdf")
DF<-data.frame(L=factor(nchar(cx)))
rownames(DF)<-rownames(N)
C<-brewer.pal(length(unique(nchar(cx))),"Blues")
names(C)<-1:max(nchar(cx))
colore<-list(L=C)

pheatmap(log10(min(outmin[outmin>0])+t(outmin[1:20,])),
         cluster_cols = FALSE,
         scale = "none",
         fontsize=8,
         color = CM(50),#cluster_rows = FALSE,
         annotation_row = DF,
         annotation_colors = colore)

dev.off()
pdf("heatmap_Asymm.pdf")


pheatmap(log10(min(outminAs1[outmin>0])+t(outminAs1)),
         cluster_cols = FALSE,
         scale = "none",
         fontsize=8,
         fontsize_row = 8,
         annotation_row = DF,
         annotation_colors = colore,
         color = CM(50))#,cluster_rows = FALSE)
dev.off()
#summary(outSuperSymm)

C<-brewer.pal(max(nchar(cx))+1,"Greens")
names(C)<-1:7
png("SuperSymm_vs_Asymm.png", width = 5, height = 5, units = 'in', res = 300)
#png(file="SuperSymm_vs_Asymm.png",res = 100)
par(mfrow=c(1,2))

plot(outSuperSymm[,1],
     min(outmin[outmin>0])+outmin[,1],
     ylim=c(min(outmin[outmin>0]),
            max(outmin)),
     log="y",type="l",
     col=C[nchar(cx[1])],
     main="Super-Symmetric",ylab="Abundance of Strings",xlab="Time")

for(i in 1:ncol(outmin)){
  lines(outSuperSymm[,1],outmin[,i],col=C[nchar(cx[i])])
}






plot(outAsymm1[,1],
     min(outminAs1[outminAs1>0])+outminAs1[,1],
     ylim=c(min(outminAs1[outminAs1>0]),
            max(outminAs1)),
     log="y",type="l",
     col=C[nchar(cx[1])],
     main="Asymmetric",ylab="Abundance of Strings",xlab="Time")

for(i in 1:ncol(outminAs1)){
  lines(outAsymm1[,1],outminAs1[,i],col=C[nchar(cx[i])])
}
dev.off()



#create the cluster
my.cluster <- parallel::makeCluster(
  8, 
  type = "FORK"
)

#check cluster definition (optional)
print(my.cluster)
#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)














startDeg<-0.1
#####
times <- seq(0, 20,2)

nrep<-100
increm<- 0.05
degra2test<-startDeg+(1:nrep) * increm
Pars_degscan<-Pars
LL_degscan <- foreach(i = 1:nrep,.combine="c") %dopar% {
  
  
  
  
  #for(i in 1:nrep){
  Pars_degscan$rates[(ncol(NnoDeg)+1):ncol(N)]<-degra2test[i]
  tmp   <- ode(yini, times, sist, Pars_degscan,method = "bdf")
  print(i)
  cbind(tmp[nrow(tmp),-1])
  
  
}


ss_degscan<-matrix(LL_degscan,nrow=nrep,ncol=2*totV,byrow = TRUE)
s_degscan<-degra2test
par(mfrow=c(1,1))
Cx<-brewer.pal(7,"Greens")
names(Cx)<-1:max(nchar(cx))
i<-1
plot(s_degscan,(ss_degscan[,i]/max(ss_degscan[,i])),
     col=Cx[which(names(Cx)==factor(nchar(cx[i])))],
     ylab="Steady State level",
     xlab="degradation rate for all strings",
     main="Steady State level",
     ylim=c(min((ss_degscan[ss_degscan>0]))/max(ss_degscan[,i]),1),
     xlim=c(0.1,max(s_degscan)),log="xy",type="b",pch=18)

for(i in 1:ncol(ss_degscan)){
  lines(s_degscan,jitter(ss_degscan[,i]/max(ss_degscan[,i]),factor = 0.01),
        col = Cx[which(names(Cx)==factor(nchar(cx[i])))],
        type = "b",pch=18,
        cex = 0.6)
}

par(mfrow=c(2,2))

target=4
plot(nchar(cx),ss_degscan[target,],
     col=Cx[which(names(Cx)==factor(nchar(cx[i])))],
     ylab="Steady State level",
     xlab="Length of strings",
     main=paste("d=",s_degscan[target],sep = ""),
     ylim=c(min((ss_degscan[ss_degscan>0]))/max(ss_degscan[,i]),1),
     xlim=c(2,6),log="xy",type="p",pch=18,cex=2)
corrLow<-cor.test(log10(nchar(cx)),log10(ss_degscan[target,]))

target=18
plot(nchar(cx),ss_degscan[target,],
     col=Cx[which(names(Cx)==factor(nchar(cx[i])))],
     ylab="Steady State level",
     xlab="Length of strings",
     main=paste("d=",s_degscan[target],sep = ""),
     ylim=c(min((ss_degscan[ss_degscan>0]))/max(ss_degscan[,i]),1),
     xlim=c(2,6),log="xy",type="p",pch=18,cex=2)
corrMedLow<-cor.test(log10(nchar(cx)),log10(ss_degscan[target,]))


target=28
plot(nchar(cx),ss_degscan[target,],
     col=Cx[which(names(Cx)==factor(nchar(cx[i])))],
     ylab="Steady State level",
     xlab="Length of strings",
     main=paste("d=",s_degscan[target],sep = ""),
     ylim=c(min((ss_degscan[ss_degscan>0]))/max(ss_degscan[,i]),1),
     xlim=c(2,6),log="xy",type="p",pch=18,cex=2)
corrMedium<-cor.test(log10(nchar(cx)),log10(ss_degscan[target,]))

target=95
plot(nchar(cx),ss_degscan[target,],
     col=Cx[which(names(Cx)==factor(nchar(cx[i])))],
     ylab="Steady State level",
     xlab="Length of strings",
     main=paste("d=",s_degscan[target],sep = ""),
     ylim=c(min((ss_degscan[ss_degscan>0]))/max(ss_degscan[,i]),1),
     xlim=c(2,6),log="xy",type="p",pch=18,cex=2)
corrStrong<-cor.test(log10(nchar(cx)),log10(ss_degscan[target,]))


allcorr<-cbind(c(corrLow$estimate,corrMedLow$estimate,corrMedium$estimate,corrStrong$estimate),
               c(corrLow$p.value,corrMedLow$p.value,corrMedium$p.value,corrStrong$p.value),
               s_degscan[c(4,18,28,95)])
rownames(allcorr)<-c("Low","Med.Low","Medium","Strong")
colnames(allcorr)<-c("R","pvalue","Deg. rate")
####varying power of s



pctg<-c(0.01,0.1,0.4)


for(m in 1:length(pctg)){
  
ParsAsymm<-Pars

intt<-intersect(which(N>0,arr.ind=TRUE)[,2],which(nchar(rownames(N))<7))
p<-sample(intt,pctg[m]*length(intt))

#####

#####
times <- seq(0, 20,.5)


mmm<-seq(-4,2,.1)

multi<-10^mmm
#multi[multi!=0]<-0
nrep<-length(multi)
ss<-matrix(0,nrow=nrep,ncol=length(cx))
s<-NULL

LL <- foreach(i = 1:nrep,.combine="c") %dopar% {
 

ParsAsymm<-Pars

#for(i in 1:nrep){
  ParsAsymm$rates[p]<-Pars$rates[p]+multi[i]
  tmp   <- ode(yini, times, sist, ParsAsymm,method = "ode45")
  print(i)
  cbind(tmp[nrow(tmp),-1])
  
  
}


ss<-matrix(LL,nrow=nrep,ncol=2*totV,byrow = TRUE)
colnames(ss)<-rownames(N)

ss<-ss[,-which(nchar(colnames(ss))==7)]
s<-multi

C<-brewer.pal(7,"Accent")
par(mfrow=c(1,2))
fname<-paste("selection_no_rep_",pctg[m],".pdf",sep="",collapse = "")
pdf(fname)
plot(s,1E-06+(ss[,1])/rowSums(ss),
     col=C[nchar(colnames(ss)[1])],
     ylab="Steady State level",
     xlab="s",

          ylim=c(min(ss[ss>0]),max(ss)),xlim=c(min(multi),max(multi)),type="l",log="xy")

for(i in 2:ncol(ss)){
  lines(s,1E-06+ss[,i]/rowSums(ss),
       col=C[nchar(colnames(ss)[i])],
       type="l",pch=18)
  
  
  
}



#####
#####

ParsAsymm<-Pars
  
  ratesx<-rates
  ratesx[1:ncol(NnoDeg)]<-sample(1:1000,ncol(NnoDeg))/1000
  changeRates<-ratesx

LLx <- foreach(i = 1:nrep,.combine="c") %dopar% {
  changeRates<-ratesx
  
  changeRates[p]<-changeRates[p]+multi[i]
  
  #for(i in 1:nrep){
  ParsAsymm$rates<-changeRates
  tmp   <- ode(yini, times, sist, ParsAsymm,method = "ode45")
  print(i)
  cbind(tmp[nrow(tmp),-1])
  
  
}


ss_diffRates<-matrix(LLx,nrow=nrep,ncol=2*totV,byrow = TRUE)

colnames(ss_diffRates)<-rownames(N)
ss_diffRates<-ss_diffRates[,-which(nchar(colnames(ss_diffRates))==7)]

plot(s,1E-06+(ss_diffRates[,1])/rowSums(ss_diffRates),
     xlim=c(min(s),max(s)),
     ylab="Steady State level",
     xlab="s",

          ylim=c(min(1E-06+(ss_diffRates[ss_diffRates>0])/rowSums(ss_diffRates)),
            max(1E-06+(ss_diffRates)/rowSums(ss_diffRates))),
     type="l",
     col=C[nchar(colnames(ss_diffRates)[1])],
     log="xy")

for(i in 2:ncol(ss_diffRates)){
  lines(s,1E-06+(ss_diffRates[,i])/rowSums(ss_diffRates),
        col=C[nchar(colnames(ss_diffRates)[i])],
        type="l")



}

dev.off()
}





par(mfrow=c(1,2))
w<-which(ss_diffRates[1,]>ss_diffRates[nrow(ss_diffRates),])

plot(s,1E-06+(ss_diffRates[,w[1]])/rowSums(ss_diffRates),
     xlim=c(min(s),max(s)),
     ylab="Steady State level",
     xlab="Rate multiplier for 0.1 of reactions",
     main="Scan of s Vs Steady State levels",
     ylim=c(min(1E-06+(ss_diffRates[ss_diffRates>0])/rowSums(ss_diffRates)),
            max(1E-06+(ss_diffRates)/rowSums(ss_diffRates))),
     type="l",
     col=C[nchar(colnames(ss_diffRates)[w[1]])],
     log="xy")

for(i in 2:length(w)){
  lines(s,1E-06+(ss_diffRates[,w[i]])/rowSums(ss_diffRates),
        col=C[nchar(colnames(ss_diffRates)[w[i]])],
        type="l")
  
  
  
}






















w<-which(ss_diffRates[1,]<ss_diffRates[nrow(ss_diffRates),])

plot(s,1E-06+(ss_diffRates[,w[1]])/rowSums(ss_diffRates),
     xlim=c(min(s),max(s)),
     ylab="Steady State level",
     xlab="s",
    # main=,
     ylim=c(min(1E-06+(ss_diffRates[ss_diffRates>0])/rowSums(ss_diffRates)),
            max(1E-06+(ss_diffRates)/rowSums(ss_diffRates))),
     type="l",
     col=C[nchar(colnames(ss_diffRates)[w[1]])],
     log="xy")

for(i in 2:length(w)){
  lines(s,1E-06+(ss_diffRates[,w[i]])/rowSums(ss_diffRates),
        col=C[nchar(colnames(ss_diffRates)[w[i]])],
        type="l")
  
  
  
}

























