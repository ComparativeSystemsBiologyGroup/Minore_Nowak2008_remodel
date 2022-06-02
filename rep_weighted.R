

sistRepl6_weighted<-function(Time,State,Pars)
{
  #State[c(1,64)]<-1
  State[which(State < 0)] <- 0
  with(as.list(c(State,Pars)),{
    #State[c(1,total+1)]<-1
    w<-which(stoich<0,arr.ind=TRUE)    
    v <- State[w[,1]]*rates[w[,2]]
    posUno<-which(rownames(stoich)=="1")
    posZero<-which(rownames(stoich)=="0")
    oioi<-c(monomers[1] * rates[which(stoich[posUno,]>0)],monomers[2] * rates[which(stoich[posZero,]>0)])
    names(v)<-NULL
    v<-c(v[1:(which(stoich[posUno,]>0)-1)],oioi,v[(which(stoich[posUno,]>0)):length(v)])
    nc<-nchar(rownames(stoich))
    n_ones <- nc - nchar(gsub(rownames(stoich),pattern = "1",replacement = ""))
    n_zeros <- nc-nchar(gsub(rownames(stoich),pattern = "0",replacement = ""))
                                          
                                          #replication rates
    sumfix<-sum(fitness*State*(n_ones*(State[posUno]/sumx)+
                                 n_zeros*(State[posZero]/sumx)))
    sumx<-sum(State)
    w<-which(stoich>0,arr.ind=TRUE)    
    if(sumx==0){
      vr<-matrix(0,nrow=ncol(stoich),ncol=1) 
    }else{
      replica <- State[w[,1]] * replicator[w[,1]] * (fitness[w[,1]]*(
        n_ones*(State[posUno]/sumx)+
        n_zeros*(State[posZero]/sumx)) 
        - sumfix/sumx)
      
      vr<-matrix(0,nrow=ncol(stoich),ncol=1)      
      
      vr[w[,2]]<-replica
    }
    x <- stoich %*% (v+vr)
    x[which(x<0)]<-0
    #x[c(1,total+1)]<-1
    list(x)
  })
}



whorep<-matrix(1,nrow=nrow(N),ncol=1)
cx<-rownames(N)
startRep<-0
#####
times <- seq(0, 50,1)

nrep<-50
increm<- .1
replic<-startRep+(1:nrep) * increm
ss_replication<-matrix(0,nrow=nrep,ncol=length(cx))
s_replication<-NULL




my.cluster <- parallel::makeCluster(
  4,
  type = "FORK"
)

#check cluster definition (optional)
print(my.cluster)
#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)
ratesx<-rates
ratesx[1:ncol(NnoDeg)]<-runif(ncol(NnoDeg))
w6<-which(nchar(rownames(N))==6)



fi<-matrix(0,nrow=nrow(N),ncol=1)

fi[w6,1]<-runif(length(w6))
#fi[w6[10],1]<-0.999

LL_replication_W <- foreach(i = 1:nrep,.combine="c") %dopar% {
  
  
  Pars<-list(stoich=N,
             rates=rates,
             total=totV,
             replicator=whorep*replic[i], 
             L=nchar(rownames(N)),
             fitness=fi,
             monomers=c(uno,zero))
  
  tmp   <- ode(yini, times, sistRepl6_weighted, Pars,method = "ode45")
  
  print(i)
  cbind(tmp[nrow(tmp),-1])
  
  
}

ss_replication_W<-matrix(LL_replication_W,nrow=nrep,ncol=2*totV,byrow = TRUE)
s_replication_W<-replic
cx<-rownames(N)
library(RColorBrewer)
Cx<- brewer.pal(7,"Dark2")
names(Cx)<-1:max(nchar(cx))
i<-1


ss_replication_W<-ss_replication_W[,-which(nchar(cx)==7)]

cx<-cx[-which(nchar(cx)==7)]
#ss_replication[ss_replication<1E-6]<-0
i=1
plot(s_replication_W,ss_replication_W[,i]/rowSums(ss_replication_W),
     col=Cx[which(names(Cx)==factor(nchar(cx[i])))],
     ylab="Steady State level",
     xlab="Replication selection strength",
     main="Prevolution to evolution",
     #ylim=c(1E-05,1),
     ylim=c(min(ss_replication_W/rowSums(ss_replication_W)),
            max(ss_replication_W/rowSums(ss_replication_W))),
     xlim=c(min(s_replication_W),max(s_replication_W)),log="xy",type="l",pch=18)

for(i in 2:ncol(ss_replication_W)){
  lines(s_replication_W,1E-07+(ss_replication_W[,i]/rowSums(ss_replication_W)),
        col = Cx[which(names(Cx)==factor(nchar(cx[i])))],
        type = "l",pch=18,
        cex = 0.6)
}



ularge<-unique(which(ss_replication>100,arr.ind=TRUE)[,2])



fi<-matrix(0,nrow=nrow(N),ncol=1)

fi[w6,1]<-sample(1:500,length(w6))/1000
fi[w6[10],1]<-0.999
LL_replication_1 <- foreach(i = 1:nrep,.combine="c") %dopar% {
  
  
  Pars<-list(stoich=N,
             rates=ratesx,
             total=totV,
             replicator=whorep*replic[i], 
             L=nchar(rownames(N)),
             fitness=fi,
             monomers=c(uno,zero))
  
  tmp   <- ode(yini, times, sistRepl6, Pars,method = "ode45")
  
  print(i)
  cbind(tmp[nrow(tmp),-1])
  
  
}

#stopCluster(my.cluster)
ss_replication_1<-matrix(LL_replication_1,nrow=nrep,ncol=2*totV,byrow = TRUE)
s_replication_1<-replic
cx<-rownames(N)

Cx<-brewer.pal(7,"Dark2")
names(Cx)<-1:max(nchar(cx))
i<-1


ss_replication_1<-ss_replication_1[,-which(nchar(cx)==7)]

cx<-cx[-which(nchar(cx)==7)]
#ss_replication_1[ss_replication_1<1E-6]<-0
plot(s_replication_1,1E-07+(ss_replication_1[,i]/rowSums(ss_replication_1)),
     col=Cx[which(names(Cx)==factor(nchar(cx[i])))],
     ylab="Steady State level",
     xlab="Replication selection strength",
     main="Prevolution to evolution",
     #ylim=c(1E-05,1),
     xlim=c(min(s_replication_1),max(s_replication_1)),#max(ss_replication_1)),
     ylim=c(min((ss_replication_1/rowSums(ss_replication_1))),max(ss_replication_1/rowSums(ss_replication_1))),log="xy",type="l",pch=18)

for(i in 2:ncol(ss_replication_1)){
  lines(s_replication_1,1E-07+(ss_replication_1[,i]/rowSums(ss_replication_1)),
        col = Cx[which(names(Cx)==factor(nchar(cx[i])))],
        type = "l",pch=18,
        cex = 0.6)
}

stopCluster(my.cluster)


