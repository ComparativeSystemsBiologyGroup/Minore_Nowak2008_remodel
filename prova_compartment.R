

#####insert replication in the model
#1


sistweighted<-function(Time,State,Pars)
{
  State[which(State < 0)] <- 0
  with(as.list(c(State,Pars)),{
    w<-which(stoich<0,arr.ind=TRUE)    
    v <- State[w[,1]]*rates[w[,2]]
    posUno<-which(rownames(stoich)=="1")
    posZero<-which(rownames(stoich)=="0")
    oioi<-c(monomers[1] * rates[which(stoich[posUno,]>0)],monomers[2] * rates[which(stoich[posZero,]>0)])
    names(v)<-NULL
    v<-c(v[1:(which(stoich[posUno,]>0)-1)],oioi,v[(which(stoich[posUno,]>0)):length(v)])
    #replication rates
    sumfix<-sum(fitness*State*(1+(State[target]/sum(State))*weights))
    sumx<-sum(State)
    w<-which(stoich>0,arr.ind=TRUE)    
    if(sumx==0){
      vr<-matrix(0,nrow=ncol(stoich),ncol=1) 
    }else{
      W<-1+(State[target]/sum(State))*weights[w[,1]]
      replica <- State[w[,1]] * replicator[w[,1]] * (fitness[w[,1]] * W - sumfix/sumx)
      vr<-matrix(0,nrow=ncol(stoich),ncol=1)      
      vr[w[,2]]<-replica
    }
    x <- stoich %*% (v+vr)
    x[which(x<0)]<-0
    list(x)
  })
}


whorep<-matrix(1,nrow=nrow(N),ncol=1)

cx<-rownames(N)
#####
times <- seq(0, 50,1)

ww<-10^seq(log10(7),log10(50),.01)
w6<-which(nchar(rownames(N))==6)
fi<-matrix(0,nrow=nrow(N),ncol=1)

#fi[w6,1]<-runif(n = length(w6),min = 0,max=1)
fi[w6,1]<-.999
LL_replication_weight <- foreach(i = 1:length(ww),.combine="c") %dopar% {
  weightx<-matrix(0,nrow=nrow(N),ncol=1)
  weightx[w6]<-ww[i]
  weightx[w6[10]] <- 10
  Pars<-list(stoich=N,
             rates=rates,
             total=totV,
             replicator=whorep*3, 
             L=nchar(rownames(N)),
             fitness=fi,
             monomers=c(uno,zero),
             target=w6[10],
             weights=weightx)
  
  tmp   <- ode(yini, times, sistweighted, Pars,method = "ode45")
  #plot.deSolve(tmp)
  print(i)
  cbind(tmp[nrow(tmp),-1])
}

ss_replication<-matrix(LL_replication_weight,nrow=length(ww),ncol=2*totV,byrow = TRUE)
s_replication<-ww
cx<-rownames(N)
Cx<-brewer.pal(7,"Dark2")
names(Cx)<-1:max(nchar(cx))
i<-1

ss_replication<-ss_replication[,-which(nchar(cx)==7)]
cx<-cx[-which(nchar(cx)==7)]
pdf(file="Catalytic_frep1_wrep10.pdf")
par(mfrow=c(1,1))
yy<-(ss_replication/rowSums(ss_replication))
miny<-min(yy[yy>0])
plot(s_replication,(ss_replication[,i]/rowSums(ss_replication)),
     col=Cx[which(names(Cx)==factor(nchar(cx[i])))],
     ylab="Steady State level",
     xlab="weight replicator",
     main="Catalytic f=1, w=10",
     #ylim=c(1E-05,1),
     ylim=c(miny,max(ss_replication/rowSums(ss_replication))),
     xlim=c(min(s_replication),
            max(s_replication)),
     log="xy",
     type="l",
     pch=18)

for(i in 2:ncol(ss_replication)){
  lines(s_replication,(ss_replication[,i]/rowSums(ss_replication)),
        col = Cx[which(names(Cx)==factor(nchar(cx[i])))],
        type = "l",pch=18,
        cex = 0.6)
}

lines(s_replication,(ss_replication[,w6[10]]/rowSums(ss_replication)),
      col = "red",
      type = "p",pch=18,
      cex = 1.25)

dev.off()







fi<-matrix(0,nrow=nrow(N),ncol=1)

fi[w6,1]<-runif(n = length(w6),min = 0,max=1)
#fi[w6,1]<-.999
LL_replication_weight <- foreach(i = 1:length(ww),.combine="c") %dopar% {
  weightx<-matrix(0,nrow=nrow(N),ncol=1)
  weightx[w6]<-ww[i]
  weightx[w6[10]] <- 10
  Pars<-list(stoich=N,
             rates=rates,
             total=totV,
             replicator=whorep*3, 
             L=nchar(rownames(N)),
             fitness=fi,
             monomers=c(uno,zero),
             target=w6[10],
             weights=weightx)
  
  tmp   <- ode(yini, times, sistweighted, Pars,method = "ode45")
  #plot.deSolve(tmp)
  print(i)
  cbind(tmp[nrow(tmp),-1])
}

ss_replication<-matrix(LL_replication_weight,nrow=length(ww),ncol=2*totV,byrow = TRUE)
s_replication<-ww
cx<-rownames(N)
Cx<-brewer.pal(7,"Dark2")
names(Cx)<-1:max(nchar(cx))
i<-1

ss_replication<-ss_replication[,-which(nchar(cx)==7)]
cx<-cx[-which(nchar(cx)==7)]
pdf(file="Catalytic_f6_in_01_wrep10.pdf")
par(mfrow=c(1,1))
plot(s_replication,(ss_replication[,i]/rowSums(ss_replication)),
     col=Cx[which(names(Cx)==factor(nchar(cx[i])))],
     ylab="Steady State level",
     xlab="weight replicator",
     main="Catalytic, f6 in [0,1], wrep=10",
     #ylim=c(1E-05,1),
     ylim=c(min(ss_replication[ss_replication>0]),max(ss_replication/rowSums(ss_replication))),
     xlim=c(min(s_replication),
            max(s_replication)),
     log="xy",
     type="l",
     pch=18)

for(i in 2:ncol(ss_replication)){
  lines(s_replication,(ss_replication[,i]/rowSums(ss_replication)),
        col = Cx[which(names(Cx)==factor(nchar(cx[i])))],
        type = "l",pch=18,
        cex = 0.6)
}

lines(s_replication,(ss_replication[,w6[10]]/rowSums(ss_replication)),
      col = "red",
      type = "p",pch=18,
      cex = 1.25)

dev.off()








fifi<-seq(from=0.3,to=1,by=0.03)

fi<-matrix(0,nrow=nrow(N),ncol=1)

fi[w6,1]<-runif(n = length(w6),min = 0,max=1)

weightx<-matrix(0,nrow=nrow(N),ncol=1)
  weightx[w6]<-runif(min=0,max=10,n=length(w6))
  weightx[w6[10]] <- 10
LL_replication_weight <- foreach(i = 1:length(fifi),.combine="c") %dopar% {
  fi[w6[10],1]<-fifi[i]
  
  Pars<-list(stoich=N,
             rates=rates,
             total=totV,
             replicator=whorep*3, 
             L=nchar(rownames(N)),
             fitness=fi,
             monomers=c(uno,zero),
             target=w6[10],
             weights=weightx)
  
  tmp   <- ode(yini, times, sistweighted, Pars,method = "ode45")
  #plot.deSolve(tmp)
  print(i)
  cbind(tmp[nrow(tmp),-1])
}

ss_replication<-matrix(LL_replication_weight,nrow=length(fifi),ncol=2*totV,byrow = TRUE)
s_replication<-ww
cx<-rownames(N)
Cx<-brewer.pal(7,"Dark2")
names(Cx)<-1:max(nchar(cx))
i<-1

ss_replication<-ss_replication[,-which(nchar(cx)==7)]
cx<-cx[-which(nchar(cx)==7)]

par(mfrow=c(1,1))

pdf(file="Catalytic_fitness_rep_scan_wrep10_w6_in_0_5.pdf")
yy<-(ss_replication/rowSums(ss_replication))
miny<-min(yy[yy>0])
xx<-fifi
plot(xx,(ss_replication[,i]/rowSums(ss_replication)),
     col=Cx[which(names(Cx)==factor(nchar(cx[i])))],
     ylab="Steady State level",
     xlab="weight replicator",
     main="Catalytic fitness scan, wrep10",
     #ylim=c(1E-05,1),
     ylim=c(miny,max(ss_replication/rowSums(ss_replication))),
     xlim=c(min(xx),
            max(xx)),
     log="y",
     type="l",
     pch=18)

for(i in 2:ncol(ss_replication)){
  lines(xx,(ss_replication[,i]/rowSums(ss_replication)),
        col = Cx[which(names(Cx)==factor(nchar(cx[i])))],
        type = "l",pch=18,
        cex = 0.6)
}

lines(xx,(ss_replication[,w6[10]]/rowSums(ss_replication)),
      col = "red",
      type = "p",pch=18,
      cex = 1.25)

dev.off()
