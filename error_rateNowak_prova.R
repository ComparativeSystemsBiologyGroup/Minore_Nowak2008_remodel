sistRepl6error<-function(Time,State,Pars)
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
    sumfix<-sum(fitness*State)
    sumx<-sum(State)
    FIT<-fitness
    if(sumx==0){
      vr<-matrix(0,nrow=ncol(stoich),ncol=1) 
    }else{
      w<-which(stoich>0,arr.ind=TRUE)
      replica <- State[w[,1]] * replicator[w[,1]] * (FIT[w[,1]] - sumfix/sumx)
      trg<-which(w[,1]==target)
        #replica[w[trg,1]]  <- State[w[trg,1]] * replicator[w[trg,1]] * FIT[w[trg,1]] * (1-error_prob)^6
        - State[w[trg,1]] * replicator[w[trg,1]] *  sumfix/sumx 
        - State[w[trg,1]] * replicator [w[trg,1]]* FIT[w[trg,1]] * (1-(1-error_prob)^6)
      
      vr<-matrix(0,nrow=ncol(stoich),ncol=1)      
      vr[w[,2]]<-replica
      vr[which(stoich[target,]>0)]<-State[target] * replicator[target] * FIT[target] * (1-error_prob)^6
      - State[target] * replicator[target] *  sumfix/sumx 
      - State[target] * replicator [target]* FIT * (1-(1-error_prob)^6)
    }
    x <- stoich %*% (v+vr)
    x[which(x<0)]<-0
    list(x)
  })
}


w6<-which(nchar(rownames(N))==3)
fi<-matrix(0,nrow=nrow(N),ncol=1)


target<-1
fi[w6[target],1]<- 1

allep<-seq(from=log10(0.001),to=log10(0.9),by=.1)
nrep<-length(allep)
repstrength<-3
whorep<-matrix(repstrength,nrow=nrow(N),ncol=1)
cx<-rownames(N)
#####
times <- seq(0, 20,1)





LL_replicationerror <- foreach(i = 1:nrep,.combine="c") %dopar% {
  Pars<-list(stoich=N,
             rates=rates,
             total=totV,
             replicator=whorep, 
             L=nchar(rownames(N)),
             fitness=fi,
             monomers=c(uno,zero),
             error_prob=10^allep[i],
             target=w6[target])
  
  tmp   <- ode(yini, times, sistRepl6error, Pars,method = "ode45")
  
  print(i)
  cbind(tmp[nrow(tmp),-1])
}

ss_replicationerror<-matrix(LL_replicationerror,nrow=nrep,ncol=2*totV,byrow = TRUE)
s_replicationerror<-allep
cx<-rownames(N)
Cx<-brewer.pal(10,"Set3")
names(Cx)<-1:max(nchar(cx))
i<-1

#ss_replicationerror<-ss_replicationerror[,-which(nchar(cx)==7)]
#cx<-cx[-which(nchar(cx)==7)]

RRR<-ss_replicationerror/rowSums(ss_replicationerror)

plot(10^s_replicationerror,(ss_replicationerror[,i]/rowSums(ss_replicationerror)),
     col=Cx[which(names(Cx)==factor(nchar(cx[i])))],
     ylab="Steady State level",
     xlab="Error probability",
     main="Error critical threshold",
     #ylim=c(1E-05,1),
     ylim=c(min(RRR[RRR>0]),max(RRR)),
     xlim=c(min(10^s_replicationerror),
            max(10^s_replicationerror)),
     log="xy",
     type="l",
     pch=18,lwd=2)

for(i in 2:ncol(ss_replicationerror)){
  lines(10^s_replicationerror,(ss_replicationerror[,i]/rowSums(ss_replicationerror)),
        col = Cx[which(names(Cx)==factor(nchar(cx[i])))],
        type = "l",pch=18,
        cex = 0.6,
        lwd=2)
}





# RRR<-1E-07+ss_replicationerror
# 
# 
# plot(10^s_replicationerror,1E-07+(ss_replicationerror[,i]),
#      col=Cx[which(names(Cx)==factor(nchar(cx[i])))],
#      ylab="Steady State level",
#      xlab="Error probability",
#      main="Error critical threshold",
#      #ylim=c(1E-05,1),
#      ylim=c(min(RRR[RRR>0]),max(RRR)),
#      xlim=c(min(10^s_replicationerror),
#             max(10^s_replicationerror)),
#      log="xy",
#      type="l",
#      pch=18,lwd=2)
# 
# for(i in 2:ncol(ss_replicationerror)){
#   lines(10^s_replicationerror,1E-07+(ss_replicationerror[,i]),
#         col = Cx[which(names(Cx)==factor(nchar(cx[i])))],
#         type = "l",pch=18,
#         cex = 0.6)
# }
# 





uc<-1-((gamma+2*ai0)/repstrength)^1/6
