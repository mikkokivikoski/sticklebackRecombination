prepare.prior=function(DATA,RANDOM){
    N=length(DATA)
	if(RANDOM){
	    p=runif(2*(N-1))+0.001
	    p=p/sum(p)
    }
    return(p)
}
MLEforP = function(DATA,ITER,PRINT,RESTRICTED,SET.P.PRIOR,P.PRIOR){
	L=c()
	n=DATA
	Co.class.in.focus=seq(0,((2*(length(DATA)-1))-1))
	N=length(DATA)
	ITER=ITER
	p=rep(round(1/(2*(N-1)),4),2*(N-1)-1)
    p=c(p,1-sum(p))
	if(RESTRICTED){
        p=rep(round(1/(2*(N-1)-1),4),2*(N-1)-1-1)
        p=c(0,p,1-sum(p))
    }
    if(SET.P.PRIOR){
	    p=P.PRIOR
	}
	Q=c()
	for(k in 1:ITER){
	sumI2N = c()
	for(i in 1:length(Co.class.in.focus)){
		sumJN=c()
		for(j in 1:N){  
		  sumS2N1 = c()  
		  for(s in 1:length(Co.class.in.focus)){
		      sumS2N1 = c(sumS2N1,((choose((s-1),(j-1)))*((0.5)^(s-1))*p[s]))
		  }
		  sumJN = c(sumJN,(((choose((i-1),(j-1)))*((0.5)^(i-1))*p[i])/sum(sumS2N1))*n[j]) 
		}
		sumI2N = c(sumI2N,sum(sumJN))   
	}
	
	E = sumI2N
	if(k%%1000==0){	
	    Q=c(Q,sum(E[which(p>0)]*log(p[which(p>0)])))
	}
	    
	for(i in 1:length(p)){
	    p[i] = E[i]/sum(n)
	}
	if(PRINT){
	if(k==1 || k==10 || k==100 ||k==ITER){
		print(p)
		plot(p[p>0],xlim=c(1,length(Co.class.in.focus)),ylim=c(0,1))
		points(p[p==0],col="red",x=which(p==0))
	}
}   
    L=c(L,calculate.likelihood.Function(generate.Qp(n,p),n,TRUE)) #Store likelihood for each iteration of p
	}
	return(p)
}
generate.Qp = function(n,p){
	K=2*(length(n)-1)-1+1
    qp = c()
    for(i in 1:length(p)){
		sumJK=c()
		for(j in 1:K){
			sumJK=c(sumJK,((choose((j-1),(i-1))*0.5^(j-1))*p[j]))
		}
		qp=c(qp,sum(sumJK))
	}
	return(qp)
}
generateDataFrom.Qp = function(M,r,Qp,DATA){
	R=data.frame()
	for(i in 1:r){
		d=as.vector(table(sample(seq(1,length(Qp)),M,prob=Qp,replace=T)))
		R=rbind(R,c(d,rep(0,(2*(length(DATA)-1)-length(d)))))
	}
	return(R)
}
MLE.From.Genrated.Data = function(R,ITER){
    All.pjs=data.frame()
    for(j in 1:nrow(R)){
	    d=as.numeric(R[j,][R[j,]>0])
		Nj=length(d)
	    p.start = prepare.prior(d,TRUE)
		pj=MLEforP(d,ITER,FALSE,FALSE,TRUE,p.start) 
		All.pjs=rbind(All.pjs,pj)
	}
	return(All.pjs)	
}
Find.Conf.Intervals=function(All.pjs,ALPHA){
    Conf.Intervals = data.frame("lowerBound"=rep(100,ncol(All.pjs)),"upperBound"=rep(100,ncol(All.pjs)))
    for(k in 1:ncol(All.pjs)){
	    Lk=quantile(All.pjs[,k],probs=c(ALPHA,1-ALPHA))[1]
        Uk=quantile(All.pjs[,k],probs=c(ALPHA,1-ALPHA))[2]
        Conf.Intervals[k,]=c(as.numeric(Lk),as.numeric(Uk))
    }
    return(Conf.Intervals)
}
calculate.likelihood.Function = function(q,data,LOG){
    L=prod((q[1:length(data)])^data)
    if(LOG){
	    l=c()
	    for(j in 1:length(data)){
	     l=c(l,(data[j]*log(q[j])))
	    }
	    L=sum(l)
	}
    return(L)
}
calculate.delta.n=function(DELTA.N, GENERATED.DATA,RESTRICTED.Q,NON.RESTRICTED.Q,LOG){
    delta.ns = c()
    L.re=c()
    L.non.re=c()
    for(i in 1:nrow(GENERATED.DATA)){
		d=as.numeric(GENERATED.DATA[i,])
		d=d[d>0]
		DELTA.n = 2*calculate.likelihood.Function(NON.RESTRICTED.Q,d,LOG) - 2*calculate.likelihood.Function(RESTRICTED.Q,d,LOG)
		delta.ns = c(delta.ns,DELTA.n)
		L.re=c(L.re,calculate.likelihood.Function(RESTRICTED.Q,d,LOG))
		L.non.re=c(L.non.re,calculate.likelihood.Function(NON.RESTRICTED.Q,d,LOG))
	}
	return(delta.ns)
}	
calculate.p.value=function(DELTA.N.ORIG, DELTA.N.SIMULATED){
    p.val=length(DELTA.N.SIMULATED[DELTA.N.SIMULATED>=DELTA.N.ORIG])/length(DELTA.N.SIMULATED)
    return(p.val)
}


GET.ps = function(DATA,ITER){
	p.unrestricted.1=MLEforP(DATA,ITER,FALSE,FALSE,T,prepare.prior(DATA,T))
	p.unrestricted.2=MLEforP(DATA,ITER,FALSE,FALSE,T,prepare.prior(DATA,T))
	p.unrestricted.3=MLEforP(DATA,ITER,FALSE,FALSE,T,prepare.prior(DATA,T))
	p.restricted=MLEforP(n,ITER,FALSE,TRUE,FALSE,3)
	
	return(list(p.unrestricted.1=p.unrestricted.1,p.unrestricted.2=p.unrestricted.2,p.unrestricted.3=p.unrestricted.3,p.restricted=p.restricted))

}

concatenate.data=function(DATA.COLUMN){
	
    tmp=table(DATA.COLUMN)
    data=c()
 	max.class=max(as.numeric(names(tmp)))
    for(i in 0:max.class){
 	   if(as.character(i) %in% names(tmp)){
		  data=c(data,as.numeric(tmp[as.character(i)]))
	   } else{
		    data=c(data,0)
	   }
    }
    return(data)
}
