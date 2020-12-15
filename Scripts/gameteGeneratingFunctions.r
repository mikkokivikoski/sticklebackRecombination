pascal.likelihoods=function(n){
	likelihoods=c()
	for(i in 0:n){
	  likelihoods=c(likelihoods,choose(n,i)*0.5^n)
	}
	return(likelihoods)
}

trunc.pois=function(lambda,k){
    p=(lambda^k)/((exp(lambda)-1)*factorial(k))
    p=p/sum(p)
    return(p)
}

generate.gamete.trunc.pois=function(lambda){
    n.cos=sample(x=1:10,size=1,prob=trunc.pois(lambda,1:10))
    obs.co=sample(x=0:n.cos,size=1,prob=pascal.likelihoods(n.cos))
    return(obs.co)
}
generate.gamete.yf=function(probs){
    n.cos=sample(x=0:9,size=1,prob=probs)
    obs.co=sample(x=0:n.cos,size=1,prob=pascal.likelihoods(n.cos))
    return(obs.co)
}

generate.gamete.pois=function(lambda){
    n.cos=1+rpois(1,lambda)
    obs.co=sample(x=0:n.cos,size=1,prob=pascal.likelihoods(n.cos))
    return(obs.co)
}

generate.gamete.geom=function(p){
    n.cos=1+rgeom(1,p)
    obs.co=sample(x=0:n.cos,size=1,prob=pascal.likelihoods(n.cos))
    return(obs.co)
}

generate.gamete.negbin=function(size,prob){
    n.cos=1+rnbinom(n=1,size=size,prob=prob)
    obs.co=sample(x=0:n.cos,size=1,prob=pascal.likelihoods(n.cos))
    return(obs.co)

}
test.fit.trunc.pois=function(iter,OBS,lambda,n.gametes){
    p.vals=c()
    p.vals2=c()
    chi=c()
    for(i in 1:iter){
        D=as.numeric(table(replicate(n=n.gametes,expr=generate.gamete.trunc.pois(lambda))))
        L=max(length(D),length(OBS))
        
        obs=c(OBS,rep(0,L-length(OBS)))
        d=c(D,rep(0,L-length(D)))
		
		test.matrix=matrix(c(d,obs),byrow=T,nrow=2)
		test.result=chisq.test(test.matrix)
		p.vals=c(p.vals,test.result$p.value)
		chi=c(chi,test.result$statistic)
		if(L>3){
			test.matrix2=test.matrix[,1:3]
			if(L>4){
			    test.matrix2[,3]=test.matrix2[,3]+apply(test.matrix[, 4:L],1,sum)
		    } else{
				test.matrix2[,3]=test.matrix2[,3]+test.matrix[,4]
			}
			test.result=chisq.test(test.matrix2)
		}
		p.vals2=c(p.vals2,test.result$p.value)	
    }
    return(list(p.vals,p.vals2,chi))
}
test.fit.geom=function(iter,OBS,p,n.gametes){
    p.vals=c()
    p.vals2=c()
    chi=c()
    for(i in 1:iter){
        D=as.numeric(table(replicate(n=n.gametes,expr=generate.gamete.geom(p))))
        L=max(length(D),length(OBS))
        
        obs=c(OBS,rep(0,L-length(OBS)))
        d=c(D,rep(0,L-length(D)))
		
		test.matrix=matrix(c(d,obs),byrow=T,nrow=2)
		test.result=chisq.test(test.matrix)
		p.vals=c(p.vals,test.result$p.value)
		chi=c(chi,test.result$statistic)
		if(L>3){
			test.matrix2=test.matrix[,1:3]
			if(L>4){
			    test.matrix2[,3]=test.matrix2[,3]+apply(test.matrix[, 4:L],1,sum)
		    } else{
				test.matrix2[,3]=test.matrix2[,3]+test.matrix[,4]
			}
			test.result=chisq.test(test.matrix2)
		}
		p.vals2=c(p.vals2,test.result$p.value)	
    }
    return(list(p.vals,p.vals2,chi))
}

test.fit.negbin=function(iter,OBS,size,p,n.gametes){
    p.vals=c()
    p.vals2=c()
    for(i in 1:iter){
        D=as.numeric(table(replicate(n=n.gametes,expr=generate.gamete.negbin(size,p))))
        L=max(length(D),length(OBS))
        
        obs=c(OBS,rep(0,L-length(OBS)))
        d=c(D,rep(0,L-length(D)))
		
		test.matrix=matrix(c(d,obs),byrow=T,nrow=2)
		test.result=chisq.test(test.matrix)
		p.vals=c(p.vals,test.result$p.value)	
        if(L>3){
			test.matrix2=test.matrix[,1:3]
			if(L>4){
			    test.matrix2[,3]=test.matrix2[,3]+apply(test.matrix[, 4:L],1,sum)
		    } else{
				test.matrix2[,3]=test.matrix2[,3]+test.matrix[,4]
			}
			test.result=chisq.test(test.matrix2)
		}
		p.vals2=c(p.vals2,test.result$p.value)	
    }
    return(p.vals)
}

test.fit.pois=function(iter,OBS,lambda,n.gametes){
    p.vals=c()
    p.vals2=c()
    chi=c()
    for(i in 1:iter){
        D=as.numeric(table(replicate(n=n.gametes,expr=generate.gamete.pois(lambda))))
        L=max(length(D),length(OBS))
        
        obs=c(OBS,rep(0,L-length(OBS)))
        d=c(D,rep(0,L-length(D)))
		
		test.matrix=matrix(c(d,obs),byrow=T,nrow=2)
		test.result=chisq.test(test.matrix)
		p.vals=c(p.vals,test.result$p.value)
		chi=c(chi,test.result$statistic)
		
		if(L>3){
			test.matrix2=test.matrix[,1:3]
			if(L>4){
			    test.matrix2[,3]=test.matrix2[,3]+apply(test.matrix[, 4:L],1,sum)
		    } else{
				test.matrix2[,3]=test.matrix2[,3]+test.matrix[,4]
			}
			test.result=chisq.test(test.matrix2)
		}
		p.vals2=c(p.vals2,test.result$p.value)	
    }
    return(list(p.vals,p.vals2,chi))
}

test.fit.yf=function(iter,OBS,probs,n.gametes){
    p.vals=c()
    p.vals2=c()
    for(i in 1:iter){
        D=as.numeric(table(replicate(n=n.gametes,expr=generate.gamete.yf(probs))))
        L=max(length(D),length(OBS))
        
        obs=c(OBS,rep(0,L-length(OBS)))
        d=c(D,rep(0,L-length(D)))
		
		test.matrix=matrix(c(d,obs),byrow=T,nrow=2)
		test.result=chisq.test(test.matrix)
		p.vals=c(p.vals,test.result$p.value)
		
		if(L>3){
			test.matrix2=test.matrix[,1:3]
			if(L>4){
			    test.matrix2[,3]=test.matrix2[,3]+apply(test.matrix[, 4:L],1,sum)
		    } else{
				test.matrix2[,3]=test.matrix2[,3]+test.matrix[,4]
			}
			test.result=chisq.test(test.matrix2)
		}
		p.vals2=c(p.vals2,test.result$p.value)	
    }
    return(list(p.vals,p.vals2))
}

