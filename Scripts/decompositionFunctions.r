ASSIGN.DATA.IN.CLASSES.CLEANING.FUNCTION=function(DATA.VECTOR.MATRIX,PRIOR.DISTRIBUTIONS,PROBABILITY.TABLE,ITERATIONS,ADDING.UNIT,PLOT,PRINT,PLOT.EVERY,LG,SEX){	
	for(I in 1:ITERATIONS){
	    if(I == 1){
	        posterior.matrix=matrix(ncol=ncol(DATA.VECTOR.MATRIX),nrow=10,0)
            f.priors=PRIOR.DISTRIBUTIONS
            if(PLOT){
	            plot('n',xlim=c(1,ncol(posterior.matrix)),ylim=c(0,1),main=paste0(LG,SEX,"priors"),ylab="probability",xlab="Chromosome bin")
	            COLORS=c("gray","blue","red","green3","brown","orange","yellow1","yellow2","yellow3","yellow4")
	            for(i in 1:nrow(posterior.matrix)){
					lines(f.priors[i,],col=COLORS[i])
				}
             } 
          } else {
			 f.priors=posterior.matrix
		 }
		for(i in 1:nrow(DATA.VECTOR.MATRIX)){
			data.vector=DATA.VECTOR.MATRIX[i,]
			co.class=sum(data.vector)
			class.priors=as.numeric(PROBABILITY.TABLE[PROBABILITY.TABLE$OBSERVED == co.class,2:11])
			
			class.posteriors=c()
			
			for(j in 0:9){
				if(class.priors[j+1]==0){
					class.posteriors=c(class.posteriors,0)
				} else{
					class.posteriors=c(class.posteriors,(class.priors[j+1]*prod(f.priors[j+1,]^(data.vector)))/sum(class.priors*t(apply(f.priors,1,function(x) prod(x^(data.vector))))))
			    }
		    }
			assigned.label=sample(seq(0,9),prob=class.posteriors,1)	
			posterior.matrix[assigned.label+1,]=posterior.matrix[assigned.label+1,]+data.vector*ADDING.UNIT	
	    }   
	for(i in 1:nrow(posterior.matrix)){
		 if(sum(posterior.matrix[i,])>0){
	         posterior.matrix[i,]=posterior.matrix[i,]/sum(posterior.matrix[i,])
		 }
	}
    if(I%%10==0 & PLOT){
         plot('n',xlim=c(1,ncol(posterior.matrix)),ylim=c(0,1),main=paste0(LG,SEX," After iteration ",I),ylab="probability",xlab="Chromosome bin")
         COLORS=c("gray","blue","red","green3","brown","orange","yellow1","yellow2","yellow3","yellow4")
         for(i in 1:nrow(posterior.matrix)){
   		    if(sum(PROBABILITY.TABLE[,i+1])>0){
				lines(posterior.matrix[i,],col=COLORS[i])
			}
		}
	}
}
	return(posterior.matrix)	
}

main.ASSIGN.DATA.IN.CLASSES.CLEANING.FUNCTION=function(d,LG,LG.length,MAX.COUNT,MLE.ESTIMATES.TABLE,NO.BINS,interval.length,SEED,ITERATIONS,ADDING.UNIT,SEX,PLOT){
	probability.table=GENERATE.PROBABILITY.TABLE(MLE.ESTIMATES.TABLE,LG,MAX.COUNT)
	data.vectors.as.matrix = CONSTRUCT.DATA.VECTOR.MATRIX(d,NO.BINS,MAX.COUNT)
	prior.distributions=GENERATE.PRIORS(NO.BINS,RANDOM=F,SEED)
	cleaned.distributions=ASSIGN.DATA.IN.CLASSES.CLEANING.FUNCTION(data.vectors.as.matrix,prior.distributions,probability.table,ITERATIONS,ADDING.UNIT,PLOT,FALSE,500,LG,SEX)
	return(cleaned.distributions)
}
####
PREPARE.DATA=function(d,SEX,INTERVAL.LENGTH){
		d$SITE1 = floor(d[,paste0(SEX,"SITE1")]/INTERVAL.LENGTH)+1
		d$SITE2 = floor(d[,paste0(SEX,"SITE2")]/INTERVAL.LENGTH)+1
		d$SITE3 = floor(d[,paste0(SEX,"SITE3")]/INTERVAL.LENGTH)+1
		d$SITE4 = floor(d[,paste0(SEX,"SITE4")]/INTERVAL.LENGTH)+1
		d$SITE5 = floor(d[,paste0(SEX,"SITE5")]/INTERVAL.LENGTH)+1
	return(d)
}

CONSTRUCT.DATA.VECTOR.MATRIX=function(DATA,No.BINS,MAX.SITES){
	data.vectors.as.matrix=matrix(nrow=nrow(DATA),ncol=No.BINS,0)
    for(i in 1:nrow(DATA)){
        sites = as.numeric(DATA[i,paste0("SITE",seq(1,MAX.SITES))])
	    sites=sites[is.na(sites)==FALSE]
	    if(length(sites)>0){
	        for(s in sites){
		        data.vectors.as.matrix[i,s] = data.vectors.as.matrix[i,s]+1 
	    }
	        }
    }
    return(data.vectors.as.matrix)
}

GENERATE.BINS = function(LG.LENGTH,BY.LENGTH,EVEN.LENGTH,BY.NUMBER,LENGTH,NUMBER){
		    if(BY.NUMBER & EVEN.LENGTH){
		        INTERVAL.LENGTH = (1+0.001/NUMBER)*(LG.LENGTH/NUMBER)
		        No.BINS=NUMBER       
		    } else if(BY.LENGTH & EVEN.LENGTH){
		        if(LG.LENGTH%%LENGTH >= 0.5*LENGTH){
		            No.BINS=length(seq(0,LG.LENGTH,LENGTH))
		            INTERVAL.LENGTH=(1+0.001/No.BINS)*(seq(0,LG.LENGTH,length.out=No.BINS+1)[2]-seq(0,LG.LENGTH,length.out=No.BINS+1)[1])
		        } else {
		            No.BINS=length(seq(0,LG.LENGTH,LENGTH))-1
		            INTERVAL.LENGTH=(1+0.001/No.BINS)*(seq(0,LG.LENGTH,length.out=No.BINS+1)[2]-seq(0,LG.LENGTH,length.out=No.BINS+1)[1])
		        }
		    } else if(BY.LENGTH){
		        INTERVAL.LENGTH=LENGTH
		        No.BINS=floor(LG.LENGTH/INTERVAL.LENGTH)+1
		    }
		    return(c(No.BINS,INTERVAL.LENGTH))
}

GENERATE.PROBABILITY.TABLE=function(MLE.ESTIMATE.TABLE,LG,MAX.COUNT){
	PROBABILITY.TABLE=data.frame(OBSERVED=seq(0,5),CO0=0,CO1=0,CO2=0,CO3=0,CO4=0,CO5=0,CO6=0,CO7=0,CO8=0,CO9=0)
	for(i in 0:MAX.COUNT){
		PROBABILITY.TABLE[i+1,2:11] = choose(seq(0,9),i)*(0.5^seq(0,9))*as.numeric(MLE.ESTIMATE.TABLE[MLE.ESTIMATE.TABLE$CHR==LG,paste0("MLEp",seq(0,9))])/sum(choose(seq(0,9),i)*(0.5^seq(0,9))*as.numeric(MLE.ESTIMATE.TABLE[MLE.ESTIMATE.TABLE$CHR==LG,paste0("MLEp",seq(0,9))]))
	}
	PROBABILITY.TABLE[PROBABILITY.TABLE$OBSERVED>MAX.COUNT,2:11]=0
	return(PROBABILITY.TABLE)
}

GENERATE.PRIORS=function(No.BINS,RANDOM,SEED){
    p=t(replicate(rep(1/No.BINS,No.BINS),n=10))
	if(RANDOM){
		set.seed(SEED)
	    p=t(replicate(runif(No.BINS),n=10))
	    p=p/rowSums(p)
    }
    return(p)
}

EM.LIKE.CLEANING.FUNCTION=function(DATA.VECTOR.MATRIX,PRIOR.DISTRIBUTIONS,PROBABILITY.TABLE,ITERATIONS,ADDING.UNIT,PLOT,PRINT,PLOT.EVERY,LG,SEX){

	for(I in 1:ITERATIONS){
		if(I==1){
			f.priors=PRIOR.DISTRIBUTIONS
			f.priors.tmp = f.priors
			if(PLOT){
				plot('n',xlim=c(1,ncol(f.priors)),ylim=c(0,1),main=paste0(LG,SEX,"priors"),ylab="probability",xlab="Chromosome bin")
	            COLORS=c("gray","blue","red","green3","brown","orange","yellow1","yellow2","yellow3","yellow4")
	            for(i in 1:nrow(f.priors)){
					lines(f.priors[i,],col=COLORS[i])
				}
		    }
		}
		for(i in 1:nrow(DATA.VECTOR.MATRIX)){
		    data.vector=DATA.VECTOR.MATRIX[i,]
		    co.class=sum(data.vector)
		    class.priors=as.numeric(PROBABILITY.TABLE[PROBABILITY.TABLE$OBSERVED == co.class,2:11])			
			class.posteriors=c()
			for(j in 0:9){
				if(class.priors[j+1]==0){
					class.posteriors=c(class.posteriors,0)
				} else{
					class.posteriors=c(class.posteriors,(class.priors[j+1]*prod(f.priors[j+1,]^(data.vector)))/sum(class.priors*t(apply(f.priors,1,function(x) prod(x^(data.vector))))))
			}
		}
			f.priors.tmp = f.priors.tmp + ADDING.UNIT*sapply(data.vector,function(x) x*class.posteriors)
			f.priors.tmp = f.priors.tmp/rowSums(f.priors.tmp)		    
	    }
	    
	    f.priors=f.priors.tmp

	    if(I%%10==0 & PLOT){
            plot('n',xlim=c(1,ncol(f.priors)),ylim=c(0,1),main=paste0(LG,SEX," After iteration ",I),ylab="probability",xlab="Chromosome bin")
            COLORS=c("gray","blue","red","green3","brown","orange","yellow1","yellow2","yellow3","yellow4")
            for(i in 2:nrow(f.priors)){
				if(sum(PROBABILITY.TABLE[,i+1])>0){
					lines(f.priors[i,],col=COLORS[i])
				}
			}
		}	    
		}
		zero.ditributions=which(colSums(PROBABILITY.TABLE[,-1])==0)
		f.priors[zero.ditributions,]=0
    return(f.priors)

}

main.EM.LIKE.CLEANING.FUNCTION=function(DATA,META.DATA,LG,NO.BINS,SEED,ITERATIONS,PLOT){
	 d=subset(DATA,CHR==LG)
	 l=list()
	 for(SEX in c("PATERNAL","MATERNAL")){
		print(paste(LG,SEX))
        max.count=max(d[,paste0(SEX,"COUNT")])
        mle.estimates=META.DATA[,c("CHR",paste0(SEX,"n",0:5),paste0(SEX,"MLEp",0:9))]
		colnames(mle.estimates) = gsub(pattern=SEX, replacement="", colnames(mle.estimates))
		LG.length = META.DATA[META.DATA$CHR==LG,"LENGTH"]
        
        interval.length=GENERATE.BINS(LG.LENGTH=LG.length,BY.LENGTH=FALSE,EVEN.LENGTH=TRUE,BY.NUMBER=TRUE,LENGTH=2000000,NUMBER=NO.BINS)[2]
	    no.bins=GENERATE.BINS(LG.LENGTH=LG.length,BY.LENGTH=FALSE,EVEN.LENGTH=TRUE,BY.NUMBER=TRUE,LENGTH=2000000,NUMBER=NO.BINS)[1]
        
        prepared.data = PREPARE.DATA(d,SEX,interval.length)
        
		ADDING.UNIT=1/938
		YLIM=c(0,0.7)
		COL=c("blue","red","green3","brown","orange")
		probability.table=GENERATE.PROBABILITY.TABLE(mle.estimates,LG,max.count)
		data.vectors.as.matrix = CONSTRUCT.DATA.VECTOR.MATRIX(prepared.data,no.bins,max.count)
		
		prior.distributions=GENERATE.PRIORS(no.bins,RANDOM=F,SEED)
    	cleaned.distributions=EM.LIKE.CLEANING.FUNCTION(data.vectors.as.matrix,prior.distributions,probability.table,ITERATIONS,ADDING.UNIT,PLOT,FALSE,500,LG,SEX)
    	l[[SEX]]=cleaned.distributions
	}
	return(l)
}

main2=function(DATA,META.DATA,LG,NO.BINS,SEED,ITERATIONS,PLOT){
	 d=subset(DATA,CHR==LG)
	 l=list()
	 for(SEX in c("PATERNAL","MATERNAL")){
		print(paste(LG,SEX))
        max.count=max(d[,paste0(SEX,"COUNT")])
        mle.estimates=META.DATA[,c("CHR",paste0(SEX,"n",0:5),paste0(SEX,"MLEp",0:9))]
		colnames(mle.estimates) = gsub(pattern=SEX, replacement="", colnames(mle.estimates))
		LG.length = META.DATA[META.DATA$CHR==LG,"LENGTH"]
        interval.length=GENERATE.BINS(LG.LENGTH=LG.length,BY.LENGTH=FALSE,EVEN.LENGTH=TRUE,BY.NUMBER=TRUE,LENGTH=2000000,NUMBER=NO.BINS)[2]
	    no.bins=GENERATE.BINS(LG.LENGTH=LG.length,BY.LENGTH=FALSE,EVEN.LENGTH=TRUE,BY.NUMBER=TRUE,LENGTH=2000000,NUMBER=NO.BINS)[1]
        
        prepared.data = PREPARE.DATA(d,SEX,interval.length)
        
		
		ADDING.UNIT=1/938
		YLIM=c(0,0.7)
		COL=c("blue","red","green3","brown","orange")	
		l[[SEX]]=main.ASSIGN.DATA.IN.CLASSES.CLEANING.FUNCTION(prepared.data,LG,LG.length,max.count,mle.estimates,no.bins,interval.length,SEED,ITERATIONS,ADDING.UNIT,SEX,PLOT)
	}
	return(l)
}

#########
assign.data=function(D,PRIOR.DISTRIBUTIONS,PROBABILITY.TABLE,SEX,INTERVAL.LENGTH,No.BINS,MAX.SITES){
	f.priors=PRIOR.DISTRIBUTIONS
	d=PREPARE.DATA(D,SEX,INTERVAL.LENGTH)
	d[,paste0("TRUECOCOUNT",SEX)]=100
	data.vector.matrix=CONSTRUCT.DATA.VECTOR.MATRIX(d,No.BINS,MAX.SITES)
    for(i in 1:nrow(data.vector.matrix)){
		data.vector=data.vector.matrix[i,]
		co.class=sum(data.vector)
		class.priors=as.numeric(PROBABILITY.TABLE[PROBABILITY.TABLE$OBSERVED == co.class,2:11])
		class.posteriors=c()
			for(j in 0:9){
				if(class.priors[j+1]==0){
					class.posteriors=c(class.posteriors,0)
				} else{
					class.posteriors=c(class.posteriors,(class.priors[j+1]*prod(f.priors[j+1,]^(data.vector)))/sum(class.priors*t(apply(f.priors,1,function(x) prod(x^(data.vector))))))
			    }
		    }
		d[i,paste0("TRUECOCOUNT",SEX)]=sample(seq(0,9),prob=class.posteriors,1)
	}
	return(d)
}
