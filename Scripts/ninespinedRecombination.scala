#!/usr/bin/env anduril

import anduril.builtin._
import anduril.tools._
import anduril.anima._
import org.anduril.runtime._
//$OPT --threads 14


object ninepinedRecombination {
  val lgs = Seq("LG1","LG2","LG3","LG4","LG5","LG6","LG7","LG8","LG9","LG10","LG11","LG12","LG13","LG14","LG15","LG16","LG17","LG18","LG19","LG20","LG21")
  val recombinationSitesHelsinki = INPUT("")//Ppungitius_Helsinki_linkagemap.txt.gz 
  val ninespinedCentromeres = INPUT("")//Centromere positions according to P. pungitius reference genome ver. 7 (Kivikoski et al. 2020)
  val v7Index = INPUT("")//Index (.fai) of the P. pungitius reference genome ver. 7 (Kivikoski et al. 2020) 
  val emAlgorithmFunctions = INPUT("emAlgortihmForCoProbabilities.r")
  val decompositionFunctions = INPUT("decompositionFunctions.r")
  val gameteGeneratingFunctions = INPUT("gameteGeneratingFunctions.r")
  val decomposedDistributions = NamedMap[REvaluate]("decomposedDistributions") 
  val decomposedDistributionsTMP = NamedMap[BinaryFolder]("decomposedDistributionsTMP") 


  val concatenateRecombinationEventsHelsinki = REvaluate(
      var1=recombinationSitesHelsinki,
      script="""
  
	   IN = read.table(var1,sep="\t",stringsAsFactors=F,skip=3,header=F)
	   header = read.table(var1,sep="\t",stringsAsFactors=F,nrows=3,header=F)
	   table.out=data.frame()
	   raw.offspring = data.frame()
	   
	   lgs=paste0("LG",1:21) 
       for(LG in lgs){
            d=subset(IN,V1==LG)
			for(i in 6:ncol(d)){
				MALE = unlist(strsplit(header[1,i],split="f"))[2]
				FEMALE = paste0(unlist(strsplit(header[1,i],"-"))[1],"-",unlist(strsplit(header[1,i],"-"))[2])
				FAMILY= header[1,i]
				OFFSPRING = header[2,i] 
				SEX = header[3,i]
				 
				tmp = data.frame(d$V1,d$V2,d$V4,d$V5,unlist(sapply(strsplit(d[,i]," "),function(x) as.numeric(x[1]))),unlist(sapply(strsplit(d[,i]," "),function(x) as.numeric(x[2]))))
				colnames(tmp) = c("CHR","SITE","PATERNALMAP","MATERNALMAP","PATERNAL","MATERNAL")
				
				sites.paternal = c()
				precise.sites.paternal = c()
				precise.sites.paternal.CM = c()
				sites.maternal = c()					
				precise.sites.maternal = c()					
				precise.sites.maternal.CM = c()					
				intervals.paternal = c()
				intervals.maternal = c()				
				
				first.paternal.non.zero=min(which(tmp$PATERNAL!=0))
				first.maternal.non.zero=min(which(tmp$MATERNAL!=0))
				
				previous.site.MATERNAL = tmp[first.maternal.non.zero,"SITE"]
				previous.site.MATERNAL.CM = tmp[first.maternal.non.zero,"MATERNALMAP"]
				previous.site.PATERNAL = tmp[first.paternal.non.zero,"SITE"]
				previous.site.PATERNAL.CM = tmp[first.paternal.non.zero,"PATERNALMAP"]
				previous.gt.PATERNAL = tmp[first.paternal.non.zero,"PATERNAL"]
				previous.gt.MATERNAL = tmp[first.maternal.non.zero,"MATERNAL"]
				haplotypes.paternal=c(previous.gt.PATERNAL)
				haplotypes.maternal=c(previous.gt.MATERNAL)
				
				for (ROW in 1:nrow(tmp)){
				  paternal.gt = tmp[ROW,"PATERNAL"]			    
				  maternal.gt = tmp[ROW,"MATERNAL"]			    
				  if(ROW>=first.paternal.non.zero){  
					  #PATERNAL			    
					  if (paternal.gt == previous.gt.PATERNAL){
						  previous.site.PATERNAL = tmp[ROW,"SITE"]
						  previous.site.PATERNAL.CM = tmp[ROW,"PATERNALMAP"]
						  
					  } else if (paternal.gt != previous.gt.PATERNAL & paternal.gt != 0){
					     precise.sites.paternal = c(precise.sites.paternal, floor(mean(c(previous.site.PATERNAL,tmp[ROW,"SITE"]))))
					     precise.sites.paternal.CM = c(precise.sites.paternal.CM, mean(c(previous.site.PATERNAL.CM,tmp[ROW,"PATERNALMAP"])))
			             previous.gt.PATERNAL = paternal.gt
						 haplotypes.paternal=c(haplotypes.paternal,previous.gt.PATERNAL)
						 previous.site.PATERNAL = tmp[ROW,"SITE"]
						 previous.site.PATERNAL.CM = tmp[ROW,"PATERNALMAP"]
					  }
				  }
			      if(ROW>=first.maternal.non.zero){
				      #MATERNAL	  
					  if (maternal.gt == previous.gt.MATERNAL){
						  previous.site.MATERNAL = tmp[ROW,"SITE"]
						  previous.site.MATERNAL.CM = tmp[ROW,"MATERNALMAP"]
						  	  
					  } else if (maternal.gt != previous.gt.MATERNAL & maternal.gt != 0){
					     precise.sites.maternal = c(precise.sites.maternal, floor(mean(c(previous.site.MATERNAL,tmp[ROW,"SITE"]))))
					     precise.sites.maternal.CM = c(precise.sites.maternal.CM, mean(c(previous.site.MATERNAL.CM,tmp[ROW,"MATERNALMAP"])))
			             previous.gt.MATERNAL = maternal.gt
			             haplotypes.maternal=c(haplotypes.maternal,previous.gt.MATERNAL)
						 previous.site.MATERNAL = tmp[ROW,"SITE"]
						 previous.site.MATERNAL.CM = tmp[ROW,"MATERNALMAP"]
					  }
				  }
	            } #FOR INDIVIDUAL ENDS
	  
	            raw.offspring.tmp = data.frame(OFFSPRING=OFFSPRING, MALE=MALE, FEMALE=FEMALE,CHR=LG,
	            MATERNALCOUNT=length(precise.sites.maternal),PATERNALCOUNT=length(precise.sites.paternal),MATERNALSITES=paste(precise.sites.maternal,collapse=","),
			    SEX=SEX,PATERNALSITES=paste(precise.sites.paternal,collapse=","), PATERNALSITESCM=paste(precise.sites.paternal.CM,collapse=","), MATERNALSITESCM=paste(precise.sites.maternal.CM,collapse=","),
			    PATERNALHAPLOTYPES=paste(haplotypes.paternal,collapse=","),MATERNALHAPLOTYPES=paste(haplotypes.maternal,collapse=","))
			    raw.offspring=rbind(raw.offspring,raw.offspring.tmp)
	  
	       }#FOR DATA ENDS
       } #FOR LG ENDS

  max.co.count=max(c(raw.offspring$MATERNALCOUNT,raw.offspring$PATERNALCOUNT))
  for(i in 1:max.co.count){
	  raw.offspring[[paste0("MATERNALSITE",i)]]=sapply(strsplit(as.character(raw.offspring$MATERNALSITES),split=","),function(x) x[i])
	  raw.offspring[[paste0("MATERNALSITE",i,"CM")]]=sapply(strsplit(as.character(raw.offspring$MATERNALSITESCM),split=","),function(x) x[i])
	  raw.offspring[[paste0("PATERNALSITE",i)]]=sapply(strsplit(as.character(raw.offspring$PATERNALSITES),split=","),function(x) x[i])
	  raw.offspring[[paste0("PATERNALSITE",i,"CM")]]=sapply(strsplit(as.character(raw.offspring$PATERNALSITESCM),split=","),function(x) x[i])
  }
  table.out=raw.offspring
  """
  )
  
  
  

  val runEmAlgrithm = REvaluate(
	    var1 = emAlgorithmFunctions,
	    table1 = concatenateRecombinationEventsHelsinki.table,
	    script="""
		    source(var1)
		    table.out=data.frame()
		    OUT.ARRAY=list()
		     
		    lgs = paste0("LG",seq(1,21))
		    lgs=lgs[1:21]
	  	    maximum.co.count=max(c(table1$MATERNALCOUNT,table1$PATERNALCOUNT))
		    OUT.PATERNAL=data.frame(CHR=lgs)
		    OUT.MATERNAL=data.frame(CHR=lgs)
		    for(i in 0:maximum.co.count){
		        OUT.PATERNAL[[paste0("n",i)]]=100
		        OUT.MATERNAL[[paste0("n",i)]]=100
		    }
		    for(i in 0:((2*maximum.co.count-1))){
		        OUT.PATERNAL[[paste0("MLEp",i)]]=100
		        OUT.MATERNAL[[paste0("MLEp",i)]]=100
		    }
		    BootstrapPvalue=100
		    LowerBound=100
		    UpperBound=100
		    OUT.PATERNAL=cbind(OUT.PATERNAL,BootstrapPvalue,LowerBound,UpperBound)
		    OUT.MATERNAL=cbind(OUT.MATERNAL,BootstrapPvalue,LowerBound,UpperBound)
		    
			ITER=5000
			ITER.CONF=5000
	 		r=1000
	 		ALPHA=0.025		    
	   	    PERFAMILY=FALSE
	   	    
		    main=function(DATA,ITER,r,ALPHA){
		        M=sum(DATA)
		        if(length(DATA)<6){
	 		        DATA.print=c(DATA,rep(0,6-length(DATA)))
	 		    } else{
	 		        DATA.print=DATA
				}
	 		    
	 		    
	 		    p.non.restricted = MLEforP(DATA,ITER,FALSE,FALSE,TRUE,prepare.prior(DATA,TRUE))
	 		    q.non.restricted = generate.Qp(DATA,p.non.restricted)
	 		    Generated.data.sets = generateDataFrom.Qp(M,r,q.non.restricted,DATA)
	 		    All.MLEs = MLE.From.Genrated.Data(Generated.data.sets,ITER.CONF)
	            conf.intervals = Find.Conf.Intervals(All.MLEs,ALPHA)
	 		    
	 		    		    
	 		    p.restricted = MLEforP(DATA,ITER,FALSE,TRUE,FALSE,3)
	 		    q.restricted = generate.Qp(DATA,p.restricted)    
	            generated.restricted.datasets = generateDataFrom.Qp(M,r,q.restricted,DATA)
	            
	            L.non.restricted=calculate.likelihood.Function(q.non.restricted,DATA,TRUE)
	            L.restricted=calculate.likelihood.Function(q.restricted,DATA,TRUE)
	            delta.N = 2*L.non.restricted - 2*L.restricted
	            delta.N.simulated = calculate.delta.n(delta.N,generated.restricted.datasets,q.restricted,q.non.restricted,TRUE)
	            P.val=calculate.p.value(delta.N,delta.N.simulated)
	            
	 		    
	 		    if(length(p.non.restricted)<10){
	 		        OUT.VECTOR=c(DATA.print,round(p.non.restricted,8),rep(0,10-length(p.non.restricted)),P.val,conf.intervals)
	 		    } else {
			        OUT.VECTOR=c(DATA.print,round(p.non.restricted[1:10],8),P.val,conf.intervals)
			    }
	 		    
		       return(OUT.VECTOR)
		       	    
		    }
		    
		    OUT.PATERNAL.PER.FAMILY=data.frame()
		    OUT.MATERNAL.PER.FAMILY=data.frame()
		    for(LG in lgs){
		        perOffspring=subset(table1,CHR==LG)
	 		    families = names(table(perOffspring$MALE)[table(perOffspring$MALE)>20])
	 		    
	 		    if(PERFAMILY){
	 		       OUT.PATERNAL.PER.FAMILY.TMP=data.frame(CHR=rep(LG,length(families)),n0=100,n1=100,n2=100,n3=100,n4=100,n5=100,MLEp0=100,MLEp1=100,MLEp2=100,MLEp3=100,MLEp4=100,MLEp5=100,MLEp6=100,MLEp7=100,MLEp8=100,MLEp9=100,BootstrapPvalue=100,LowerBound=100,UpperBound=100,LG.LENGTH=rep(LG.length,length(families)),FAMILY=families)
	 		       OUT.MATERNAL.PER.FAMILY.TMP=data.frame(CHR=rep(LG,length(families)),n0=100,n1=100,n2=100,n3=100,n4=100,n5=100,MLEp0=100,MLEp1=100,MLEp2=100,MLEp3=100,MLEp4=100,MLEp5=100,MLEp6=100,MLEp7=100,MLEp8=100,MLEp9=100,BootstrapPvalue=100,LowerBound=100,UpperBound=100,LG.LENGTH=rep(LG.length,length(families)),FAMILY=families)
	 		       for(p in families){
	 		           maternal.data=concatenate.data(subset(perOffspring, MALE==p)$MATERNALCOUNT)
	 		           paternal.data=concatenate.data(perOffspring[perOffspring$MALE==p,"PATERNALCOUNT"])
	 		           
	 		           OUT.PATERNAL.PER.FAMILY.TMP[OUT.PATERNAL.PER.FAMILY.TMP$FAMILY==p,2:20] = main(paternal.data,ITER,r,ALPHA)
	 		           OUT.MATERNAL.PER.FAMILY.TMP[OUT.MATERNAL.PER.FAMILY.TMP$FAMILY==p,2:20] = main(maternal.data,ITER,r,ALPHA)
	 		        
	 		       }
	 		       
	 		     OUT.PATERNAL.PER.FAMILY=rbind(OUT.PATERNAL.PER.FAMILY,OUT.PATERNAL.PER.FAMILY.TMP)  
	 		     OUT.MATERNAL.PER.FAMILY=rbind(OUT.MATERNAL.PER.FAMILY,OUT.MATERNAL.PER.FAMILY.TMP)
	 		       
			    } 
			    
			    if(TRUE){
		 		    maternal.data=concatenate.data(perOffspring$MATERNALCOUNT)
		 		    paternal.data=concatenate.data(perOffspring$PATERNALCOUNT)
		 		    
		 		    OUT.PATERNAL[OUT.PATERNAL$CHR==LG,2:20]=main(paternal.data,ITER,r,ALPHA)
		 		    OUT.MATERNAL[OUT.MATERNAL$CHR==LG,2:20]=main(maternal.data,ITER,r,ALPHA)
		 		    
		 	    }
		    }
	 		
	 		OUT.ARRAY[["PATERNAL"]]=OUT.PATERNAL
	 		OUT.ARRAY[["MATERNAL"]]=OUT.MATERNAL
	 		if(PERFAMILY){
	 		    OUT.ARRAY[["PATERNAL.PER.FAMILY"]]=OUT.PATERNAL.PER.FAMILY
	 		    OUT.ARRAY[["MATERNAL.PER.FAMILY"]]=OUT.MATERNAL.PER.FAMILY
			}
			array.out=OUT.ARRAY
			
			colnames(OUT.PATERNAL)=paste0("PATERNAL",colnames(OUT.PATERNAL))
			colnames(OUT.MATERNAL)=paste0("MATERNAL",colnames(OUT.MATERNAL))
				    
		    """
		    )


      
      val chromosomeMetadata=REvaluate(
      var1=v7Index,
      var2=recombinationSitesHelsinki,
      table1=concatenateRecombinationEventsHelsinki.table,
      table2=ninespinedCentromeres,
      inArray=runEmAlgrithm.outArray,
      script="""
		
		lgs=paste0("LG",1:21)
		v7Index=read.table(var1,sep="\t",header=F,stringsAsFactors=F)
		map.data=read.table(var2,sep="\t",header=F,stringsAsFactors=F,skip=3)
		
		LENGTH=c()
		CENTROMERESTART=c()
		CENTROMERE=c()
		CENTROMEREEND=c()
		CM.male=c()
		CM.female=c()
		centromere.maternal=c()
		centromere.paternal=c()
		centromere.start.paternal=c()
		centromere.start.maternal=c()
		centromere.end.paternal=c()
		centromere.end.maternal=c()
		PARMLENGTH=c()
		QARMLENGTH=c()
		PARMLENGTHCMMATERNAL=c()
		PARMLENGTHCMPATERNAL=c()
		QARMLENGTHCMMATERNAL=c()
		QARMLENGTHCMPATERNAL=c()
		MINCOCOUNTMATERNAL=c()
		MINCOCOUNTPATERNAL=c()
		MAXCOCOUNTMATERNAL=c()
		MAXCOCOUNTPATERNAL=c()
		MEANCOCOUNTMATERNAL=c()
		MEANCOCOUNTPATERNAL=c()
		
		for(LG in lgs){
		    lg.length=max(v7Index[v7Index$V1==LG,2])
		    map=subset(map.data, V1==LG)
		    centromere=table2[table2$CHR==LG,"CENTROMERE"]
		    centromere.start=table2[table2$CHR==LG,"CENTROMERESTART"]
		    centromere.end=table2[table2$CHR==LG,"CENTROMEREEND"]
		    map.length.paternal=max(map[,4])
		    map.length.maternal=max(map[,5])
		    co.data=subset(table1,CHR==LG)
		    
		    LENGTH=c(LENGTH,lg.length)
		    CENTROMERESTART=c(CENTROMERESTART,centromere.start)
		    CENTROMERE=c(CENTROMERE,centromere)
		    CENTROMEREEND=c(CENTROMEREEND,centromere.end)
		    CM.male=c(CM.male,map.length.paternal)
		    CM.female=c(CM.female,map.length.maternal)
		    centromere.paternal=c(centromere.paternal,map[which.min(abs(map[[2]]-centromere)),4])
		    centromere.maternal=c(centromere.maternal,map[which.min(abs(map[[2]]-centromere)),5])
		    centromere.start.paternal=c(centromere.start.paternal,map[which.min(abs(map[[2]]-centromere.start)),4])
		    centromere.start.maternal=c(centromere.start.maternal,map[which.min(abs(map[[2]]-centromere.start)),5])
		    centromere.end.paternal=c(centromere.end.paternal,map[which.min(abs(map[[2]]-centromere.end)),4])
		    centromere.end.maternal=c(centromere.end.maternal,map[which.min(abs(map[[2]]-centromere.end)),5])
		    PARMLENGTH=c(PARMLENGTH,min(centromere,lg.length-centromere))
		    QARMLENGTH=c(QARMLENGTH,max(centromere,lg.length-centromere))
		    PARMLENGTHCMMATERNAL=c(PARMLENGTHCMMATERNAL,min(map[which.min(abs(map[[2]]-centromere)),5],map.length.maternal-map[which.min(abs(map[[2]]-centromere)),5]))
		    PARMLENGTHCMPATERNAL=c(PARMLENGTHCMPATERNAL,min(map[which.min(abs(map[[2]]-centromere)),4],map.length.paternal-map[which.min(abs(map[[2]]-centromere)),4]))
		    QARMLENGTHCMMATERNAL=c(QARMLENGTHCMMATERNAL,max(map[which.min(abs(map[[2]]-centromere)),5],map.length.maternal-map[which.min(abs(map[[2]]-centromere)),5]))
		    QARMLENGTHCMPATERNAL=c(QARMLENGTHCMPATERNAL,max(map[which.min(abs(map[[2]]-centromere)),4],map.length.paternal-map[which.min(abs(map[[2]]-centromere)),4]))
		    MINCOCOUNTMATERNAL=c(MINCOCOUNTMATERNAL,min(co.data$MATERNALCOUNT))
		    MINCOCOUNTPATERNAL=c(MINCOCOUNTPATERNAL,min(co.data$PATERNALCOUNT))
		    MAXCOCOUNTMATERNAL=c(MAXCOCOUNTMATERNAL,max(co.data$MATERNALCOUNT))
		    MAXCOCOUNTPATERNAL=c(MAXCOCOUNTPATERNAL,max(co.data$PATERNALCOUNT))    
		    MEANCOCOUNTMATERNAL=c(MEANCOCOUNTMATERNAL,mean(co.data$MATERNALCOUNT))
		    MEANCOCOUNTPATERNAL=c(MEANCOCOUNTPATERNAL,mean(co.data$PATERNALCOUNT))
		}
		chromosome.meta.data=data.frame(CHR=lgs,LENGTH,CHRLENGTHCMMATERNAL=CM.female,CHRLENGTHCMPATERNAL=CM.male,CENTROMERESTART,CENTROMERE,CENTROMEREEND,PARMLENGTH,QARMLENGTH,PARMLENGTHCMMATERNAL,PARMLENGTHCMPATERNAL,QARMLENGTHCMMATERNAL,QARMLENGTHCMPATERNAL,CENTROMEREMATERNALCM=centromere.maternal,CENTROMEREPATERNALCM=centromere.paternal,CENTROMERESTARTMATERNALCM=centromere.start.maternal,CENTROMERESTARTPATERNALCM=centromere.start.paternal,CENTROMEREENDMATERNALCM=centromere.end.maternal,CENTROMEREENDPATERNALCM=centromere.end.paternal,MINCOCOUNTMATERNAL,MINCOCOUNTPATERNAL,MAXCOCOUNTMATERNAL,MAXCOCOUNTPATERNAL,MEANCOCOUNTMATERNAL,MEANCOCOUNTPATERNAL)
		paternal=array[["PATERNAL"]]
		maternal=array[["MATERNAL"]]
		colnames(paternal)=paste0("PATERNAL",colnames(paternal))
		colnames(maternal)=paste0("MATERNAL",colnames(maternal))
		table.out=cbind(chromosome.meta.data,paternal[,-1],maternal[,-1])
      """
      )
      
      val fitModels=REvaluate(
          table1=concatenateRecombinationEventsHelsinki.table,
          table2=chromosomeMetadata.table,
          var1=gameteGeneratingFunctions,
          inArray=runEmAlgrithm.outArray,
          script="""
		  library(lme4)
		  library(xoi)
		  ##################
	      table.out=data.frame()
	      rm('optOut1')
		  fig.dir <- get.output(cf, 'optOut1')
		  dir.create(fig.dir, recursive=TRUE)
		  setwd(fig.dir)
		  source(var1)
		  ##################			
		  h=9
			pascalTriangle <- function(h) {
			   lapply(0:h, function(i) choose(i, 0:i))
			}
			P=pascalTriangle(h)
			sample.co.count=function(CHR,PROB.TABLE){  
			   probs=as.numeric(PROB.TABLE[PROB.TABLE$CHR==CHR,8:17])
			   co.count=sample(x=seq(0,9),size=1,prob=probs)
			   observed.co.count=sample(x=seq(0,co.count),size=1,prob=P[[co.count+1]]*0.5^co.count)
			   return(observed.co.count)
			}
		  ##################
		  lgs = paste0("LG",seq(1,21))
		  
		  
          p.val.stahl.maternal=c()
          p.val.stahl.paternal=c()
          p.val.gamma.maternal=c()
          p.val.gamma.paternal=c()
          p.val.yf.maternal=c()
          p.val.yf.paternal=c()
         
          maternal.gamma.estimate=c()
          maternal.gamma.estimate.sd=c()
          maternal.gamma.logLik=c()
          maternal.stahl.estimate=c()
          maternal.stahl.logLik=c()
          maternal.sd=c()
          maternal.p.sd=c()
          maternal.gamma.sd=c()
         
          paternal.gamma.estimate=c()
          paternal.gamma.estimate.sd=c()
          paternal.gamma.logLik=c()
          paternal.stahl.estimate=c()
          paternal.stahl.estimate.sd=c()
          paternal.stahl.logLik=c()
          paternal.sd=c()
          paternal.gamma.sd=c()
          estimates=data.frame()
          model.fit=data.frame()
          l=list()
          assess.second.derivative=function(h.values,nuhat,phat,data,LOGLIK,CHRLEN,FORNU,FORP,MAX.CONV){	
			 if(FORNU){
				 logLik.1=stahlLoglik(xoloc=data, chrlen = CHRLEN, nu=nuhat+h.values, p=phat, max.conv = MAX.CONV)
				 logLik.2=stahlLoglik(xoloc=data, chrlen = CHRLEN, nu=nuhat-h.values, p=phat, max.conv = MAX.CONV)
		     }   
	         if(FORP){
				 logLik.1=stahlLoglik(xoloc=data, chrlen = CHRLEN, nu=nuhat, p=phat+h.values, max.conv = MAX.CONV)
				 logLik.2=stahlLoglik(xoloc=data, chrlen = CHRLEN, nu=nuhat, p=phat-h.values, max.conv = MAX.CONV)
		     }   
	         second.derivative=((logLik.1)-(2*LOGLIK)+(logLik.2))/(h.values^2)  
	         out.list=list()
	         out.list[[1]]=second.derivative
	         out.list[[2]]=logLik.1
	         out.list[[3]]=logLik.2
	         return(out.list)
		 }
         
         pooled.maternal=list()
         pooled.paternal=list()
         for(LG in lgs) {
            for(SEX in c("MATERNAL","PATERNAL")){
            CHR.LENGTH.CM=table2[table2$CHR==LG,paste0("CHRLENGTHCM",SEX)]
            MAX.CONV=max(which(as.numeric(table2[table2$CHR==LG,paste0(SEX,"MLEp",1:9)])>0))
            PROB.TABLE=array[[SEX]]
            chr.length.cm.maternal=table2[table2$CHR==LG,"CHRLENGTHCMMATERNAL"]
            chr.length.cm.paternal=table2[table2$CHR==LG,"CHRLENGTHCMPATERNAL"]
            d=subset(table1,CHR==LG)
		    #####
		    
		    ##################################
		    ##-QUANTIFY INTERFERENCE STARTS-##
		    ##################################
		    #STAHL MODEL
		    		    
		    #################
		    specialD=as.matrix(d[d[[paste0(SEX,"COUNT")]]>1,paste0(SEX,"SITE",1:5,"CM")])

		    xoloc.list = as.list(apply(d[,paste0(SEX,"SITE",1:5,"CM")],1,function(x) sort(as.numeric(x[1:5])[is.na(as.numeric(x[1:5]))==F])))
		    
		    if(SEX=="MATERNAL"){pooled.maternal=append(pooled.maternal,xoloc.list)}
		    if(SEX=="PATERNAL"){pooled.paternal=append(pooled.paternal,xoloc.list)}
		  
		    censor1=c(as.numeric(d[d[[paste0(SEX,"COUNT")]]==1,paste0(SEX,"SITE1CM")]),as.numeric(unlist(apply(specialD,1,function(x) min(as.numeric(x[1:5])[is.na(as.numeric(x[1:5]))==F])))))
		    censor0=as.numeric(unlist(apply(specialD,1,function(x) diff(sort(as.numeric(x[1:5])[is.na(as.numeric(x[1:5]))==F])))))
		    censor2=c(CHR.LENGTH.CM-as.numeric(d[d[[paste0(SEX,"COUNT")]]==1,paste0(SEX,"SITE1CM")]),as.numeric(unlist(apply(specialD,1,function(x) CHR.LENGTH.CM-max(as.numeric(x[1:5])[is.na(as.numeric(x[1:5]))==F])))))
		    censor3=rep(CHR.LENGTH.CM,nrow(d[d[[paste0(SEX,"COUNT")]]==0,]))
		    interCO.distances=c(censor1,censor0,censor2,censor3)
		  
			simulate.yf=replicate(10000,sample.co.count(LG,PROB.TABLE))
			stahl.out=fitStahl(xoloc.list, CHR.LENGTH.CM,nu=c(1,20),p=0.01,max.conv=MAX.CONV,max.subd = 5000,min.subd = 100)
			simulated.stahl=simStahl(n.sim=10000,nu=median(c(1,round(stahl.out[1]),5000)),p=max(10^(-6),stahl.out[2]),L=CHR.LENGTH.CM,obligate_chiasma=T)
			simulated.gamma=simStahl(n.sim=10000,nu=max(c(1,round(stahl.out[4]))),p=0,L=CHR.LENGTH.CM,obligate_chiasma=T)
			  
			        
		    MAX.COUNT=max(c(d[[paste0(SEX,"COUNT")]],unlist(lapply(simulated.stahl,function(x) length(x)))))
		        TEST.MATRIX=matrix(c(as.numeric(table(factor(d[[paste0(SEX,"COUNT")]],levels=seq(0,MAX.COUNT)))),as.numeric(table(factor(unlist(lapply(simulated.stahl,function(x) length(x))),levels=seq(0,MAX.COUNT))))),byrow=T,nrow=2)
		        Xsq.stahl3 <- chisq.test(x=TEST.MATRIX)
		        MAX.COUNT=max(c(d[[paste0(SEX,"COUNT")]],unlist(lapply(simulated.gamma,function(x) length(x)))))
		        TEST.MATRIX=matrix(c(as.numeric(table(factor(d[[paste0(SEX,"COUNT")]],levels=seq(0,MAX.COUNT)))),as.numeric(table(factor(unlist(lapply(simulated.gamma,function(x) length(x))),levels=seq(0,MAX.COUNT))))),nrow=2,byrow=T)
		        Xsq.gamma3 <- chisq.test(x=TEST.MATRIX)
		        MAX.COUNT=max(c(d[[paste0(SEX,"COUNT")]],unlist(lapply(simulated.stahl,function(x) length(x))),unlist(lapply(simulated.gamma,function(x) length(x)))))
			    MAX.COUNT=max(c(d[[paste0(SEX,"COUNT")]],simulate.yf))
		        TEST.MATRIX=matrix(c(as.numeric(table(factor(d[[paste0(SEX,"COUNT")]],levels=seq(0,MAX.COUNT)))),as.numeric(table(factor(simulate.yf,levels=seq(0,MAX.COUNT))))),nrow=2,byrow=T)
		        Xsq.yf2 <- chisq.test(x=TEST.MATRIX)
		        
		        l[[paste0(SEX,"modelFit")]]=rbind(l[[paste0(SEX,"modelFit")]],c(as.numeric(Xsq.stahl3[1:3]),as.numeric(Xsq.gamma3[1:3]),as.numeric(Xsq.yf2[1:3])))
                
                h.values=seq(0.5,0.01,length.out=50)*stahl.out[1]
                L=assess.second.derivative(h.values,nuhat=stahl.out[1],phat=stahl.out[2],data=xoloc.list,LOGLIK=stahl.out[3],CHRLEN=CHR.LENGTH.CM,FORNU=T,FORP=F,MAX.CONV)
                second.derivative=L[[1]]
                likelihood.function=c(L[[3]],stahl.out[3],rev(L[[2]]))
                
                
                STAHLVSD=sqrt(-1/second.derivative[50])
                
                h.values=seq(0.5,0.01,length.out=50)*stahl.out[2]
                if(stahl.out[2]>0.01){
                L=assess.second.derivative(h.values,nuhat=stahl.out[1],phat=stahl.out[2],data=xoloc.list,LOGLIK=stahl.out[3],CHRLEN=CHR.LENGTH.CM,FORNU=F,FORP=T,MAX.CONV)
                second.derivative=L[[1]]
                likelihood.function=c(L[[3]],stahl.out[3],rev(L[[2]]))
                STAHLPSD=sqrt(-1/second.derivative[50])
			    } else {STAHLPSD=0} 
                h.values=seq(0.5,0.01,length.out=50)*stahl.out[4]
                L=assess.second.derivative(h.values,nuhat=stahl.out[4],phat=0,data=xoloc.list,LOGLIK=stahl.out[5],CHRLEN=CHR.LENGTH.CM,FORNU=T,FORP=F,MAX.CONV)
                second.derivative.gamma=L[[1]]
                GAMMASD=sqrt(-1/second.derivative.gamma[50])
                E=as.numeric(stahl.out)
                l[[paste0(SEX,"estimatesStahl")]]=rbind(l[[paste0(SEX,"estimatesStahl")]],as.numeric(c(E[1],STAHLVSD,E[2],STAHLPSD,E[3:4],GAMMASD,E[5:6])))
                
		  
	  }#FORSEX ENDS
	  
  }#FOR LG ENDS
		  
		  
		 estimates=data.frame(cbind(l[["MATERNALestimatesStahl"]],l[["PATERNALestimatesStahl"]]))
		 estimates=cbind(lgs,estimates)
	     colnames(estimates)=c("CHR",paste0("MATERNAL",c("vStahl","vStahlSd","pStahl","pStahlSd","logLikStahl","vGamma","vGammaSd","logLikGamma","likelihood ratio")),paste0("PATERNAL",c("vStahl","vStahlSd","pStahl","pStahlSd","logLikStahl","vGamma","vGammaSd","logLikGamma","likelihood ratio")))
	     model.fit=data.frame(cbind(lgs,l[["MATERNALmodelFit"]],l[["PATERNALmodelFit"]]))
	     colnames(model.fit)=c("CHR",paste0("MATERNAL",c("Chisq.Stahl","dfStahl","pStahl","Chisq.Gamma","dfGamma","pGamma","Chisq.YF","dfYF","pYF")),paste0("PATERNAL",c("Chisq.Stahl","dfStahl","pStahl","Chisq.Gamma","dfGamma","pGamma","Chisq.YF","dfYF","pYF")))
	     out.list=list()
	     out.list[["MODELFIT"]]=model.fit
	     out.list[["ESTIMATES"]]=estimates
	      		  
	    test.results=list()
		test.numbers=seq(0.001,5,length.out=5000)
		for(SEX in c("PATERNAL","MATERNAL")){
			for(i in 1:21){
				   n.gametes=10000
			       CHR=table2[i,"CHR"]
			       PROBS=as.numeric(table2[i,paste0(SEX,"MLEp",0:9)])
			       test.data=as.numeric(table2[i,paste0(SEX,"n",0:5)])[which(table2[i,paste0(SEX,"n",0:5)]>0)]
			       LAMBDA.POIS=2*(sum(table2[i,paste0(SEX,"n",0:5)]*0:5)/sum(table2[i,paste0(SEX,"n",0:5)]))-1
			       if(LAMBDA.POIS<0){LAMBDA.POIS=0.001}
			       EXP=2*(sum(table2[i,paste0(SEX,"n",0:5)]*0:5)/sum(table2[i,paste0(SEX,"n",0:5)]))
			       LAMBDA.trunc.pois=test.numbers[which.min(abs((test.numbers/(1-exp(-test.numbers)))-EXP))]
			       p_hat.geom=1/(1+2*(sum(test.data*(0:(length(test.data)-1)))/sum(table2[i,paste0(SEX,"n",0:5)]))-1)
			       if(p_hat.geom>1){p_hat.geom=0.999}
			       t1=test.fit.pois(1,test.data,lambda=LAMBDA.POIS,n.gametes)
			       t2=test.fit.trunc.pois(1,test.data,lambda=LAMBDA.trunc.pois,n.gametes)
				   t3=test.fit.geom(1,test.data,p_hat.geom,n.gametes)
				   t4=test.fit.yf(1,test.data,probs=PROBS,n.gametes)
			       
			       test.results[[SEX]]=rbind(test.results[[SEX]],c(t1[[3]],t1[[1]],LAMBDA.POIS,t2[[3]],t2[[1]],LAMBDA.trunc.pois,t3[[3]],t3[[1]],p_hat.geom,t4[[1]]))
				}
				
			}
		collect.results=data.frame(CHR=paste0("LG",1:21))
		collect.results=cbind(collect.results,as.data.frame(test.results[["MATERNAL"]]),as.data.frame(test.results[["PATERNAL"]]))
		colnames(collect.results)=c("CHR",paste0("MATERNAL",c("ChiPOISSON","pPOISSON","LAMBDA_hat","ChiTRUNCPOISSON","pTRUNCPOISSON","TRUNC_LAMBDA_hat","ChiGEOMETRIC","pGEOMETRIC","P_hat","pYF")),paste0("PATERNAL",c("ChiPOISSON","pPOISSON","LAMBDA_hat","ChiTRUNCPOISSON","pTRUNCPOISSON","TRUNC_LAMBDA_hat","ChiGEOMETRIC","pGEOMETRIC","P_hat","pYF")))
		out.list[["ALTERNATIVEMODELS"]]=collect.results
		array.out=out.list
		  
      """
      )
      
      val estimateHeritability=REvaluate(
          table1=concatenateRecombinationEventsHelsinki.table,
          table2=chromosomeMetadata.table,
          inArray=runEmAlgrithm.outArray,
          script="""
		      table.out=data.frame()
	          rm('optOut1')
			  fig.dir <- get.output(cf, 'optOut1')
			  dir.create(fig.dir, recursive=TRUE)
			  setwd(fig.dir)
			  ##################	
		    
				library(MCMCglmm)
				library(lme4)
				library(QGglmm)
				d=table1 
				d1=subset(d,CHR!="LG12" & !(MALE %in% c("88-m-2","72-m-2","69-m-1","66-m-2")))
				d1$OFFSPRING2=paste(d1$OFFSPRING,d1$FEMALE,d1$MALE,sep=":")
				paternalAspect=aggregate(d1$PATERNALCOUNT,by=list(d1$OFFSPRING2),sum)
				paternalAspect$animal=sapply(strsplit(paternalAspect$Group.1,split=":"), function(x) x[1])
				paternalAspect$mother=sapply(strsplit(paternalAspect$Group.1,split=":"),function(x) x[2])
				paternalAspect$father=sapply(strsplit(paternalAspect$Group.1,split=":"),function(x) x[3])
				paternalAspect=paternalAspect[,c("animal","mother","father","x")]
				colnames(paternalAspect)[4]="maleR"
				maternalAspect=aggregate(d1$MATERNALCOUNT,by=list(d1$OFFSPRING2),sum)
				Data = cbind(paternalAspect,maternalAspect[[2]])
				colnames(Data)[5]="femaleR"
				Data$totalR=Data$maleR+Data$femaleR
				Ped=data.frame(animal=c(unique(Data$mother),unique(Data$father)),mother=NA,father=NA)
				Ped=rbind(Ped,Data[,c("animal","mother","father")])
						
				prior1.1 <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1,nu=0.002))
				
				NITT=10100000
		        THIN=1000
		        BURNIN=100000
				out=data.frame()
				FAMILY=c("gaussian","poisson")
				for(f in FAMILY){
				    table.vector=c()
					m1 <- MCMCglmm(totalR ~ 1, family=f, random = ~animal, pedigree = Ped,data = Data, nitt = NITT, thin = THIN, burnin = BURNIN, prior = prior1.1,verbose = FALSE)
					m2 <- MCMCglmm(femaleR ~ 1, family=f, random = ~animal, pedigree = Ped,data = Data, nitt = NITT, thin = THIN, burnin = BURNIN, prior = prior1.1,verbose = FALSE)
					m3 <- MCMCglmm(maleR ~ 1, family=f, random = ~animal, pedigree = Ped,data = Data, nitt = NITT, thin = THIN, burnin = BURNIN, prior = prior1.1,verbose = FALSE)
					for(m in list(m1,m2,m3)){
					    V=as.mcmc(rowSums(m$VCV))
					    V.pheno=posterior.mode(V)
		                V.pheno.CI=HPDinterval(V)
		                V.additive=posterior.mode(as.mcmc(m$VCV[, "animal"]))
		                V.additive.CI=HPDinterval(as.mcmc(m$VCV[, "animal"]))
		                V.resid=posterior.mode(as.mcmc(m$VCV[, "units"]))
		                V.resid.CI=HPDinterval(as.mcmc(m$VCV[, "units"]))
		                
					    posterior.heritability <- m$VCV[, "animal"]/(m$VCV[,"animal"] + m$VCV[, "units"])
					    posterior.heritability.estimate <- posterior.mode(posterior.heritability)
					    posterior.heritability.CI=HPDinterval(posterior.heritability, 0.95)
					    table.vector=c(table.vector,posterior.heritability.estimate,posterior.heritability.CI[1],posterior.heritability.CI[2],V.additive,V.additive.CI[1],V.additive.CI[2],V.resid,V.resid.CI[1],V.resid.CI[2],V.pheno,V.pheno.CI[1],V.pheno.CI[2])
				    }
				    out=rbind(out,signif(table.vector,3))
				}
				
				out=cbind(c("gaussian_empirical_data","poisson_empirical_data"),out)
				colnames(out)=c("model",paste0("totalR",c("h2","h2CIlow","h2CIup","Va","VaCIlow","VaCIup","Vr","VrCIlow","VrCIup","Vp","VpCIlow","VpCIup")),paste0("femaleR",c("h2","h2CIlow","h2CIup","Va","VaCIlow","VaCIup","Vr","VrCIlow","VrCIup","Vp","VpCIlow","VpCIup")),paste0("maleR",c("h2","h2CIlow","h2CIup","Va","VaCIlow","VaCIup","Vr","VrCIlow","VrCIup","Vp","VpCIlow","VpCIup")))
				table.out=out	
		      	
		      """
      )
    
          val estimateHeritabilitySIMULATED=REvaluate(
          table1=concatenateRecombinationEventsHelsinki.table,
          table2=chromosomeMetadata.table,
          inArray=runEmAlgrithm.outArray,
          script="""
		      table.out=data.frame()
	          rm('optOut1')
			  fig.dir <- get.output(cf, 'optOut1')
			  dir.create(fig.dir, recursive=TRUE)
			  setwd(fig.dir)
			  ##################	
		    PROB.TABLE.MALE=array[["PATERNAL"]]
		    PROB.TABLE.FEMALE=array[["MATERNAL"]]
			h=9
			pascalTriangle <- function(h) {
			   lapply(0:h, function(i) choose(i, 0:i))
			}
			P=pascalTriangle(h)
			sample.co.count=function(CHR,PROB.TABLE){  
			   probs=as.numeric(PROB.TABLE[PROB.TABLE$CHR==CHR,8:17])
			   probs.variant=probs
			   probs.variant[3]=probs.variant[3]+sum(probs[4:10])
			   probs.variant[4:10]=0
			   co.count=sample(x=seq(0,9),size=1,prob=probs.variant)
			   observed.co.count=sample(x=seq(0,co.count),size=1,prob=P[[co.count+1]]*0.5^co.count)
			   return(observed.co.count)
			}
		
		###########
		
		library(MCMCglmm)
		library(lme4)
		
		d=table1 
		d1=subset(d,CHR!="LG12" & !(MALE %in% c("88-m-2","72-m-2","69-m-1","66-m-2")))
		d1$OFFSPRING2=paste(d1$OFFSPRING,d1$FEMALE,d1$MALE,sep=":")
		paternalAspect=aggregate(d1$PATERNALCOUNT,by=list(d1$OFFSPRING2),sum)
		paternalAspect$animal=sapply(strsplit(paternalAspect$Group.1,split=":"), function(x) x[1])
		paternalAspect$mother=sapply(strsplit(paternalAspect$Group.1,split=":"),function(x) x[2])
		paternalAspect$father=sapply(strsplit(paternalAspect$Group.1,split=":"),function(x) x[3])
		paternalAspect=paternalAspect[,c("animal","mother","father","x")]
		colnames(paternalAspect)[4]="maleR"
		maternalAspect=aggregate(d1$MATERNALCOUNT,by=list(d1$OFFSPRING2),sum)
		Data = cbind(paternalAspect,maternalAspect[[2]])
		colnames(Data)[5]="femaleR"
		Data$totalR=Data$maleR+Data$femaleR
		Ped=data.frame(animal=c(unique(Data$mother),unique(Data$father)),mother=NA,father=NA)
		Ped=rbind(Ped,Data[,c("animal","mother","father")])
		
		simulatedFemaleR=matrix(ncol=nrow(Data))
		simulatedMaleR=matrix(ncol=nrow(Data))
		CHRS=paste0("LG",c(1:11,13:21))
		for(CHR in CHRS){
		    simulatedFemaleR=rbind(simulatedFemaleR,replicate(nrow(Data),sample.co.count(CHR,PROB.TABLE.FEMALE)))
		    simulatedMaleR=rbind(simulatedMaleR,replicate(nrow(Data),sample.co.count(CHR,PROB.TABLE.MALE)))
		}
		
		Data$simulatedFemaleR=colSums(simulatedFemaleR[-1,])
		Data$simulatedMaleR=colSums(simulatedMaleR[-1,])
		Data$simulatedTotalR=Data$simulatedMaleR+Data$simulatedFemaleR
		
		prior1.1 <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1,nu=0.002))
		
		NITT=10100000
		THIN=1000
		BURNIN=100000
		out=data.frame()
		FAMILY=c("gaussian","poisson")
		for(f in FAMILY){
		    table.vector=c()
			m1 <- MCMCglmm(simulatedTotalR ~ 1, family=f, random = ~animal, pedigree = Ped,data = Data, nitt = NITT, thin = THIN, burnin = BURNIN, prior = prior1.1,verbose = FALSE)
			m2 <- MCMCglmm(simulatedFemaleR ~ 1, family=f, random = ~animal, pedigree = Ped,data = Data, nitt = NITT, thin = THIN, burnin = BURNIN, prior = prior1.1,verbose = FALSE)
			m3 <- MCMCglmm(simulatedMaleR ~ 1, family=f, random = ~animal, pedigree = Ped,data = Data, nitt = NITT, thin = THIN, burnin = BURNIN, prior = prior1.1,verbose = FALSE)
			for(m in list(m1,m2,m3)){
			    V=as.mcmc(rowSums(m$VCV))
			    V.pheno=posterior.mode(V)
                V.pheno.CI=HPDinterval(V)
                V.additive=posterior.mode(as.mcmc(m$VCV[, "animal"]))
                V.additive.CI=HPDinterval(as.mcmc(m$VCV[, "animal"]))
                V.resid=posterior.mode(as.mcmc(m$VCV[, "units"]))
                V.resid.CI=HPDinterval(as.mcmc(m$VCV[, "units"]))
			    posterior.heritability <- m$VCV[, "animal"]/(m$VCV[,"animal"] + m$VCV[, "units"])
			    posterior.heritability.estimate <- posterior.mode(posterior.heritability)
			    posterior.heritability.CI=HPDinterval(posterior.heritability, 0.95) 
			    table.vector=c(table.vector,posterior.heritability.estimate,posterior.heritability.CI[1],posterior.heritability.CI[2],V.additive,V.additive.CI[1],V.additive.CI[2],V.resid,V.resid.CI[1],V.resid.CI[2],V.pheno,V.pheno.CI[1],V.pheno.CI[2])	
		    }
		    out=rbind(out,signif(table.vector,3))
		}
		
		out=cbind(c("gaussian_simulated_data","poisson_simulated_data"),out)
		colnames(out)=c("model",paste0("totalR",c("h2","h2CIlow","h2CIup","Va","VaCIlow","VaCIup","Vr","VrCIlow","VrCIup","Vp","VpCIlow","VpCIup")),paste0("femaleR",c("h2","h2CIlow","h2CIup","Va","VaCIlow","VaCIup","Vr","VrCIlow","VrCIup","Vp","VpCIlow","VpCIup")),paste0("maleR",c("h2","h2CIlow","h2CIup","Va","VaCIlow","VaCIup","Vr","VrCIlow","VrCIup","Vp","VpCIlow","VpCIup")))
		table.out=out	
      """
      )
      

	for (LG <- lgs) {
	  decomposedDistributions(LG) = REvaluate(
	    table1 = concatenateRecombinationEventsHelsinki.table,
	    table2 = chromosomeMetadata.table,
	    var1 = decompositionFunctions,
	    param1=LG,
	    script="""
	     source(var1) 	    
	     table.out = data.frame()
	     rm('optOut1')
		 fig.dir <- get.output(cf, 'optOut1')
	     dir.create(fig.dir, recursive=TRUE)
	     setwd(fig.dir)
				
		 l=list()
         LG=param1
            ITERATIONS=200
            decomposed.distributions.1=main.EM.LIKE.CLEANING.FUNCTION(table1,table2,LG,15,13,ITERATIONS,FALSE)
 	        decomposed.distributions.2=main2(table1,table2,LG,15,13,ITERATIONS,FALSE)
            
            for(SEX in c("PATERNAL","MATERNAL")){
                l[[paste0(LG,SEX,"CLEANEM")]] = decomposed.distributions.1[[SEX]]
                l[[paste0(LG,SEX,"CLEANV")]] = decomposed.distributions.2[[SEX]]
			}
		array.out=l		
	    """
	  )
	  decomposedDistributionsTMP(LG) = Array2Folder(in=decomposedDistributions(LG).outArray,fileMode="@file@")
	}
	
	val decomposedDistributionsALLTMP = Array2Folder(in = decomposedDistributionsTMP,fileMode=".") 
	val decomposedDistributionsALL = Folder2Array(folder1 = decomposedDistributionsALLTMP,filePattern="(.*).csv")
    
    	val assignCorrectCoCount=REvaluate(
	    table1 = concatenateRecombinationEventsHelsinki.table,
	    table2 = chromosomeMetadata.table,
	    var1 = decompositionFunctions,
	    inArray = decomposedDistributionsALL,
	    script="""
	     source(var1) 	    
	     table.out = data.frame()
	     lgs = paste0("LG",seq(1,21))
		 
         for (LG in lgs) {
            centromere = table2[table2$CHR==LG,"CENTROMERE"]
		    lg.length = table2[table2$CHR==LG,"LENGTH"]
		    shape = centromere/lg.length
		    q.arm.length=table2[table2$CHR==LG,"QARMLENGTH"]
		    p.arm.length=table2[table2$CHR==LG,"PARMLENGTH"]
		    IN=subset(table1,CHR==LG)
		    ###########
		    for(SEX in c("PATERNAL","MATERNAL")){
				max.count=max(IN[,paste0(SEX,"COUNT")])
				mle.estimates = table2[,c("CHR",paste0(SEX,"n",0:5),paste0(SEX,"MLEp",0:9))]
			    colnames(mle.estimates) = gsub(SEX,"",colnames(mle.estimates))
			    interval.length=GENERATE.BINS(LG.LENGTH=lg.length,BY.LENGTH=FALSE,EVEN.LENGTH=TRUE,BY.NUMBER=TRUE,LENGTH=2000000,NUMBER=15)[2]
            	no.bins=GENERATE.BINS(LG.LENGTH=lg.length,BY.LENGTH=FALSE,EVEN.LENGTH=TRUE,BY.NUMBER=TRUE,LENGTH=2000000,NUMBER=15)[1]
                probability.table = GENERATE.PROBABILITY.TABLE(mle.estimates,LG,max.count)
                
                f.priors.1=as.matrix(array[[paste0(LG,SEX,"CLEANEM")]])
                f.priors.2=as.matrix(array[[paste0(LG,SEX,"CLEANV")]])
            	sites.all = as.numeric(unlist(strsplit(IN[,paste0(SEX,"SITES")],split=",")))
	            sites.all = sites.all[is.na(sites.all)==F]
                if(shape<=0.5){
		            sites.parm = sites.all[sites.all<=centromere]
		            sites.qarm = sites.all[sites.all>centromere]
		        } else {
		            sites.parm = sites.all[sites.all>=centromere]
		            sites.qarm = sites.all[sites.all<centromere]
		        }	
            	d=PREPARE.DATA(IN,SEX,interval.length)
            	d[,paste0("TRUECOCOUNT.EM",SEX)]=100
            	d[,paste0("PROBTRUECOCOUNT.EM",SEX)]=100
            	d[,paste0("PRIORTRUECOCOUNT.EM",SEX)]=100
            	d[,paste0("TRUECOCOUNT.V",SEX)]=100
            	d[,paste0("PROBTRUECOCOUNT.V",SEX)]=100
            	d[,paste0("PRIORTRUECOCOUNT.V",SEX)]=100
            	d[,paste0("ParmCOUNT",SEX)]=100
            	d[,paste0("QarmCOUNT",SEX)]=100
            	   
            	data.vector.matrix=CONSTRUCT.DATA.VECTOR.MATRIX(d,no.bins,max.count)
            	
                for(i in 1:nrow(data.vector.matrix)){
            		data.vector=data.vector.matrix[i,]
            		co.class=sum(data.vector)
            		class.priors=as.numeric(probability.table[probability.table$OBSERVED == co.class,2:11])
            		class.posteriors=c()
            		class.posteriors=c()
            			for(j in 0:9){
            				if(class.priors[j+1]==0){
            					class.posteriors=c(class.posteriors,0)
            				} else{
            					class.posteriors=c(class.posteriors,(class.priors[j+1]*prod(f.priors.1[j+1,]^(data.vector)))/sum(class.priors*t(apply(f.priors.1,1,function(x) prod(x^(data.vector))))))
            			    }
            		    }
            		ASSIGN=sample(seq(0,9),prob=class.posteriors,1)   
            		d[i,paste0("TRUECOCOUNT.EM",SEX)]=ASSIGN
            		d[i,paste0("PROBTRUECOCOUNT.EM",SEX)]=class.posteriors[ASSIGN+1]
            		d[i,paste0("PRIORTRUECOCOUNT.EM",SEX)]=class.priors[co.class+1]
            		d[i,paste0("ParmCOUNT",SEX)]=length(which(unlist(strsplit(d[i,paste0(SEX,"SITES")],split=",")) %in% sites.parm))
		            d[i,paste0("QarmCOUNT",SEX)]=length(which(unlist(strsplit(d[i,paste0(SEX,"SITES")],split=",")) %in% sites.qarm))    	
            	
            	}
            	
                for(i in 1:nrow(data.vector.matrix)){
            		data.vector=data.vector.matrix[i,]
            		co.class=sum(data.vector)
            		class.priors=as.numeric(probability.table[probability.table$OBSERVED == co.class,2:11])
            		class.posteriors=c()
            		class.posteriors=c()
            			for(j in 0:9){
            				if(class.priors[j+1]==0){
            					class.posteriors=c(class.posteriors,0)
            				} else{
            					class.posteriors=c(class.posteriors,(class.priors[j+1]*prod(f.priors.2[j+1,]^(data.vector)))/sum(class.priors*t(apply(f.priors.2,1,function(x) prod(x^(data.vector))))))
            			    }
            		    }
            		ASSIGN=sample(seq(0,9),prob=class.posteriors,1)
            		d[i,paste0("TRUECOCOUNT.V",SEX)]=ASSIGN
            		d[i,paste0("PROBTRUECOCOUNT.V",SEX)]=class.posteriors[ASSIGN+1]
            		d[i,paste0("PRIORTRUECOCOUNT.V",SEX)]=class.priors[co.class+1]
            	}
            IN=cbind(IN,d[,paste0(c("TRUECOCOUNT.EM","PROBTRUECOCOUNT.EM","PRIORTRUECOCOUNT.EM","TRUECOCOUNT.V","PROBTRUECOCOUNT.V","PRIORTRUECOCOUNT.V","ParmCOUNT","QarmCOUNT"),SEX)])
            
            }
         
            
         table.out=rbind(table.out,IN)   
		} 	  	
	    """
	)
 
 val plotCleanedDistributions = REvaluate(
    	table1 = assignCorrectCoCount.table,
	    table2 = chromosomeMetadata.table,
	    inArray = decomposedDistributionsALL,
	    script="""
	   #######################################
	    table.out = data.frame()
	     rm('optOut1')
		 fig.dir <- get.output(cf, 'optOut1')
	     dir.create(fig.dir, recursive=TRUE)
	     setwd(fig.dir)
		FORM.BIN.LENGTH = function(LG.LENGTH,BY.LENGTH,EVEN.LENGTH,BY.NUMBER,LENGTH,NUMBER){
		
		    if(BY.NUMBER & EVEN.LENGTH){
		        INTERVAL.LENGTH = (1+0.001/NUMBER)*(LG.LENGTH/NUMBER)       
		    } else if(BY.LENGTH & EVEN.LENGTH){
		        if(LG.LENGTH%%LENGTH >= 0.5*LENGTH){
		            BIN.COUNT=length(seq(0,LG.LENGTH,LENGTH))
		            INTERVAL.LENGTH=(1+0.001/BIN.COUNT)*(seq(0,LG.LENGTH,length.out=BIN.COUNT+1)[2]-seq(0,LG.LENGTH,length.out=BIN.COUNT+1)[1])
		        } else {
		            BIN.COUNT=length(seq(0,LG.LENGTH,LENGTH))-1
		            INTERVAL.LENGTH=(1+0.001/BIN.COUNT)*(seq(0,LG.LENGTH,length.out=BIN.COUNT+1)[2]-seq(0,LG.LENGTH,length.out=BIN.COUNT+1)[1])
		        }
		    } else if(BY.LENGTH){
		        INTERVAL.LENGTH=LENGTH
		    }
		    return(INTERVAL.LENGTH)
		}
		####################################
		distribution.plotting.table=data.frame()
		OUT.LIST=list()
		 lgs = paste0("LG",seq(1,21))
		 NO.BINS=15
		 
		 for (LG in lgs) {
		    LG.length = table2[table2$CHR == LG,"LENGTH"]
		    BIN.LENGTH=FORM.BIN.LENGTH(LG.LENGTH=LG.length,BY.LENGTH=FALSE,EVEN.LENGTH=TRUE,BY.NUMBER=TRUE,LENGTH=0,NUMBER=NO.BINS)
		    CENTROMERE = table2[table2$CHR==LG,"CENTROMERE"]
		    CENTROMERE.BIN = floor(CENTROMERE/BIN.LENGTH)
		    BINS=1:(floor(LG.length/BIN.LENGTH)+1)
		    		    		    
		    d=subset(table1,CHR==LG)    
	        
	        p.sites.all=as.numeric(unlist(strsplit(d$PATERNALSITES,split=",")))
	        p.sites.all=p.sites.all[is.na(p.sites.all)==F]
	        p.sites.all.in.bins=floor(p.sites.all/BIN.LENGTH)+1
	        
	        m.sites.all=as.numeric(unlist(strsplit(d$MATERNALSITES,split=",")))
	        m.sites.all=m.sites.all[is.na(m.sites.all)==F]
	        m.sites.all.in.bins=floor(m.sites.all/BIN.LENGTH)+1
	        
	        p.sites.1=as.numeric(unlist(strsplit(subset(d,PATERNALCOUNT==1)$PATERNALSITES,split=",")))
	        p.sites.1=p.sites.1[is.na(p.sites.1)==F]
	        m.sites.1=as.numeric(unlist(strsplit(subset(d,MATERNALCOUNT==1)$MATERNALSITES,split=",")))
	        m.sites.1=m.sites.1[is.na(m.sites.1)==F]
	        
	        p.sites.2.1=as.numeric(sapply(strsplit(subset(d,PATERNALCOUNT==2)$PATERNALSITES,split=","), function(x) x[1]))
	        p.sites.2.1=p.sites.2.1[is.na(p.sites.2.1)==F]
	        p.sites.2.2=as.numeric(sapply(strsplit(subset(d,PATERNALCOUNT==2)$PATERNALSITES,split=","), function(x) x[2]))
	        p.sites.2.2=p.sites.2.2[is.na(p.sites.2.2)==F]
	        m.sites.2.1=as.numeric(sapply(strsplit(subset(d,MATERNALCOUNT==2)$MATERNALSITES,split=","), function(x) x[1]))
	        m.sites.2.1=m.sites.2.1[is.na(m.sites.2.1)==F]
	        m.sites.2.2=as.numeric(sapply(strsplit(subset(d,MATERNALCOUNT==2)$MATERNALSITES,split=","), function(x) x[2]))
	        m.sites.2.2=m.sites.2.2[is.na(m.sites.2.2)==F]
	        
	        m.sites.3.1=as.numeric(sapply(strsplit(subset(d,MATERNALCOUNT==3)$MATERNALSITES,split=","), function(x) x[1]))
	        m.sites.3.1=m.sites.3.1[is.na(m.sites.3.1)==F]
	        m.sites.3.2=as.numeric(sapply(strsplit(subset(d,MATERNALCOUNT==3)$MATERNALSITES,split=","), function(x) x[2]))
	        m.sites.3.2=m.sites.3.2[is.na(m.sites.3.2)==F]
	        m.sites.3.3=as.numeric(sapply(strsplit(subset(d,MATERNALCOUNT==3)$MATERNALSITES,split=","), function(x) x[3]))
	        m.sites.3.3=m.sites.3.3[is.na(m.sites.3.3)==F]
	         
	       for(i in 1:5){
	        d[[paste0("MATERNALSITE",i)]] = floor(d[[paste0("MATERNALSITE",i)]]/BIN.LENGTH)+1
	        d[[paste0("PATERNALSITE",i)]] = floor(d[[paste0("PATERNALSITE",i)]]/BIN.LENGTH)+1
	       }  
	       for(SEX in c("PATERNAL","MATERNAL")){
	           p=table2[,c("CHR",paste0(SEX,"n",0:5),paste0(SEX,"MLEp",0:9))]
	           colnames(p) = gsub(SEX,"",colnames(p))
	           raw.tmp=as.data.frame(d[,c(paste0(SEX,c("COUNT","SITE1","SITE2","SITE3","SITE4","SITE5")),paste0(c("TRUECOCOUNT.V","TRUECOCOUNT.EM"),SEX))])
	           colnames(raw.tmp)=c("COUNT","SITE1","SITE2","SITE3","SITE4","SITE5","TRUECOCOUNT.V","TRUECOCOUNT.EM")
	           COLORS=c("blue","red","green3","brown","orange")
		       
	           clean = array[[paste0(LG,SEX,"CLEANEM")]]
		       
	           clean = array[[paste0(LG,SEX,"CLEANV")]]
		       
	           
		   }
		   p=table2[,c("CHR",paste0("MATERNAL","n",0:5),paste0("MATERNAL","MLEp",0:9))]
	       colnames(p) = gsub("MATERNAL","",colnames(p))
	       clean.paternal = array[[paste0(LG,"PATERNALCLEANEM")]]
		   paternal.total = apply(clean.paternal*as.numeric(p[p$CHR==LG,paste0("MLEp",seq(0,9))]),2,sum)/sum(apply(clean.paternal*as.numeric(p[p$CHR==LG,paste0("MLEp",seq(0,9))]),2,sum))
		   paternal.total.E=sum(paternal.total*seq(1,length(paternal.total)))
		   raw.tmp=as.data.frame(d[,c(paste0("PATERNAL",c("COUNT","SITE1","SITE2","SITE3","SITE4","SITE5")),paste0(c("TRUECOCOUNT.V","TRUECOCOUNT.EM"),"PATERNAL"))])
	       colnames(raw.tmp)=c("COUNT","SITE1","SITE2","SITE3","SITE4","SITE5","TRUECOCOUNT.V","TRUECOCOUNT.EM")
	       raw.paternal=as.data.frame(table(factor(c(raw.tmp$SITE1,raw.tmp$SITE2,raw.tmp$SITE3,raw.tmp$SITE4,raw.tmp$SITE5),levels=BINS))/sum(table(factor(c(raw.tmp$SITE1,raw.tmp$SITE2,raw.tmp$SITE3,raw.tmp$SITE4,raw.tmp$SITE5),levels=BINS))))
		   
		   p=table2[,c("CHR",paste0("PATERNAL","n",0:5),paste0("PATERNAL","MLEp",0:9))]
	       colnames(p) = gsub("PATERNAL","",colnames(p))
	       clean.maternal = array[[paste0(LG,"MATERNALCLEANEM")]]
		   maternal.total = apply(clean.maternal*as.numeric(p[p$CHR==LG,paste0("MLEp",seq(0,9))]),2,sum)/sum(apply(clean.maternal*as.numeric(p[p$CHR==LG,paste0("MLEp",seq(0,9))]),2,sum))
		   maternal.total.E=sum(maternal.total*seq(1,length(maternal.total)))
		   raw.tmp=as.data.frame(d[,c(paste0("MATERNAL",c("COUNT","SITE1","SITE2","SITE3","SITE4","SITE5")),paste0(c("TRUECOCOUNT.V","TRUECOCOUNT.EM"),"MATERNAL"))])
	       colnames(raw.tmp)=c("COUNT","SITE1","SITE2","SITE3","SITE4","SITE5","TRUECOCOUNT.V","TRUECOCOUNT.EM")
	       raw.maternal=as.data.frame(table(factor(c(raw.tmp$SITE1,raw.tmp$SITE2,raw.tmp$SITE3,raw.tmp$SITE4,raw.tmp$SITE5),levels=BINS))/sum(table(factor(c(raw.tmp$SITE1,raw.tmp$SITE2,raw.tmp$SITE3,raw.tmp$SITE4,raw.tmp$SITE5),levels=BINS))))
		   
   		   distribution.plotting.table=rbind(distribution.plotting.table,data.frame(CHR=LG,BIN=1:NO.BINS,DENSITY=c(paternal.total,raw.paternal[,2],maternal.total,raw.maternal[,2]),SEX=rep(c("PATERNAL","MATERNAL"),each=2*NO.BINS),
   		   TYPE=rep(c("PREDICTED","EMPIRICAL","PREDICTED","EMPIRICAL"),each=NO.BINS),COMPARISON=rep(c("COMP1","COMP2","COMP2","COMP1"),each=NO.BINS),
   		   EXPECTED=rep(c(sum(paternal.total*seq(1,NO.BINS)),sum(raw.paternal[,2]*seq(1,NO.BINS)),sum(maternal.total*seq(1,NO.BINS)),sum(raw.maternal[,2]*seq(1,NO.BINS))),each=NO.BINS),
   		   VARIANCE=rep(c(sum(((seq(1,NO.BINS)-sum(paternal.total*seq(1,NO.BINS)))^2)*paternal.total),sum(((seq(1,NO.BINS)-sum(raw.paternal[,2]*seq(1,NO.BINS)))^2)*raw.paternal[,2]),
		   sum(((seq(1,NO.BINS)-sum(maternal.total*seq(1,NO.BINS)))^2)*maternal.total),sum(((seq(1,NO.BINS)-sum(raw.maternal[,2]*seq(1,NO.BINS)))^2)*raw.maternal[,2])),each=NO.BINS)))


   		   
		   #######################################
		   p=table2[,c("CHR",paste0("MATERNAL","n",0:5),paste0("MATERNAL","MLEp",0:9))]
	       colnames(p) = gsub("MATERNAL","",colnames(p))
	       clean = array[[paste0(LG,"PATERNALCLEANV")]]
		   paternal.total = apply(clean*as.numeric(p[p$CHR==LG,paste0("MLEp",seq(0,9))]),2,sum)/sum(apply(clean*as.numeric(p[p$CHR==LG,paste0("MLEp",seq(0,9))]),2,sum))
		   paternal.total.E=sum(paternal.total*seq(1,length(paternal.total)))
		   raw.tmp=as.data.frame(d[,c(paste0("PATERNAL",c("COUNT","SITE1","SITE2","SITE3","SITE4","SITE5")),paste0(c("TRUECOCOUNT.V","TRUECOCOUNT.EM"),"PATERNAL"))])
	       colnames(raw.tmp)=c("COUNT","SITE1","SITE2","SITE3","SITE4","SITE5","TRUECOCOUNT.V","TRUECOCOUNT.EM")
	       raw.paternal=as.data.frame(table(factor(c(raw.tmp$SITE1,raw.tmp$SITE2,raw.tmp$SITE3,raw.tmp$SITE4,raw.tmp$SITE5),levels=BINS))/sum(table(factor(c(raw.tmp$SITE1,raw.tmp$SITE2,raw.tmp$SITE3,raw.tmp$SITE4,raw.tmp$SITE5),levels=BINS))))
		   
		   p=table2[,c("CHR",paste0("PATERNAL","n",0:5),paste0("PATERNAL","MLEp",0:9))]
	       colnames(p) = gsub("PATERNAL","",colnames(p))
	       clean = array[[paste0(LG,"MATERNALCLEANV")]]
		   maternal.total = apply(clean*as.numeric(p[p$CHR==LG,paste0("MLEp",seq(0,9))]),2,sum)/sum(apply(clean*as.numeric(p[p$CHR==LG,paste0("MLEp",seq(0,9))]),2,sum))
		   maternal.total.E=sum(maternal.total*seq(1,length(maternal.total)))
		   raw.tmp=as.data.frame(d[,c(paste0("MATERNAL",c("COUNT","SITE1","SITE2","SITE3","SITE4","SITE5")),paste0(c("TRUECOCOUNT.V","TRUECOCOUNT.EM"),"MATERNAL"))])
	       colnames(raw.tmp)=c("COUNT","SITE1","SITE2","SITE3","SITE4","SITE5","TRUECOCOUNT.V","TRUECOCOUNT.EM")
	       raw.maternal=as.data.frame(table(factor(c(raw.tmp$SITE1,raw.tmp$SITE2,raw.tmp$SITE3,raw.tmp$SITE4,raw.tmp$SITE5),levels=BINS))/sum(table(factor(c(raw.tmp$SITE1,raw.tmp$SITE2,raw.tmp$SITE3,raw.tmp$SITE4,raw.tmp$SITE5),levels=BINS))))
		   
		   simulated.clean.maternal=sample(seq(1,NO.BINS),prob=maternal.total,size=500,replace=T)
		   simulated.clean.paternal=sample(seq(1,NO.BINS),prob=paternal.total,size=500,replace=T)
		   
		  
	   }
	   OUT.LIST[["empricalAndInferredDistributions"]]=distribution.plotting.table
      array.out=OUT.LIST
    
    """
    
    )

      val plots=REvaluate(
      table1=concatenateRecombinationEventsHelsinki.table,
      table2=chromosomeMetadata.table,
   	  table3=assignCorrectCoCount.table,
   	  inArray=fitModels.outArray,
      script="""
          ##################
	      table.out=data.frame()
	      rm('optOut1')
		  fig.dir <- get.output(cf, 'optOut1')
		  dir.create(fig.dir, recursive=TRUE)
		  setwd(fig.dir)
		  ##################			
		  lgs = paste0("LG",seq(1,21))
		  fig.index=0
          library(lme4)
          library(glmmTMB)
          OUT.LIST=list() 
          
		  ##################
          COLLECT.SHAPE.PLOTS = data.frame()
          COLLECT.DISTRIBUTION.PLOTS=data.frame()
          for(LG in lgs){
              tmp=table1[table1$CHR==LG,]
              BIN.COUNT=25
              BIN.SIZE=ceiling(table2[table2$CHR==LG,"LENGTH"]/BIN.COUNT)
              YLIM=c(0,0.3)
              paternal.co=as.numeric(prop.table(table(factor(ceiling(as.numeric(unlist(strsplit(tmp[tmp$PATERNALCOUNT>0,"PATERNALSITES"],split=",")))/BIN.SIZE),levels=1:BIN.COUNT))))
              maternal.co=as.numeric(prop.table(table(factor(ceiling(as.numeric(unlist(strsplit(tmp[tmp$MATERNALCOUNT>0,"MATERNALSITES"],split=",")))/BIN.SIZE),levels=1:BIN.COUNT))))
              COLLECT.DISTRIBUTION.PLOTS=rbind(COLLECT.DISTRIBUTION.PLOTS,data.frame(BIN=1:BIN.COUNT,DENSITY=c(paternal.co,maternal.co),SEX=rep(c("PATERNAL","MATERNAL"),each=BIN.COUNT),CO=1,CO.GROUP=1,CHR=LG))
              paternal.co.1=as.numeric(prop.table(table(factor(ceiling(as.numeric(tmp[tmp$PATERNALCOUNT==2,"PATERNALSITE1"])/BIN.SIZE),levels=1:BIN.COUNT))))
              paternal.co.2=as.numeric(prop.table(table(factor(ceiling(as.numeric(tmp[tmp$PATERNALCOUNT==2,"PATERNALSITE2"])/BIN.SIZE),levels=1:BIN.COUNT))))
              maternal.co.1=as.numeric(prop.table(table(factor(ceiling(as.numeric(tmp[tmp$MATERNALCOUNT==2,"MATERNALSITE1"])/BIN.SIZE),levels=1:BIN.COUNT))))
              maternal.co.2=as.numeric(prop.table(table(factor(ceiling(as.numeric(tmp[tmp$MATERNALCOUNT==2,"MATERNALSITE2"])/BIN.SIZE),levels=1:BIN.COUNT))))
              COLLECT.DISTRIBUTION.PLOTS=rbind(COLLECT.DISTRIBUTION.PLOTS,data.frame(BIN=1:BIN.COUNT,DENSITY=c(paternal.co.1,paternal.co.2,maternal.co.1,maternal.co.2),SEX=rep(c("PATERNAL","MATERNAL"),each=2*BIN.COUNT),CO=rep(c(1,2,1,2),each=BIN.COUNT),CO.GROUP=2,CHR=LG))
          
              maternal.co.1=as.numeric(prop.table(table(factor(ceiling(as.numeric(tmp[tmp$MATERNALCOUNT==3,"MATERNALSITE1"])/BIN.SIZE),levels=1:BIN.COUNT))))
              maternal.co.2=as.numeric(prop.table(table(factor(ceiling(as.numeric(tmp[tmp$MATERNALCOUNT==3,"MATERNALSITE2"])/BIN.SIZE),levels=1:BIN.COUNT))))
              maternal.co.3=as.numeric(prop.table(table(factor(ceiling(as.numeric(tmp[tmp$MATERNALCOUNT==3,"MATERNALSITE3"])/BIN.SIZE),levels=1:BIN.COUNT))))
              COLLECT.DISTRIBUTION.PLOTS=rbind(COLLECT.DISTRIBUTION.PLOTS,data.frame(BIN=1:BIN.COUNT,DENSITY=c(maternal.co.1,maternal.co.2,maternal.co.3),SEX="MATERNAL",CO=rep(c(1,2,3),each=BIN.COUNT),CO.GROUP=3,CHR=LG))
              
	      }
	      BIN.COUNT=25
	      YLIM=c(0,0.25)    
	      table2$TYPE="TELOCENTRIC"
	      for(LG in lgs){
	          q.arm.length=table2[table2$CHR==LG,"QARMLENGTH"]
	          p.arm.length=table2[table2$CHR==LG,"PARMLENGTH"]
	          ratio=q.arm.length/p.arm.length
	          if(ratio<1.7){
	              table2[table2$CHR==LG,"TYPE"]="METACENTRIC"
	          } else if(ratio<3){
	              table2[table2$CHR==LG,"TYPE"]="SUBMETACENTRIC"
	          } else if(ratio<7){
	              table2[table2$CHR==LG,"TYPE"]="SUBTELOCENTRIC"
	          } else if(ratio<80){
	              table2[table2$CHR==LG,"TYPE"]="ACROCENTRIC"
	          }
	               
	      }
	      for(i in paste0(c("META","SUBMETA","SUBTELO","ACRO","TELO"),"CENTRIC")){
	      for(SEX in c("MATERNAL","PATERNAL")){
	          if(SEX=="MATERNAL"){COL=rgb(1,0,0,0.5)}else{COL=rgb(0,0,1,0.5)}
	          all=matrix(ncol=BIN.COUNT)
		      for(LG in table2[table2$TYPE==i,"CHR"]){
			      tmp=table1[table1$CHR==LG,]
		          lg.length=table2[table2$CHR==LG,"LENGTH"]
	              centromere=table2[table2$CHR==LG,"CENTROMERE"]
		          shape=centromere/lg.length
		          BIN.SIZE=ceiling(lg.length/BIN.COUNT)
		          co=as.numeric(prop.table(table(factor(ceiling(as.numeric(unlist(strsplit(tmp[tmp[[paste0(SEX,"COUNT")]]>0,paste0(SEX,"SITES")],split=",")))/BIN.SIZE),levels=1:BIN.COUNT))))
		          if(shape > 0.5){co=rev(co)}
	              all=rbind(all,co)
	              COLLECT.SHAPE.PLOTS=rbind(COLLECT.SHAPE.PLOTS,data.frame(co,SEX,CHR=LG,SHAPE=i,QARMCOUNT="all",BIN=1:BIN.COUNT,COLOR=SEX))
		      }
		      COLLECT.SHAPE.PLOTS=rbind(COLLECT.SHAPE.PLOTS,data.frame(co=apply(all[-1,],2,median),SEX,CHR="MEDIAN",SHAPE=i,QARMCOUNT="all",BIN=1:BIN.COUNT,COLOR="MEDIAN"))
		  }
		  }
		  
		  for(i in paste0(c("SUBMETA","SUBTELO","ACRO","TELO"),"CENTRIC")){
		      for(j in 2:3){
		          all=matrix(ncol=BIN.COUNT)
		          for(LG in table2[table2$TYPE==i,"CHR"]){
				      tmp=table3[table3$CHR==LG,]
			          lg.length=table2[table2$CHR==LG,"LENGTH"]
		              centromere=table2[table2$CHR==LG,"CENTROMERE"]
			          shape=centromere/lg.length
			          BIN.SIZE=ceiling(lg.length/BIN.COUNT)
			          maternal.co=as.numeric(prop.table(table(factor(ceiling(as.numeric(unlist(strsplit(tmp[tmp$QarmCOUNTMATERNAL==j,"MATERNALSITES"],split=",")))/BIN.SIZE),levels=1:BIN.COUNT))))
		              if(nrow(tmp[tmp$QarmCOUNTMATERNAL==j,]) > 20){ ##Previous version: if(nrow(tmp[tmp$MATERNALCOUNT==j,])/j > 10)
			              if(shape > 0.5){paternal.co=rev(paternal.co);maternal.co=rev(maternal.co)}
			              all=rbind(all,maternal.co)
			              COLLECT.SHAPE.PLOTS=rbind(COLLECT.SHAPE.PLOTS,data.frame(co=maternal.co,SEX="MATERNAL",CHR=LG,SHAPE=i,QARMCOUNT=as.character(j),BIN=1:BIN.COUNT,COLOR="MATERNAL"))
					  }    
			      }
			      if(nrow(all)>2){
			          COLLECT.SHAPE.PLOTS=rbind(COLLECT.SHAPE.PLOTS,data.frame(co=apply(all[-1,],2,median),SEX="MATERNAL",CHR="MEDIAN",SHAPE=i,QARMCOUNT=as.character(j),BIN=1:BIN.COUNT,COLOR="MEDIAN"))
			      }
			  }
		  }
		  BIN.COUNT=25
	      for(LG in table2[table2$TYPE=="ACROCENTRIC","CHR"]){
		      tmp=table3[table3$CHR==LG,]
	          lg.length=table2[table2$CHR==LG,"QARMLENGTH"]
              centromere=table2[table2$CHR==LG,"CENTROMERE"]
	          shape=centromere/table2[table2$CHR==LG,"LENGTH"]
	          BIN.SIZE=ceiling(lg.length/BIN.COUNT)
	          if(shape<=0.5){
		          maternal.co=as.numeric(unlist(strsplit(tmp[tmp$QarmCOUNTMATERNAL==2,"MATERNALSITES"],split=",")))
		          maternal.q=maternal.co[which(maternal.co>centromere)]
		          maternal.freq=as.numeric(prop.table(table(factor(ceiling((maternal.q-centromere)/BIN.SIZE),levels=1:BIN.COUNT))))
		        } else{
	              maternal.co=as.numeric(unlist(strsplit(tmp[tmp$QarmCOUNTMATERNAL==2,"MATERNALSITES"],split=",")))
		          maternal.q=maternal.co[which(maternal.co<centromere)]
		          maternal.freq=rev(as.numeric(prop.table(table(factor(ceiling((centromere-maternal.q)/BIN.SIZE),levels=1:BIN.COUNT)))))
		        }
		      
		  }    
          ##############################
          ##DISTANCE TO CENTEROMERE
           lme.model.table=data.frame()
           for (LG in lgs[-12]) {
                d=subset(table3,CHR==LG)
                d$QARM1stPaternalPERQARMLENGTH=100
                d$QARM1stMaternalPERQARMLENGTH=100
                centromere=table2[table2$CHR==LG,"CENTROMERE"]
                q.arm.length=table2[table2$CHR==LG,"QARMLENGTH"]
                shape=table2[table2$CHR==LG,"CENTROMERE"]/table2[table2$CHR==LG,"LENGTH"]
                if(shape>0.5){
	                d$MATERNALSITE1 = sapply(strsplit(d$MATERNALSITES,split=","),function(x) as.numeric(x)[order(as.numeric(x),decreasing=T)][1])
			        d$MATERNALSITE2 = sapply(strsplit(d$MATERNALSITES,split=","),function(x) as.numeric(x)[order(as.numeric(x),decreasing=T)][2])
		            d$MATERNALSITE3 = sapply(strsplit(d$MATERNALSITES,split=","),function(x) as.numeric(x)[order(as.numeric(x),decreasing=T)][3])
		            d$MATERNALSITE4 = sapply(strsplit(d$MATERNALSITES,split=","),function(x) as.numeric(x)[order(as.numeric(x),decreasing=T)][4])
		            d$MATERNALSITE5 = sapply(strsplit(d$MATERNALSITES,split=","),function(x) as.numeric(x)[order(as.numeric(x),decreasing=T)][5])
			        d$PATERNALSITE1 = sapply(strsplit(d$PATERNALSITES,split=","),function(x) as.numeric(x)[order(as.numeric(x),decreasing=T)][1])
			        d$PATERNALSITE2 = sapply(strsplit(d$PATERNALSITES,split=","),function(x) as.numeric(x)[order(as.numeric(x),decreasing=T)][2])
		            d$PATERNALSITE3 = sapply(strsplit(d$PATERNALSITES,split=","),function(x) as.numeric(x)[order(as.numeric(x),decreasing=T)][3])
		            d$PATERNALSITE4 = sapply(strsplit(d$PATERNALSITES,split=","),function(x) as.numeric(x)[order(as.numeric(x),decreasing=T)][4])
		            d$PATERNALSITE5 = sapply(strsplit(d$PATERNALSITES,split=","),function(x) as.numeric(x)[order(as.numeric(x),decreasing=T)][5])
                }
				tmp.maternal=subset(d,CHR==LG & PRIORTRUECOCOUNT.EMMATERNAL>0 & QarmCOUNTMATERNAL>0)
				tmp.paternal=subset(d,CHR==LG & PRIORTRUECOCOUNT.EMPATERNAL>0 & QarmCOUNTPATERNAL>0)
				for(i in (0:4)){
				    tmp.maternal[tmp.maternal$ParmCOUNTMATERNAL==i,"QARM1stMaternalPERQARMLENGTH"]=abs((tmp.maternal[tmp.maternal$ParmCOUNTMATERNAL==i,paste0("MATERNALSITE",i+1)])-centromere)/q.arm.length
				    tmp.paternal[tmp.paternal$ParmCOUNTPATERNAL==i,"QARM1stPaternalPERQARMLENGTH"]=abs((tmp.paternal[tmp.paternal$ParmCOUNTPATERNAL==i,paste0("PATERNALSITE",i+1)])-centromere)/q.arm.length
				}
				
				tmp.maternal=subset(tmp.maternal,select=c("CHR","OFFSPRING","QARM1stMaternalPERQARMLENGTH","QarmCOUNTMATERNAL"))
				tmp.paternal=subset(tmp.paternal,select=c("CHR","OFFSPRING","QARM1stPaternalPERQARMLENGTH","QarmCOUNTPATERNAL"))
				colnames(tmp.maternal)=c("CHR","OFFSPRING","DISTANCETOCENTROMERE","QARMCOCOUNT")
				colnames(tmp.paternal)=c("CHR","OFFSPRING","DISTANCETOCENTROMERE","QARMCOCOUNT")
				if(nrow(tmp.maternal)>0){tmp.maternal$SEX="MATERNAL"}
				if(nrow(tmp.paternal)>0){tmp.paternal$SEX="PATERNAL"}
				lme.model.table=rbind(lme.model.table,tmp.maternal)
				lme.model.table=rbind(lme.model.table,tmp.paternal)
		}
			
		
		lme.model.table$fQARMCOCOUNT = as.factor(lme.model.table$QARMCOCOUNT)
		lme.model.table$fSEX = as.factor(lme.model.table$SEX)	
		lme.model.table$fCHR = as.factor(lme.model.table$CHR)
                lme.model.table$fOFFSPRING = as.factor(lme.model.table$OFFSPRING)
                lme.model.table.2=subset(lme.model.table,!(CHR %in% paste0("LG",c(1,12,13,16,2,4,8,9)))) ### WHY THIS EXISTS??
		OUT.LIST[["lme.model.table"]]=lme.model.table
        #########
		median.distance.to.centromere.male=c()
		median.distance.to.centromere.female=c()
		
        OUT.LIST[["shapePlots"]] = COLLECT.SHAPE.PLOTS
        OUT.LIST[["distributionPlots"]] = COLLECT.DISTRIBUTION.PLOTS
        array.out=OUT.LIST
      """
      )

val blindToCentromere = REvaluate(
	    table1 = assignCorrectCoCount.table,
	    table2 = chromosomeMetadata.table,
	    script="""
		 library(lme4)
		 library(xoi)
		 table.out = data.frame()
	     rm('optOut1')
		 fig.dir <- get.output(cf, 'optOut1')
	     dir.create(fig.dir, recursive=TRUE)
	     setwd(fig.dir)
				
		 lgs = paste0("LG",seq(1,21))
		 COLLECT.ALL.COMPARISONS=vector("list",10)
         OUT=data.frame()
         
         for (LG in lgs) {
            meta.data=table2[table2$CHR==LG,]
            for(SEX in c("MATERNAL","PATERNAL")){
            lg.length.cm.maternal=meta.data$CHRLENGTHCMMATERNAL
            lg.length.cm.paternal=meta.data$CHRLENGTHCMPATERNAL
            lg.length.cm=meta.data[[paste0("CHRLENGTHCM",SEX)]]
            centromere = meta.data$CENTROMERE
            centromere.cm = meta.data[[paste0("CENTROMERE",SEX,"CM")]]
		    lg.length = meta.data$LENGTH
		    shape = centromere/lg.length
		    q.arm.length=meta.data$QARMLENGTH
		    p.arm.length=meta.data$PARMLENGTH
		    q.arm.length.cm=meta.data[[paste0("QARMLENGTHCM",SEX)]]
		    p.arm.length.cm=meta.data[[paste0("PARMLENGTHCM",SEX)]]
		    d=subset(table1,CHR==LG)
		    if(shape>0.5){
		        d$MATERNALSITE1 = sapply(strsplit(d$MATERNALSITES,split=","),function(x) as.numeric(x)[order(as.numeric(x),decreasing=T)][1])
		        d$MATERNALSITE2 = sapply(strsplit(d$MATERNALSITES,split=","),function(x) as.numeric(x)[order(as.numeric(x),decreasing=T)][2])
	            d$MATERNALSITE3 = sapply(strsplit(d$MATERNALSITES,split=","),function(x) as.numeric(x)[order(as.numeric(x),decreasing=T)][3])
	            d$MATERNALSITE4 = sapply(strsplit(d$MATERNALSITES,split=","),function(x) as.numeric(x)[order(as.numeric(x),decreasing=T)][4])
	            d$MATERNALSITE5 = sapply(strsplit(d$MATERNALSITES,split=","),function(x) as.numeric(x)[order(as.numeric(x),decreasing=T)][5])
		        d$PATERNALSITE1 = sapply(strsplit(d$PATERNALSITES,split=","),function(x) as.numeric(x)[order(as.numeric(x),decreasing=T)][1])
		        d$PATERNALSITE2 = sapply(strsplit(d$PATERNALSITES,split=","),function(x) as.numeric(x)[order(as.numeric(x),decreasing=T)][2])
	            d$PATERNALSITE3 = sapply(strsplit(d$PATERNALSITES,split=","),function(x) as.numeric(x)[order(as.numeric(x),decreasing=T)][3])
	            d$PATERNALSITE4 = sapply(strsplit(d$PATERNALSITES,split=","),function(x) as.numeric(x)[order(as.numeric(x),decreasing=T)][4])
	            d$PATERNALSITE5 = sapply(strsplit(d$PATERNALSITES,split=","),function(x) as.numeric(x)[order(as.numeric(x),decreasing=T)][5])
		    }
		    d=subset(d,PRIORTRUECOCOUNT.EMMATERNAL>0.8 | PRIORTRUECOCOUNT.EMPATERNAL>0.8)
		    p.sites.all=as.numeric(unlist(strsplit(d$PATERNALSITES,split=",")))
	        p.sites.all=p.sites.all[is.na(p.sites.all)==F]
	        m.sites.all=as.numeric(unlist(strsplit(d$MATERNALSITES,split=",")))
	        m.sites.all=m.sites.all[is.na(m.sites.all)==F]
		    p.sites.all.cm=as.numeric(unlist(strsplit(d$PATERNALSITESCM,split=",")))
	        p.sites.all.cm=p.sites.all.cm[is.na(p.sites.all.cm)==F]
	        m.sites.all.cm=as.numeric(unlist(strsplit(d$MATERNALSITESCM,split=",")))
	        m.sites.all.cm=m.sites.all.cm[is.na(m.sites.all.cm)==F]
	        if(shape<=0.5){
	            p.sites.parm = p.sites.all[p.sites.all<=centromere]
	            p.sites.qarm = p.sites.all[p.sites.all>centromere]
	            m.sites.parm = m.sites.all[m.sites.all<=centromere]
	            m.sites.qarm = m.sites.all[m.sites.all>centromere]
	        } else {
	            p.sites.parm = p.sites.all[p.sites.all>=centromere]
	            p.sites.qarm = p.sites.all[p.sites.all<centromere]
	            m.sites.parm = m.sites.all[m.sites.all>=centromere]
	            m.sites.qarm = m.sites.all[m.sites.all<centromere]
	        }
	        
	       ############################################################ 
	        if(SEX=="MATERNAL"){
		        m=subset(d,MATERNALCOUNT==2 |MATERNALCOUNT==3|MATERNALCOUNT==4)
		        m$ParmCOUNT=0
		        for(i in 1:nrow(m)){
		            m[i,"ParmCOUNT"]=length(which(unlist(strsplit(m[i,"MATERNALSITES"],split=",")) %in% m.sites.parm))
		            m[i,"QarmCOUNT"]=length(which(unlist(strsplit(m[i,"MATERNALSITES"],split=",")) %in% m.sites.qarm))    	
			    }
		    } else{
		        m=subset(d,PATERNALCOUNT==2 |PATERNALCOUNT==3|PATERNALCOUNT==4)
		        m$ParmCOUNT=0
		        for(i in 1:nrow(m)){
		            m[i,"ParmCOUNT"]=length(which(unlist(strsplit(m[i,"PATERNALSITES"],split=",")) %in% p.sites.parm))
		            m[i,"QarmCOUNT"]=length(which(unlist(strsplit(m[i,"PATERNALSITES"],split=",")) %in% p.sites.qarm))    	
			    }
		    }
		    if(nrow(subset(m,ParmCOUNT==1))<4){next} 
		    CoClass=m[,paste0(SEX,"COUNT")]
	        IND=m$OFFSPRING
	        FEMALE=m$FEMALE
	        CoCountInParm=m$ParmCOUNT
	        CoCountInQarm=m$QarmCOUNT
	        SITE1=m[,paste0(SEX,"SITE1")]
	        SITE2=m[,paste0(SEX,"SITE2")]
	        SITE1CM=m[,paste0(SEX,"SITE1CM")]
	        SITE2CM=m[,paste0(SEX,"SITE2CM")]
	        DIST1TOC = abs(m[,paste0(SEX,"SITE1")]-centromere)
	        DIST1TOCCM = abs(m[,paste0(SEX,"SITE1CM")]-centromere.cm)
	        DIST1TOT = m[,paste0(SEX,"SITE1")]
	        if(shape>0.5){DIST1TOT = lg.length-m[,paste0(SEX,"SITE1")]}
	        DIST2TOC = abs(m[,paste0(SEX,"SITE2")]-centromere)
		    DIST3TOC = abs(m[,paste0(SEX,"SITE3")]-centromere)
		    DIST4TOC = abs(m[,paste0(SEX,"SITE4")]-centromere)
		    DIST1TO2 = abs(SITE1-SITE2)
		    
		    DIST2TOCCM = abs(m[,paste0(SEX,"SITE2CM")]-centromere.cm)
		    DIST3TOCCM = abs(m[,paste0(SEX,"SITE3CM")]-centromere.cm)
		    DIST4TOCCM = abs(m[,paste0(SEX,"SITE4CM")]-centromere.cm)
		    DIST1TO2CM = abs(m[,paste0(SEX,"SITE2CM")]-m[,paste0(SEX,"SITE1CM")])
		
	        
	        DIST1TOCPERLGLENGTH = DIST1TOC/lg.length
	        DIST2TOCPERLGLENGTH = DIST2TOC/lg.length
	        DIST3TOCPERLGLENGTH = DIST3TOC/lg.length
	        DIST4TOCPERLGLENGTH = DIST4TOC/lg.length
	        DIST1TO2PERLGLENGTH = DIST1TO2/lg.length
	        
	        DIST1TOCPERQARMLENGTH = DIST1TOC/q.arm.length
	        DIST1TOCPERPARMLENGTH = DIST1TOC/p.arm.length
	        DIST2TOCPERQARMLENGTH = DIST2TOC/q.arm.length
	        DIST2TOCPERPARMLENGTH = DIST2TOC/p.arm.length
	        DIST3TOCPERQARMLENGTH = DIST3TOC/q.arm.length
	        DIST4TOCPERQARMLENGTH = DIST4TOC/q.arm.length
	        DIST1TO2PERQARMLENGTH = DIST1TO2/q.arm.length
	        TRUECOCOUNT.V=m[,paste0("TRUECOCOUNT.V",SEX)]
	        PLOTTING.TABLE=data.frame(IND,FEMALE,CoClass,TRUECOCOUNT.V,CoCountInParm,CoCountInQarm,SITE1,SITE2,SITE1CM,SITE2CM,DIST1TOT,DIST1TOC,DIST2TOC,DIST3TOC,DIST4TOC,DIST1TO2CM,DIST1TOCCM,DIST2TOCCM,DIST3TOCCM,DIST4TOCCM,DIST1TO2CM,DIST1TOCPERLGLENGTH,
	        DIST2TOCPERLGLENGTH,DIST3TOCPERLGLENGTH,DIST4TOCPERLGLENGTH,DIST1TO2PERLGLENGTH,DIST1TOCPERQARMLENGTH,DIST2TOCPERQARMLENGTH,DIST3TOCPERQARMLENGTH,DIST4TOCPERQARMLENGTH,DIST1TO2PERQARMLENGTH,DIST1TOCPERPARMLENGTH,DIST2TOCPERPARMLENGTH)
	        
	       
		
		if(LG!="LG12"){
		OUT.tmp=subset(PLOTTING.TABLE,CoCountInParm==1)
		if(nrow(OUT.tmp)>=1){
			OUT.tmp$P.ARM.LENGTH=p.arm.length
			OUT.tmp$Q.ARM.LENGTH=q.arm.length
			OUT.tmp$LG.LENGTH=lg.length
			OUT.tmp$LG=LG
			OUT.tmp$P.ARM.LENGTH.CM=p.arm.length.cm
			OUT.tmp$Q.ARM.LENGTH.CM=q.arm.length.cm
			OUT.tmp$LG.LENGTH.CM=lg.length.cm
			OUT=rbind(OUT,OUT.tmp)
		}
	    table.out=OUT
	}
}
}
	"""
)
    
    
}
