#!/usr/bin/env anduril
import anduril.builtin._
import anduril.tools._
import anduril.anima._
import org.anduril.runtime._
//$OPT --threads 14


object jointFigures {

    val chromosomeMetadata9spn = INPUT("result_ninespinedRecombination/chromosomeMetadata/table.csv")
    val chromosomeMetadata3spn = INPUT("result_threespinedRecombination/chromosomeMetadata/table.csv")
    val crossoverDistributionTables9spn = INPUT("result_ninespinedRecombination/plots/outArray/shapePlots.csv")
    val crossoverDistributionTables3spn = INPUT("result_threespinedRecombination/plots/outArray/shapePlots.csv")
    val crossoverDistributionTables9spn2 = INPUT("result_ninespinedRecombination/plots/outArray/distributionPlots.csv")
    val crossoverDistributionTables3spn2 = INPUT("result_threespinedRecombination/plots/outArray/distributionPlots.csv")
    val betaRegressionTables9spn = INPUT("result_ninespinedRecombination/plots/outArray/lme.model.table.csv")
    val betaRegressionTables3spn = INPUT("result_threespinedRecombination/plots/outArray/lme.model.table.csv")
    val modelFits9spn = INPUT("result_ninespinedRecombination/fitModels/outArray/MODELFIT.csv")
    val modelFits3spn = INPUT("result_threespinedRecombination/fitModels/outArray/MODELFIT.csv")
    val modelEstimates9spn = INPUT("result_ninespinedRecombination/fitModels/outArray/ESTIMATES.csv")
    val modelEstimates3spn = INPUT("result_threespinedRecombination/fitModels/outArray/ESTIMATES.csv")
    val blindToCentromere9spn = INPUT("result_ninespinedRecombination/blindToCentromere/table.csv")
    val blindToCentromere3spn = INPUT("result_threespinedRecombination/blindToCentromere/table.csv")
    val concatenateCrossovers9spn = INPUT("result_ninespinedRecombination/concatenateRecombinationEventsHelsinki/table.csv")
    val concatenateCrossovers3spn = INPUT("result_threespinedRecombination/concatenateRecombinationEventsThreespine/table.csv")
    val cleanedDistributions9spn = INPUT("result_ninespinedRecombination/plotCleanedDistributions/outArray/empricalAndInferredDistributions.csv")
    val cleanedDistributions3spn = INPUT("result_threespinedRecombination/plotCleanedDistributions/outArray/empricalAndInferredDistributions.csv")
    
    val plotMetadata = REvaluate(
    table1 = chromosomeMetadata9spn,
    table2 = chromosomeMetadata3spn,
    script = """
     table.out=data.frame()
     rm('optOut1')
	 fig.dir <- get.output(cf, 'optOut1')
	 dir.create(fig.dir, recursive=TRUE)
	 setwd(fig.dir)
     library(ggplot2)
     library("gridExtra")
     get_legend<-function(myggplot){
	 tmp <- ggplot_gtable(ggplot_build(myggplot))
	 leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
	 legend <- tmp$grobs[[leg]]
	 return(legend)
     }
     ##
    
     Nspn.plot.table=table1[,c("CHR","LENGTH","MEANCOCOUNTMATERNAL")]
     colnames(Nspn.plot.table)[3] = "MEANCOCOUNTPATERNAL"
     Nspn.plot.table=rbind(Nspn.plot.table,table1[,c("CHR","LENGTH","MEANCOCOUNTPATERNAL")])
     colnames(Nspn.plot.table)[3] = "COCOUNT"
     Nspn.plot.table$SEX=rep(c("FEMALE","MALE"),each=21)
     Nspn.plot.table$SEXCODE=rep(c("SEX2","SEX1"),each=21)
     Nspn.plot.table$SPECIES="P. pungitius"
     Nspn.plot.table$COCOUNTPERbp = Nspn.plot.table$COCOUNT/Nspn.plot.table$LENGTH
     Tspn.plot.table=table2[,c("CHR","LENGTH","MEANCOCOUNTMATERNAL")]
     colnames(Tspn.plot.table)[3] = "MEANCOCOUNTPATERNAL"
     Tspn.plot.table=rbind(Tspn.plot.table,table2[,c("CHR","LENGTH","MEANCOCOUNTPATERNAL")])
     colnames(Tspn.plot.table)[3] = "COCOUNT"
     Tspn.plot.table$SEX=rep(c("FEMALE","MALE"),each=21)
     Tspn.plot.table$SEXCODE=rep(c("SEX2","SEX1"),each=21)
     Tspn.plot.table$SPECIES="G. aculeatus"
     Tspn.plot.table$COCOUNTPERbp = Tspn.plot.table$COCOUNT/Tspn.plot.table$LENGTH
     
     print("Fit first model, male as base line P.pun")
     model.9spn=lm(COCOUNT ~ LENGTH*SEXCODE,data=subset(Nspn.plot.table,CHR!="LG12"))
     print(summary(model.9spn))
     print("Fit first model, male as base line G.acu")
     model.3spn=lm(COCOUNT ~ LENGTH*SEXCODE,data=subset(Tspn.plot.table,CHR!="chrXIX"))
     print(summary(model.3spn))
     
     plot.table=rbind(Nspn.plot.table,Tspn.plot.table)
     plot.table$SPECIES=factor(plot.table$SPECIES,levels=c("P. pungitius","G. aculeatus"))
     print(coef(model.9spn))
     lm.table=data.frame(INTERCEPT = c(coef(model.9spn)[1],(coef(model.9spn)[1]+coef(model.9spn)[3]),coef(model.3spn)[1],(coef(model.3spn)[1]+coef(model.3spn)[3])),
     SLOPE= c(coef(model.9spn)[2],(coef(model.9spn)[2]+coef(model.9spn)[4]),coef(model.3spn)[2],(coef(model.3spn)[2]+coef(model.3spn)[4])),
     SEX=c("MALE","FEMALE","MALE","FEMALE") ,SPECIES=c("P. pungitius","P. pungitius","G. aculeatus","G. aculeatus"))
     lm.table$TEXT=paste0(lm.table$SPECIES,": y = ",signif(lm.table$INTERCEPT,2),"+",signif(lm.table$SLOPE,2),"*x")
     lm.table$SPECIES=factor(lm.table$SPECIES,levels=c("P. pungitius","G. aculeatus"))
     print(lm.table)
     plot.table$COL=paste0(plot.table$SEX,":",plot.table$SPECIES)
     plot.table$COL=factor(plot.table$COL, levels=c("FEMALE:P. pungitius","FEMALE:G. aculeatus","MALE:P. pungitius","MALE:G. aculeatus"))
     lm.table$COL=paste0(lm.table$SEX,":",lm.table$SPECIES)
     lm.table$COL=factor(lm.table$COL, levels=c("FEMALE:P. pungitius","FEMALE:G. aculeatus","MALE:P. pungitius","MALE:G. aculeatus"))
     p1=ggplot(plot.table,aes(x=LENGTH,y=COCOUNT,col=COL,fill=COL,shape=COL)) + 
     geom_abline(data=lm.table,inherit.aes=F,mapping=aes(intercept=INTERCEPT,slope=SLOPE,colour=COL,linetype=SPECIES),size=0.9) +
     geom_point(aes(shape=COL),color="black",size=1.5,stroke=0.2) + ylim(0.4,1.9) +
     ylab(label="Average number of crossovers per offspring") +xlab(label="Chromosome length (bp)")+
     scale_color_manual(values=c("tan2","tan3","steelblue2","steelblue3")) +
     scale_fill_manual(values=c("tan1","tan3","steelblue1","steelblue3"),labels = c("G. aculeatus female", "P. pungitius female","G. aculeatus male", "P. pungitius male")) +
     scale_shape_manual(values=c(21,24,21,24),labels=c("G. aculeatus female", "P. pungitius female","G. aculeatus male", "P. pungitius male")) +
     annotate("rect",xmin=14.5*10^6,xmax=26.8*10^6,ymin=1.46,ymax=1.85,fill="white") +
     annotate("text",size=2, hjust=0,x = 14.65*10^6, y = 1.8, label = "italic(`P. pun. female`) : y == 0.14 + 4.7 %*% 10 ^ -8 * x",parse = TRUE) +
     annotate("text",size=2, hjust=0,x = 14.65*10^6, y = 1.7, label = "italic(`G. acu. female`) : y == 0.22 + 1.6 %*% 10 ^ -8 * x",parse = TRUE) +
     annotate("text",size=2, hjust=0,x = 14.65*10^6, y = 1.6, label = "italic(`P. pun. male`) : y == 0.42 + 7.3 %*% 10 ^ -9 * x",parse = TRUE) +
     annotate("text",size=2, hjust=0,x = 14.65*10^6, y = 1.5, label = "italic(`G. acu. male`) : y == 0.13 + 4.1 %*% 10 ^ -8 * x",parse = TRUE) +
     theme(legend.position="none",text=element_text(size=6),axis.text=element_text(size=4),panel.grid.minor=element_blank())
     
     Nspn.plot.table$SEXCODE=rep(c("SEX1","SEX2"),each=21)
     Tspn.plot.table$SEXCODE=rep(c("SEX1","SEX2"),each=21)
     print(summary(lm(COCOUNTPERbp ~ LENGTH,data=subset(Nspn.plot.table,CHR!="LG12" & SEX=="FEMALE"))))
     print("Fit the second model, female as base line P.pun")
     model.9spn=lm(COCOUNTPERbp ~ LENGTH*SEXCODE,data=subset(Nspn.plot.table,CHR!="LG12"))
     print(summary(model.9spn))
     print("Fit the second model, female as base line G.acu")
     model.3spn=lm(COCOUNTPERbp ~ LENGTH*SEXCODE,data=Tspn.plot.table,CHR!="chrXIX")
     print(summary(model.3spn))
     print(coef(model.9spn))
     lm.table=data.frame(INTERCEPT = c(coef(model.9spn)[1],(coef(model.9spn)[1]+coef(model.9spn)[3]),coef(model.3spn)[1],(coef(model.3spn)[1]+coef(model.3spn)[3])),
     SLOPE= c(coef(model.9spn)[2],(coef(model.9spn)[2]+coef(model.9spn)[4]),coef(model.3spn)[2],(coef(model.3spn)[2]+coef(model.3spn)[4])),
     SEX=c("FEMALE","MALE","FEMALE","MALE") ,SPECIES=c("P. pungitius","P. pungitius","G. aculeatus","G. aculeatus"))
     lm.table$TEXT=paste0(lm.table$SPECIES,": y = ",signif(lm.table$INTERCEPT,2),signif(lm.table$SLOPE,2),"*x")
     lm.table$SPECIES=factor(lm.table$SPECIES,levels=c("P. pungitius","G. aculeatus"))
     lm.table$COL=paste0(lm.table$SEX,":",lm.table$SPECIES)
     lm.table$COL=factor(lm.table$COL, levels=c("FEMALE:P. pungitius","FEMALE:G. aculeatus","MALE:P. pungitius","MALE:G. aculeatus"))
     p2=ggplot(plot.table,aes(x=LENGTH,y=COCOUNTPERbp,col=COL,fill=COL,shape=COL)) + ylim(0.65*10^(-8),6.5*10^(-8)) +
     geom_abline(data=lm.table,inherit.aes=F,mapping=aes(intercept=INTERCEPT,slope=SLOPE,colour=COL,linetype=SPECIES),size=0.9,show.legend=F) +
     annotate("rect",xmin=14.5*10^6,xmax=29*10^6,ymin=0.655*10^(-8), ymax = 1.9*10^(-8),fill="white") +
     annotate("text",size=2, hjust=0,x = 14.65*10^6, y = 1.7*10^(-8), label = "italic(`P. pun. female`) : y == 6.0 %*% 10 ^ {-8} - 2.6 %*% 10 ^ {-16} * x",parse = TRUE) +
     annotate("text",size=2, hjust=0,x = 14.65*10^6, y = 1.4*10^(-8), label = "italic(`G. acu. female`) : y == 5.4 %*% 10 ^ {-8} - 2.8 %*% 10 ^ {-16} * x",parse = TRUE) +
     annotate("text",size=2, hjust=0,x = 14.65*10^6, y = 1.1*10^(-8), label = "italic(`P. pun. male`) : y == 4.6 %*% 10 ^ {-8} - 8.6 %*% 10 ^ {-16} * x",parse = TRUE) +
     annotate("text",size=2, hjust=0,x = 14.65*10^6, y = 0.8*10^(-8), label = "italic(`G. acu. male`) : y == 3.6 %*% 10 ^ {-8} - 4.3 %*% 10 ^ {-16} * x",parse = TRUE) +
     geom_point(aes(shape=COL),color="black",size=1.5,stroke=0.2) +
     ylab(label="Average number of crossovers per base pair")+xlab(label="Chromosome length (bp)") + 
     scale_color_manual(values=c("tan1","tan3","steelblue1","steelblue3"),labels = c("P. pungitius \n female","G. aculeatus \n female", "P. pungitius \n male","G. aculeatus \n male")) +
     scale_fill_manual(values=c("tan1","tan3","steelblue1","steelblue3"),labels = c("P. pungitius \n female","G. aculeatus \n female", "P. pungitius \n male","G. aculeatus \n male")) +
     scale_shape_manual(values=c(21,24,21,24),labels=c("P. pungitius \n female","G. aculeatus \n female", "P. pungitius \n male","G. aculeatus \n male")) +
     theme(text=element_text(size=6),axis.text=element_text(size=4),panel.grid.minor=element_blank(),legend.title=element_blank(),legend.position="bottom")
     
     legend <- get_legend(p2)
   	 p2=p2+theme(legend.position="none")
   	 png("regression3.png",width = 154, height = 75, units = "mm",res=1000)
     grid.arrange(p1, p2,legend,nrow=2, ncol=2,widths=c(2.3, 2.3),heights=c(2,0.2),layout_matrix = rbind(c(1, 2),c(3, 3))) 
     dev.off()
   	 pdf("regression3.pdf",width = 6, height = 3)
     grid.arrange(p1, p2,legend,nrow=2, ncol=2,widths=c(2.3, 2.3),heights=c(2,0.2),layout_matrix = rbind(c(1, 2),c(3, 3))) 
     dev.off()
     
     #########################################################################################
     #9-spine
     tmp=table1[order(table1$LENGTH),]
     maternal.tmp.ppun=data.frame(CHR=tmp$CHR,CHR.LENGTH=rep(tmp$LENGTH,10),MLE=c(tmp$MATERNALMLEp0,tmp$MATERNALMLEp1,tmp$MATERNALMLEp2,tmp$MATERNALMLEp3,tmp$MATERNALMLEp4,tmp$MATERNALMLEp5,tmp$MATERNALMLEp6,tmp$MATERNALMLEp7,tmp$MATERNALMLEp8,tmp$MATERNALMLEp9),CO.CLASS=rep(c("0","1","2","3","4","5","6","7","8","9"),each=21),SEX="Maternal",RESULT="Inferred",SPECIES="P. pungitius")
   	 paternal.tmp.ppun=data.frame(CHR=tmp$CHR,CHR.LENGTH=rep(tmp$LENGTH,10),MLE=c(tmp$PATERNALMLEp0,tmp$PATERNALMLEp1,tmp$PATERNALMLEp2,tmp$PATERNALMLEp3,tmp$PATERNALMLEp4,tmp$PATERNALMLEp5,tmp$PATERNALMLEp6,tmp$PATERNALMLEp7,tmp$PATERNALMLEp8,tmp$PATERNALMLEp9),CO.CLASS=rep(c("0","1","2","3","4","5","6","7","8","9"),each=21),SEX="Paternal",RESULT="Inferred",SPECIES="P. pungitius")   
     ppun.table=rbind(maternal.tmp.ppun,paternal.tmp.ppun)
     maternal.tmp.ppun=data.frame(CHR=tmp$CHR,CHR.LENGTH=rep(tmp$LENGTH,6),MLE=c(tmp$MATERNALn0,tmp$MATERNALn1,tmp$MATERNALn2,tmp$MATERNALn3,tmp$MATERNALn4,tmp$MATERNALn5)/sum(tmp[1,c(paste0("MATERNALn",0:5))]),CO.CLASS=rep(c("0","1","2","3","4","5"),each=21),SEX="Maternal",RESULT="Observed",SPECIES="P. pungitius")
   	 paternal.tmp.ppun=data.frame(CHR=tmp$CHR,CHR.LENGTH=rep(tmp$LENGTH,6),MLE=c(tmp$PATERNALn0,tmp$PATERNALn1,tmp$PATERNALn2,tmp$PATERNALn3,tmp$PATERNALn4,tmp$PATERNALn5)/sum(tmp[1,c(paste0("MATERNALn",0:5))]),CO.CLASS=rep(c("0","1","2","3","4","5"),each=21),SEX="Paternal",RESULT="Observed",SPECIES="P. pungitius")   
     ppun.table=rbind(ppun.table,maternal.tmp.ppun,paternal.tmp.ppun)
     ppun.table$RESULT=factor(ppun.table$RESULT,levels=c("Observed","Inferred"))     
     ppun.table$CHRNUMBER=as.numeric(sapply(strsplit(as.character(ppun.table$CHR),split="LG"), function(x) x[2]))
     ppun.table$CO.CLASS2=factor(ppun.table$CO.CLASS,levels=c("0","1","2","3","4","5","6","7","8","9","6-9"))
     ppun.table[ppun.table$CO.CLASS %in% c(6,7,8,9),"CO.CLASS2"]="6-9"
     
     #3-spine
     tmp=table2[order(table2$LENGTH),]
     maternal.tmp.gacu=data.frame(CHR=tmp$CHR,CHR.LENGTH=rep(tmp$LENGTH,10),MLE=c(tmp$MATERNALMLEp0,tmp$MATERNALMLEp1,tmp$MATERNALMLEp2,tmp$MATERNALMLEp3,tmp$MATERNALMLEp4,tmp$MATERNALMLEp5,tmp$MATERNALMLEp6,tmp$MATERNALMLEp7,tmp$MATERNALMLEp8,tmp$MATERNALMLEp9),CO.CLASS=rep(c("0","1","2","3","4","5","6","7","8","9"),each=21),SEX="Maternal",RESULT="Inferred",SPECIES="G. aculeatus")
   	 paternal.tmp.gacu=data.frame(CHR=tmp$CHR,CHR.LENGTH=rep(tmp$LENGTH,10),MLE=c(tmp$PATERNALMLEp0,tmp$PATERNALMLEp1,tmp$PATERNALMLEp2,tmp$PATERNALMLEp3,tmp$PATERNALMLEp4,tmp$PATERNALMLEp5,tmp$PATERNALMLEp6,tmp$PATERNALMLEp7,tmp$PATERNALMLEp8,tmp$PATERNALMLEp9),CO.CLASS=rep(c("0","1","2","3","4","5","6","7","8","9"),each=21),SEX="Paternal",RESULT="Inferred",SPECIES="G. aculeatus")
   	 gacu.table=rbind(maternal.tmp.gacu,paternal.tmp.gacu)
     maternal.tmp.gacu=data.frame(CHR=tmp$CHR,CHR.LENGTH=rep(tmp$LENGTH,6),MLE=c(tmp$MATERNALn0,tmp$MATERNALn1,tmp$MATERNALn2,tmp$MATERNALn3,tmp$MATERNALn4,tmp$MATERNALn5)/sum(tmp[1,c(paste0("MATERNALn",0:5))]),CO.CLASS=rep(c("0","1","2","3","4","5"),each=21),SEX="Maternal",RESULT="Observed",SPECIES="G. aculeatus")
   	 paternal.tmp.gacu=data.frame(CHR=tmp$CHR,CHR.LENGTH=rep(tmp$LENGTH,6),MLE=c(tmp$PATERNALn0,tmp$PATERNALn1,tmp$PATERNALn2,tmp$PATERNALn3,tmp$PATERNALn4,tmp$PATERNALn5)/sum(tmp[1,c(paste0("MATERNALn",0:5))]),CO.CLASS=rep(c("0","1","2","3","4","5"),each=21),SEX="Paternal",RESULT="Observed",SPECIES="G. aculeatus") 
     gacu.table=rbind(gacu.table,maternal.tmp.gacu,paternal.tmp.gacu)
     gacu.table$RESULT=factor(gacu.table$RESULT,levels=c("Observed","Inferred"))
     gacu.table$CHRNUMBER=as.numeric(as.roman(sapply(strsplit(as.character(gacu.table$CHR),split="chr"), function(x) x[2])))
     gacu.table$CO.CLASS2=factor(gacu.table$CO.CLASS,levels=c("0","1","2","3","4","5","6","7","8","9","6-9"))
     gacu.table[gacu.table$CO.CLASS %in% c(6,7,8,9),"CO.CLASS2"]="6-9"
     combined.plot.table=rbind(ppun.table,gacu.table)
     
     myPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","black","beige")

     p1=ggplot(ppun.table, aes(factor(CHR.LENGTH),fill=CO.CLASS2))+geom_bar(position="stack",aes(weight=MLE))+ scale_x_discrete("Chromosome", labels = ppun.table$CHRNUMBER)+ylab(label="Proportion") + scale_y_continuous(breaks=c(0,0.5,1)) +
   		  theme(axis.text.x=element_text(angle=60, hjust = 1,size=8),legend.position="none",strip.text.y=element_blank(),strip.text.x=element_text(face="italic")) + facet_grid(SEX+RESULT ~ SPECIES,scales="free_x") +
   		  scale_fill_manual(values=myPalette)
   		  
     p2=ggplot(gacu.table, aes(factor(CHR.LENGTH),fill=CO.CLASS2))+geom_bar(position="stack",aes(weight=MLE))+ scale_x_discrete("Chromosome", labels = gacu.table$CHRNUMBER) + labs(fill = "Crossover count")+
   		  theme(axis.text.x=element_text(angle = 60, hjust = 1,size=8),axis.text.y = element_blank(),axis.title.y = element_blank(),axis.ticks.y=element_blank(),legend.position="bottom",strip.text.x=element_text(face="italic")) + 
   		  facet_grid(SEX+RESULT ~ SPECIES,scales="free_x") + guides(fill = guide_legend(nrow = 1)) +
   		  scale_fill_manual(values=myPalette)
   		  
   	legend <- get_legend(p2)
   	p2=p2+theme(legend.position="none")
    png("barplot4.png",width = 200, height = 120, units = "mm",res=300)
   	grid.arrange(p1, p2, legend, ncol=2,nrow=2,widths=c(2.3, 2.3),heights=c(2,0.3),layout_matrix = rbind(c(1, 2),c(3, 3)))
    dev.off()
   	pdf("barplot4.pdf",width = 8, height = 4.5)
   	grid.arrange(p1, p2, legend, ncol=2,nrow=2,widths=c(2.3, 2.3),heights=c(2,0.3),layout_matrix = rbind(c(1, 2),c(3, 3)))
    dev.off()     
    
    ##ARM RATIO
    
    arm.ratio.ppun=data.frame(CHROMOSOME=table1$CHR,RATIO.bp=table1$PARMLENGTH/table1$QARMLENGTH,RATIO.cm=c(table1$PARMLENGTHCMMATERNAL/table1$QARMLENGTHCMMATERNAL,table1$PARMLENGTHCMPATERNAL/table1$QARMLENGTHCMPATERNAL,(table1$PARMLENGTHCMPATERNAL+table1$PARMLENGTHCMMATERNAL)/(table1$QARMLENGTHCMPATERNAL+table1$QARMLENGTHCMMATERNAL)),SEX=rep(c("MATERNAL","PATERNAL","AVERAGE"),each=21),SPECIES="P. pungitius")
    arm.ratio.gacu=data.frame(CHROMOSOME=table2$CHR,RATIO.bp=table2$PARMLENGTH/table2$QARMLENGTH,RATIO.cm=c(table2$PARMLENGTHCMMATERNAL/table2$QARMLENGTHCMMATERNAL,table2$PARMLENGTHCMPATERNAL/table2$QARMLENGTHCMPATERNAL,(table2$PARMLENGTHCMPATERNAL+table2$PARMLENGTHCMMATERNAL)/(table2$QARMLENGTHCMPATERNAL+table2$QARMLENGTHCMMATERNAL)),SEX=rep(c("MATERNAL","PATERNAL","AVERAGE"),each=21),SPECIES="G. aculeatus")
    
    
    plot.table=rbind(arm.ratio.ppun,arm.ratio.gacu)
    plot.table$SPECIES=factor(plot.table$SPECIES, levels=c("P. pungitius","G. aculeatus"))
    plot.table$SEX=factor(plot.table$SEX,labels=c("Maternal","Paternal","Average"), levels=c("MATERNAL","PATERNAL","AVERAGE"))
    plot.table$GROUP=paste(plot.table$SPECIES,plot.table$SEX)
    plot.table$GROUP=factor(plot.table$GROUP,levels=c(paste("P. pungitius",c("Maternal","Paternal","Average")),paste("G. aculeatus",c("Maternal","Paternal","Average"))))
    sex.chrom.annotation=subset(plot.table,CHROMOSOME %in% c("LG12","chrXIX"))
    ggplot(subset(plot.table,!(CHROMOSOME %in% c("LG12","chrXIX"))),aes(x=RATIO.bp,y=RATIO.cm,fill=GROUP,shape=GROUP)) +
    geom_abline(intercept=0,slope=1) + geom_point(aes(shape=GROUP,fill=GROUP),color="black",size=2,stroke=0.2) +
    ylab(label="p(cM)/q(cM)") +xlab(label="p(bp)/q(bp)")+ ylim(0,1) +xlim(0,1) +
    scale_fill_manual(values=c("tan1","steelblue1","grey40","tan3","steelblue3","grey50"))+
    scale_shape_manual(values=c(21,21,21,24,24,24))+
    geom_text(data=sex.chrom.annotation,label="sex chrom.",hjust = 0,vjust=1,nudge_y=-0.02,nudge_x=+0.02,size=2) +
    geom_point(data=sex.chrom.annotation,shape=22,fill="firebrick",col="firebrick",show.legend=F) +
    facet_grid(SPECIES ~SEX) +
    theme(legend.position="bottom",legend.title=element_blank(),strip.text.y=element_text(face="italic"),panel.grid.minor=element_blank(),aspect.ratio=1)
    ggsave("ratio.png",device="png")
    ggsave("ratio.pdf",device="pdf",width=170,height=130,units="mm")
    
    """
    )
    
    val plotDistributions=REvaluate(
        table1=crossoverDistributionTables9spn,
        table2=crossoverDistributionTables3spn,
        table3=crossoverDistributionTables9spn2,
        table4=crossoverDistributionTables3spn2,
        script="""
          table.out=data.frame()
	      rm('optOut1')
		  fig.dir <- get.output(cf, 'optOut1')
		  dir.create(fig.dir, recursive=TRUE)
		  setwd(fig.dir)
	      library(ggplot2)
	      library("gridExtra")
	      get_legend<-function(myggplot){
		  tmp <- ggplot_gtable(ggplot_build(myggplot))
		  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
		  legend <- tmp$grobs[[leg]]
		  return(legend)
	      }
	      ##
          myPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","black","beige")[2:4]
          table1$SPECIES="P. pungitius"
          table2$SPECIES="G. aculeatus"     
          plot.table=rbind(table1,table2)
          plot.table$COLOR=factor(plot.table$COLOR,levels=c("MATERNAL","PATERNAL","MEDIAN"))
          plot.table$COLOR=factor(plot.table$COLOR,labels=c("Maternal","Paternal","Median"))
          plot.table$SEX=factor(plot.table$SEX,levels=c("MATERNAL","PATERNAL"))
          plot.table$SEX=factor(plot.table$SEX,labels=c("Maternal","Paternal"))
          plot.table$SHAPE=factor(plot.table$SHAPE,levels=c("METACENTRIC","SUBMETACENTRIC","SUBTELOCENTRIC","ACROCENTRIC","TELOCENTRIC"))
          plot.table$SHAPE=factor(plot.table$SHAPE,labels=c("Metacentric","Submetacentric","Subtelocentric","Acrocentric","Telocentric"))
          plot.table$CENTROMERE = c(13,9,6,3,1)
          plot.table$CENTROMERE2 = c(9,6,3,1,1)
          plot.table[plot.table$SHAPE=="Metacentric","CENTROMERE"]=13
          plot.table[plot.table$SHAPE=="Submetacentric","CENTROMERE"]=9
          plot.table[plot.table$SHAPE=="Subtelocentric","CENTROMERE"]=6
          plot.table[plot.table$SHAPE=="Acrocentric","CENTROMERE"]=3
          plot.table[plot.table$SHAPE=="Telocentric","CENTROMERE"]=1
          plot.table[plot.table$SHAPE=="Metacentric","CENTROMERE2"]=9
          plot.table[plot.table$SHAPE=="Submetacentric","CENTROMERE2"]=6
          plot.table[plot.table$SHAPE=="Subtelocentric","CENTROMERE2"]=3
          plot.table[plot.table$SHAPE=="Acrocentric","CENTROMERE2"]=1
          plot.table[plot.table$SHAPE=="Telocentric","CENTROMERE2"]=1
          plot.table$SPECIES=factor(plot.table$SPECIES,levels=c("P. pungitius","G. aculeatus"))
          
          plot.table.1=subset(plot.table,QARMCOUNT=="all")
          plot.table.2=subset(plot.table,QARMCOUNT!="all")
          p1=ggplot(plot.table.1,aes(x=BIN,y=co,color=COLOR,group=CHR))+
          geom_rect(aes(xmin=CENTROMERE,xmax=CENTROMERE2,ymin=0,ymax=0.30),color= "grey60", fill="grey70",inherit.aes=F) +
          geom_line(aes(size=COLOR)) + scale_size_manual(values=c(0.2,0.2,0.3)) +
          facet_grid(SHAPE ~ SEX+SPECIES) + scale_color_manual(values=c("orange", "skyblue2", "purple")) +
          scale_x_continuous(breaks=c(1,25))+ ylab(label="Density") + xlab("Chromosome bin") +
          annotate("segment",x=13, xend=13,y=0.30,yend=0.24,arrow=arrow(length = unit(0.15, "cm"))) +
          theme(strip.text.y=element_text(size=5),strip.text.x=element_text(size=5),text=element_text(size=6),legend.position="bottom",legend.title=element_blank(),panel.grid=element_blank())
          
          p2=ggplot(subset(plot.table.2,QARMCOUNT==2),aes(x=BIN,y=co,color=COLOR,group=CHR))+
          geom_rect(aes(xmin=CENTROMERE,xmax=CENTROMERE2,ymin=0,ymax=0.15),color= "grey60", fill="grey70",inherit.aes=F) +
          geom_line(aes(size=COLOR))+scale_size_manual(values=c(0.2,0.3)) +
          facet_grid(SHAPE ~ SEX+SPECIES)+ scale_color_manual(values=c("orange","purple")) +
          scale_x_continuous(breaks=c(1,25)) + ylab(label="Density") + xlab("Chromosome bin") +
          theme(panel.grid=element_blank(),legend.position="none") + 
          annotate("segment",x=13, xend=13,y=0.15,yend=0.12,arrow=arrow(length = unit(0.15, "cm"))) +
          theme(strip.text.y=element_text(size=5),strip.text.x=element_text(size=5),text=element_text(size=6))
          
          table3$SPECIES="P. pungitius"
          table4$SPECIES="G. aculeatus"
          plot.table=rbind(table3,table4)
          plot.table$SPECIES=factor(plot.table$SPECIES,levels=c("P. pungitius","G. aculeatus"))
          plot.table$SEX=factor(plot.table$SEX,labels=c("Maternal","Paternal"))
          plot.table$CHR=factor(plot.table$CHR,levels=c(paste0("LG",1:21),paste0("chr",as.roman(1:21))))
          plot.table$CENTROMERE=rep(c(8,6,21,18,15,24,3,11,10,17,18,14,10,0,24,17,16,25,13,14,17,17,1,5,17,3,3,11,0,5,1,2,23,2,21,21,0,15,7,5,23,6),each=9*25)
          plot.table$CO.GROUP=factor(plot.table$CO.GROUP,labels=c("1 crossover","2 crossovers","3 crossovers"))
          plot.chromosomes=c(paste0("LG",c(1,2,4)),paste0("chr",as.roman(c(1,4,7))))
          
          p=ggplot(subset(plot.table,CHR %in% c("LG1","chrI")),aes(x=BIN,y=DENSITY,fill=as.factor(CO),col=as.factor(CO)))+
          geom_segment(aes(x = CENTROMERE,xend = CENTROMERE),inherit.aes=F, y = 0, yend = 0.3, colour = "grey60", fill="grey70") +  
          geom_line()+ geom_ribbon(aes(x=BIN,ymax=DENSITY),ymin=0,alpha=0.3) +
          scale_color_manual("Ordinal number of crossover",values=myPalette)+
          scale_fill_manual("Ordinal number of crossover",values=myPalette)+ 
          facet_grid(CHR+SPECIES ~ SEX + CO.GROUP) +
          ylab("Density") + xlab("Chromosome bin") + scale_x_continuous(breaks=c(1,25)) + scale_y_continuous(breaks=c(0,0.3)) +
          annotate("segment",x=13, xend=13,y=0.35,yend=0.25,arrow=arrow(length = unit(0.15, "cm"))) +
          theme(legend.position="top",legend.key.size=unit(5,units="mm"),panel.grid=element_blank(),strip.text.y=element_text(size=6),strip.text.x=element_text(size=6.6),text=element_text(size=6))
          
          png("combined2.png",height=185,width=120,units="mm",res=300)
          grid.arrange(p,p1,p2,layout_matrix=rbind(c(1,1),c(2,3)),heights=c(1,1.9),widths=c(3,2),nrow=2,ncol=2)  
          dev.off()
          
          pdf("combined2.pdf",height=7.28,width=4.64,onefile=F)
          legend=get_legend(p1)
          p1=p1+theme(legend.position="none")
          grid.arrange(p,p1,p2,legend,layout_matrix=rbind(c(1,1),c(2,3),c(2,4)),heights=c(1,1.9,0.15),widths=c(3,2),nrow=3,ncol=2)  
          dev.off()
          
          p=ggplot(subset(plot.table,SPECIES=="P. pungitius"),aes(x=BIN,y=DENSITY,fill=as.factor(CO),col=as.factor(CO)))+
          geom_line()+ geom_ribbon(aes(x=BIN,ymax=DENSITY),ymin=0,alpha=0.3) +
          geom_segment(aes(x = CENTROMERE,xend = CENTROMERE),inherit.aes=F, y = 0, yend = 0.3, colour = "white", fill="white",) +  
          scale_color_manual("Ordinal number of crossover",values=myPalette)+
          scale_fill_manual("Ordinal number of crossover",values=myPalette)+ 
          facet_grid(CHR ~ SEX + CO.GROUP,scales="free_y") +
          ylab("Density") + xlab("Chromosome bin") + scale_x_continuous(breaks=c(1,25)) + scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1.0)) +
          annotate("segment",x=13, xend=13,y=0.35,yend=0.25,arrow=arrow(length = unit(0.05, "cm"))) +
          theme(legend.position="top",legend.key.size=unit(5,units="mm"),panel.grid=element_blank(),strip.text.y=element_text(size=4.5),strip.text.x=element_text(size=4.5),text=element_text(size=4.5))
          ggsave("supplementppun.pdf",device="pdf",plot=p,width=80,height=213,units="mm")

          p=ggplot(subset(plot.table,SPECIES=="G. aculeatus"),aes(x=BIN,y=DENSITY,fill=as.factor(CO),col=as.factor(CO)))+
          geom_line()+ geom_ribbon(aes(x=BIN,ymax=DENSITY),ymin=0,alpha=0.3) +
          geom_segment(aes(x = CENTROMERE,xend = CENTROMERE),inherit.aes=F, y = 0, yend = 0.3, colour = "white", fill="white",) +  
          scale_color_manual("Ordinal number of crossover",values=myPalette)+
          scale_fill_manual("Ordinal number of crossover",values=myPalette)+ 
          facet_grid(CHR ~ SEX + CO.GROUP,scales="free_y") +
          ylab("Density") + xlab("Chromosome bin") + scale_x_continuous(breaks=c(1,25)) + scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1.0)) +
          annotate("segment",x=13, xend=13,y=0.35,yend=0.25,arrow=arrow(length = unit(0.05, "cm"))) +
          theme(legend.position="top",legend.key.size=unit(5,units="mm"),panel.grid=element_blank(),strip.text.y=element_text(size=4.5),strip.text.x=element_text(size=4.5),text=element_text(size=4.5))
          ggsave("supplementgacu.pdf",device="pdf",plot=p,width=80,height=213,units="mm")
        """
    )
    val plotBetaRegression=REvaluate(
        table1=betaRegressionTables9spn,
        table2=betaRegressionTables3spn,
        script="""
          table.out=data.frame()
	      rm('optOut1')
		  fig.dir <- get.output(cf, 'optOut1')
		  dir.create(fig.dir, recursive=TRUE)
		  setwd(fig.dir)
	      library(ggplot2)
	      library("gridExtra")
	      library(lme4)
          library(glmmTMB)
	      ##
          table1$SPECIES="P. pungitius"
          table2$SPECIES="G. aculeatus"
          table2=table2[,colnames(table1)]
          print(head(table2))
          plot.table=rbind(table1,table2)
          plot.table$SPECIES=factor(plot.table$SPECIES,levels=c("P. pungitius","G. aculeatus"))
          plot.table$SEX=factor(plot.table$SEX,labels=c("Maternal","Paternal"))
          p=ggplot(plot.table, aes(as.factor(QARMCOCOUNT),DISTANCETOCENTROMERE)) + 
          geom_boxplot(aes(fill=SEX),position = position_dodge(preserve = "single")) + 
          xlab("No. Crossovers in q-arm") + ylab("Distance to centromere") + 
          facet_grid(~SPECIES) + scale_fill_manual(values=c("tan2", "steelblue2")) +
          theme(panel.grid=element_blank(),legend.title=element_blank(),legend.position="bottom",strip.text.x=element_text(face="italic"))   
		  ggsave("dist2centromere1.png", plot = last_plot(), device = "png", scale = 1, width = 150, height = 75, units = "mm", dpi = 300)
		  ggsave("dist2centromere1.pdf", plot = last_plot(), device = "pdf", scale = 1, width = 150, height = 75, units = "mm")
		  
		  table1$SEX=factor(table1$SEX,levels=c("PATERNAL","MATERNAL"))
		  m1.ppun=glmmTMB(DISTANCETOCENTROMERE ~ fQARMCOCOUNT + SEX + (1|fOFFSPRING) + (1|fCHR), data=table1, family=beta_family) 
          m2.ppun=glmmTMB(DISTANCETOCENTROMERE ~ fQARMCOCOUNT*SEX + (1|fOFFSPRING) + (1|fCHR), data=table1, family=beta_family)
		  print(summary(m1.ppun))
		  print(summary(m2.ppun))
		  table2$SEX=factor(table2$SEX,levels=c("PATERNAL","MATERNAL"))
		  m1.gacu=glmmTMB(DISTANCETOCENTROMERE ~ fQARMCOCOUNT + SEX + (1|fOFFSPRING) + (1|fCHR), data=table2, family=beta_family) 
		  m2.gacu=glmmTMB(DISTANCETOCENTROMERE ~ fQARMCOCOUNT*SEX + (1|fOFFSPRING) + (1|fCHR), data=table2, family=beta_family)
		  print(summary(m1.gacu))
		  print(summary(m2.gacu))
        """
    )
    
    val plotModelFits=REvaluate(
        table1=modelEstimates9spn,
        table2=modelEstimates3spn,
        table3=modelFits9spn,
        table4=modelFits3spn,
        script="""
          table.out=data.frame()
	      rm('optOut1')
		  fig.dir <- get.output(cf, 'optOut1')
		  dir.create(fig.dir, recursive=TRUE)
		  setwd(fig.dir)
	      library(ggplot2)
	      library("gridExtra")
	      
	      ppun.maternal=data.frame(CHR=table1$CHR,SEX="MATERNAL",ESTIMATE.TYPE=rep(c("vStahl","pStahl","vGamma"),each=21),
	      VALUE=c(table1$MATERNALvStahl,table1$MATERNALpStahl,table1$MATERNALvGamma),
	      UPPER=c(table1$MATERNALvStahl+2*table1$MATERNALvStahlSd,table1$MATERNALpStahl+2*table1$MATERNALpStahlSd,table1$MATERNALvGamma+2*table1$MATERNALvGammaSd),
	      LOWER=c(table1$MATERNALvStahl-2*table1$MATERNALvStahlSd,table1$MATERNALpStahl-2*table1$MATERNALpStahlSd,table1$MATERNALvGamma-2*table1$MATERNALvGammaSd)
	      )
	      ppun.paternal=data.frame(CHR=table1$CHR,SEX="PATERNAL",ESTIMATE.TYPE=rep(c("vStahl","pStahl","vGamma"),each=21),
	      VALUE=c(log10(table1$PATERNALvStahl),table1$PATERNALpStahl,table1$PATERNALvGamma),
	      UPPER=c(log10(table1$PATERNALvStahl+2*table1$PATERNALvStahlSd),table1$PATERNALpStahl+2*table1$PATERNALpStahlSd,table1$PATERNALvGamma+2*table1$PATERNALvGammaSd),
	      LOWER=c(log10(table1$PATERNALvStahl-2*table1$PATERNALvStahlSd),table1$PATERNALpStahl-2*table1$PATERNALpStahlSd,table1$PATERNALvGamma-2*table1$PATERNALvGammaSd)
	      )
	      ppun.plot.table=rbind(ppun.maternal,ppun.paternal)
	      ppun.plot.table$SPECIES="PPUN"
	      ppun.plot.table$CHR=factor(ppun.plot.table$CHR,levels=paste0("LG",1:21))
	      ppun.plot.table$ESTIMATE.TYPE=factor(ppun.plot.table$ESTIMATE.TYPE,labels=c("v gamma-sprinkling","p gamma-sprinkling","v gamma"),levels=c("vStahl","pStahl","vGamma"))
	      
	      gacu.maternal=data.frame(CHR=table2$CHR,SEX="MATERNAL",ESTIMATE.TYPE=rep(c("vStahl","pStahl","vGamma"),each=21),
	      VALUE=c(table2$MATERNALvStahl,table2$MATERNALpStahl,table2$MATERNALvGamma),
	      UPPER=c(table2$MATERNALvStahl+2*table2$MATERNALvStahlSd,table2$MATERNALpStahl+2*table2$MATERNALpStahlSd,table2$MATERNALvGamma+2*table2$MATERNALvGammaSd),
	      LOWER=c(table2$MATERNALvStahl-2*table2$MATERNALvStahlSd,table2$MATERNALpStahl-2*table2$MATERNALpStahlSd,table2$MATERNALvGamma-2*table2$MATERNALvGammaSd)
	      )
	      gacu.paternal=data.frame(CHR=table2$CHR,SEX="PATERNAL",ESTIMATE.TYPE=rep(c("vStahl","pStahl","vGamma"),each=21),
	      VALUE=c(log10(table2$PATERNALvStahl),table2$PATERNALpStahl,table2$PATERNALvGamma),
	      UPPER=c(log10(table2$PATERNALvStahl+2*table2$PATERNALvStahlSd),table2$PATERNALpStahl+2*table2$PATERNALpStahlSd,table2$PATERNALvGamma+2*table2$PATERNALvGammaSd),
	      LOWER=c(log10(table2$PATERNALvStahl-2*table2$PATERNALvStahlSd),table2$PATERNALpStahl-2*table2$PATERNALpStahlSd,table2$PATERNALvGamma-2*table2$PATERNALvGammaSd)
	      )
	      gacu.plot.table=rbind(gacu.maternal,gacu.paternal)
	      gacu.plot.table$SPECIES="GACU"
	      gacu.plot.table$CHR=factor(gacu.plot.table$CHR,levels=paste0("chr",as.roman(1:21)))
	      gacu.plot.table$ESTIMATE.TYPE=factor(gacu.plot.table$ESTIMATE.TYPE,labels=c("v gamma-sprinkling","p gamma-sprinkling","v gamma"),levels=c("vStahl","pStahl","vGamma"))
	      
	      plot.table=rbind(ppun.plot.table,gacu.plot.table)
	      plot.table$SPECIES = factor(plot.table$SPECIES,labels=c("P. pungitius","G. aculeatus"),levels=c("PPUN","GACU"))
	      plot.table$SEX = factor(plot.table$SEX,labels=c("Maternal","Paternal"),levels=c("MATERNAL","PATERNAL"))
	      plot.table$COL=paste0(plot.table$SEX,":",plot.table$SPECIES)
	      plot.table$CHRN=1:21
	      p1=ggplot(plot.table,aes(CHRN,VALUE)) + geom_point(aes(shape=SPECIES,colour=COL),size=1.5) + 
	      geom_errorbar(aes(ymin=LOWER,ymax=UPPER),width=0.2)+
	      facet_grid(ESTIMATE.TYPE ~ SEX+SPECIES,scales="free") + xlab(label="Chromosome") + ylab(label="Estimate") +
	      theme(axis.text.x=element_text(size=5),legend.position="none",panel.grid=element_blank()) + 
	      scale_color_manual(values=c("tan2","tan3","steelblue2","steelblue3"))
	      ggsave("estimates1.png")
	      ggsave("estimates1.pdf",device="pdf",height=140,width=117,units="mm")
	      
	      print(colnames(table3))
	      print(colnames(table4))
	      plot.table=data.frame(CHR=1:21,PVALUE=c(table3$MATERNALpStahl,table3$MATERNALpGamma,table3$MATERNALpYF,table3$PATERNALpStahl,table3$PATERNALpGamma,table3$PATERNALpYF,table4$MATERNALpStahl,table4$MATERNALpGamma,table4$MATERNALpYF,table4$PATERNALpStahl,table4$PATERNALpGamma,table4$PATERNALpYF),MODEL=rep(c("Gamma-sprinkling","Gamma","Yu & Feingold"),each=21),SEX=rep(c("Maternal","Paternal"),each=3*21),SPECIES=rep(c("P. pungitius","G. aculeatus"),each=6*21),
	      SIGN=c("1","2"))
	      plot.table$SPECIES = factor(plot.table$SPECIES,levels=c("P. pungitius", "G. aculeatus"))
	      print(plot.table)
	      ggplot(plot.table,aes(x=CHR,y=PVALUE)) + geom_point(size=0.5) + facet_grid(MODEL~SPECIES+SEX) +
	       scale_color_manual(labels = c("p<0.05","p>=0.05"), values = c("salmon","skyblue3")) + 
	       ylab("p-value") + xlab("Chromosome") + geom_hline(yintercept=0.05,lty=2) +
	       theme(legend.title=element_blank(),legend.position="bottom",panel.grid=element_blank(),text=element_text(size=5))
	      ggsave("p-values.png",width=100,height=130,units="mm",dpi=300)
	      
	      
        """
        )
        
     val plotBlindToCentromere=REvaluate(
         table1=blindToCentromere9spn,
         table2=blindToCentromere3spn,
         script="""
          table.out=data.frame()
	      rm('optOut1')
		  fig.dir <- get.output(cf, 'optOut1')
		  dir.create(fig.dir, recursive=TRUE)
		  setwd(fig.dir)
	      library(ggplot2)
	      library("gridExtra")
          
          ##P.pun 
          d.vector=seq(5,100,5)
		  cor.estimates=c()
		  upper.conf=c()
		  lower.conf=c()
		  used.distances=c()
		  for(d in d.vector){
		    tmp=table1[table1$DIST1TOCCM<=d & table1$DIST2TOCCM<=d &  table1$P.ARM.LENGTH.CM>=d & table1$Q.ARM.LENGTH.CM >= d & table1$CoCountInParm==1,]
	        if(nrow(tmp)>=4){
	            COR=cor.test(x=tmp$DIST2TOCCM,y=tmp$DIST1TOCCM,method="pearson")
	            cor.estimates=c(cor.estimates,as.numeric(COR$estimate))
	            lower.conf=c(lower.conf,as.numeric(COR$conf.int[1]))
	            upper.conf=c(upper.conf,as.numeric(COR$conf.int[2]))
	            used.distances=c(used.distances,d)
			}
	      }
	      ppun.plot.table=data.frame(VALUE=c(cor.estimates,lower.conf,upper.conf),D=used.distances,GROUP=rep(c("ESTIMATE","LOWER","UPPER"),each=length(used.distances)),LTYPE=rep(c("ESTIMATE","CONFIDENCE","CONFIDENCE"),each=length(used.distances)),SPECIES="PPUN")  
          ppun.plot.table$LTYPE=factor(ppun.plot.table$LTYPE,levels=c("ESTIMATE","CONFIDENCE"))
          ##G.acu
          d.vector=seq(5,100,5)
		  cor.estimates=c()
		  upper.conf=c()
		  lower.conf=c()
		  used.distances=c()
		  for(d in d.vector){
		    tmp=table2[table2$DIST1TOCCM<=d & table2$DIST2TOCCM<=d &  table2$P.ARM.LENGTH.CM>=d & table2$Q.ARM.LENGTH.CM >= d & table2$CoCountInParm==1,]
	        if(nrow(tmp)>=4){
	            COR=cor.test(x=tmp$DIST2TOCCM,y=tmp$DIST1TOCCM,method="pearson")
	            cor.estimates=c(cor.estimates,as.numeric(COR$estimate))
	            lower.conf=c(lower.conf,as.numeric(COR$conf.int[1]))
	            upper.conf=c(upper.conf,as.numeric(COR$conf.int[2]))
	            used.distances=c(used.distances,d)
			}
	    }
	    gacu.plot.table=data.frame(VALUE=c(cor.estimates,lower.conf,upper.conf),D=used.distances,GROUP=rep(c("ESTIMATE","LOWER","UPPER"),each=length(used.distances)),LTYPE=rep(c("ESTIMATE","CONFIDENCE","CONFIDENCE"),each=length(used.distances)),SPECIES="GACU")  
        gacu.plot.table$LTYPE=factor(gacu.plot.table$LTYPE,levels=c("ESTIMATE","CONFIDENCE"))
        
        plot.table=rbind(ppun.plot.table,gacu.plot.table)
        plot.table$SPECIES=factor(plot.table$SPECIES,labels=c("P. pungitius","G. aculeatus"),levels=c("PPUN","GACU"))
        p=ggplot(plot.table,aes(x=D,y=VALUE,group=GROUP)) + geom_hline(yintercept=0,color="firebrick") +
        geom_line(aes(linetype=LTYPE)) + facet_grid(~SPECIES) +
        ylab(label="Correlation coefficient") + 
        xlab(label="Maximum distance of crossovers from the centromere (cM)") + 
        theme(legend.position="none",panel.grid.minor=element_blank(),text=element_text(size=6),strip.text.x=element_text(face="italic",size=6))
        ggsave("blindToCentromere.png")    
        ggsave("blindToCentromere.pdf",device="pdf",height=40,width=85,units="mm")    
         """
     )
     
     val plotDistanceToCentromere=REvaluate( 
         table1=concatenateCrossovers9spn,
         table2=concatenateCrossovers3spn,
         table3=chromosomeMetadata9spn,
         table4=chromosomeMetadata3spn,
         script="""
          table.out=data.frame()
	      rm('optOut1')
		  fig.dir <- get.output(cf, 'optOut1')
		  dir.create(fig.dir, recursive=TRUE)
		  setwd(fig.dir)
	      library(ggplot2)
	      library("gridExtra")
          ##P.pun 
          mean.distance.to.centromere.male=c()
		  mean.distance.to.centromere.female=c()
		  lgs=paste0("LG",1:21)
		  for(LG in lgs){
		    d=subset(table1,CHR==LG)
		    centromere=table3[table3$CHR==LG,"CENTROMERE"]
			mean.distance.to.centromere.male=c(mean.distance.to.centromere.male,mean(abs(centromere-as.numeric(unlist(strsplit(d$PATERNALSITES,split=","))))))
			mean.distance.to.centromere.female=c(mean.distance.to.centromere.female,mean(abs(centromere-as.numeric(unlist(strsplit(d$MATERNALSITES,split=","))))))
	      }
	      ppun.table=data.frame(CHR=lgs,CHR2=rep(c("autosome","sex","autosome"),times=c(11,1,9)),
	      PATERNAL=mean.distance.to.centromere.male,MATERNAL=mean.distance.to.centromere.female,SPECIES="PPUN")
          ##G.acu 
          mean.distance.to.centromere.male=c()
		  mean.distance.to.centromere.female=c()
		  chrs=paste0("chr",as.roman(1:21))
		  for(CHROM in chrs){
		    d=subset(table2,CHR==CHROM)
		    centromere=table4[table4$CHR==CHROM,"CENTROMERE"]
			mean.distance.to.centromere.male=c(mean.distance.to.centromere.male,mean(abs(centromere-as.numeric(unlist(strsplit(d$PATERNALSITES,split=","))))))
			mean.distance.to.centromere.female=c(mean.distance.to.centromere.female,mean(abs(centromere-as.numeric(unlist(strsplit(d$MATERNALSITES,split=","))))))
	      }
	      gacu.table=data.frame(CHR=chrs,CHR2=rep(c("autosome","sex","autosome"),times=c(18,1,2)),
	      PATERNAL=mean.distance.to.centromere.male,MATERNAL=mean.distance.to.centromere.female,SPECIES="GACU")
	      ##
	      plot.table=rbind(ppun.table,gacu.table)
	      plot.table$SPECIES=factor(plot.table$SPECIES,labels=c("P. pungitius","G. aculeatus"),levels=c("PPUN","GACU"))
	      plot.table$MATERNAL = plot.table$MATERNAL/10^6
	      plot.table$PATERNAL = plot.table$PATERNAL/10^6
	      p=ggplot(plot.table,aes(x=MATERNAL,y=PATERNAL)) + geom_point(aes(shape=CHR2,col=CHR2),show.legend=F) + geom_abline(intercept=0,slope=1) + 
	      xlim(range(c(plot.table$PATERNAL,plot.table$MATERNAL))) + ylim(range(c(plot.table$PATERNAL,plot.table$MATERNAL))) +
	      xlab("Maternal")+ylab("Paternal")+labs(title="Crossover distance from centromere (Mbp)")+coord_fixed(ratio=1) + facet_grid(~SPECIES) + 
	      ggsave("distanceToCentromere.png")
	      
	      p=ggplot(plot.table,aes(x=MATERNAL,y=PATERNAL,color=SPECIES)) + geom_point(aes(shape=CHR2),show.legend=T) + 
	      geom_abline(intercept=0,slope=1,size=0.3) + 
	      xlim(range(c(plot.table$PATERNAL,plot.table$MATERNAL))) + ylim(range(c(plot.table$PATERNAL,plot.table$MATERNAL))) +
	      xlab("Maternal")+ylab("Paternal")+labs(title="Crossover distance from centromere (Mbp)")+
	      annotate("text",x=plot.table[21+19,"MATERNAL"]+2.5,y=plot.table[21+19,"PATERNAL"],label="G. aculeatus sex chromosome",size=1.5)+
	      annotate("text",x=plot.table[12,"MATERNAL"]+2.5,y=plot.table[12,"PATERNAL"],label="P. pungitius sex chromosome",size=1.5)+
	      guides(shape = FALSE) + labs(color = "Species") +
	      theme(text=element_text(size=6),legend.position=c(0.8,y=0.15),panel.grid.minor=element_blank())
	      ggsave("distanceToCentromere2.png",width=85,height=85,units="mm",dpi=1000)
	      ggsave("distanceToCentromere2.pdf",width=85,height=85,units="mm")
	      
	      
         """
     )
     
     val plotSimulatedDistributions = REvaluate(
     table1=cleanedDistributions9spn,
     table2=cleanedDistributions3spn,
     script="""
          table.out=data.frame()
	      rm('optOut1')
		  fig.dir <- get.output(cf, 'optOut1')
		  dir.create(fig.dir, recursive=TRUE)
		  setwd(fig.dir)
	      library(ggplot2)
	      library("gridExtra")
          ##P.pun
          use.chromosomes=paste0("LG",c(7,10))
          table1$GROUP=paste0(table1$SEX,table1$TYPE)
          table1$SPECIES="PPUN"
          ##G.acu
          use.chromosomes=c(use.chromosomes,paste0("chr",as.roman(c(8,12))))
          table2$GROUP=paste0(table2$SEX,table2$TYPE)
          table2$SPECIES="GACU"
          ##
          plot.table=rbind(table2,table1)
          plot.table$SPECIES=factor(plot.table$SPECIES,labels=c("P. pungitius", "G. aculeatus"),levels=c("PPUN","GACU"))
          plot.table$SEX=factor(plot.table$SEX,labels=c("Maternal", "Paternal"),levels=c("MATERNAL","PATERNAL"))
          plot.table$TYPE=factor(plot.table$TYPE,labels=c("Empirical", "Predicted"),levels=c("EMPIRICAL","PREDICTED"))
          plot.table$CHR=factor(plot.table$CHR,levels=c("chrVIII","chrXII","LG7","LG10"))
          plot.table=subset(plot.table,CHR %in% use.chromosomes)
          
          text.table=plot.table[seq(1,nrow(plot.table),by=15),]
          text.table$EXPLABEL=paste0("E=",signif(as.numeric(text.table$EXPECTED),4))
          text.table$VARLABEL=paste0("var=",signif(as.numeric(text.table$VARIANCE),5))
        
          p=ggplot(plot.table,aes(x=BIN,y=DENSITY,group=GROUP,color=SEX,fill=SEX)) + geom_line(aes(linetype=TYPE)) +  
          geom_ribbon(aes(x=BIN,ymax=DENSITY),linetype=0,ymin=0,alpha=0.3) +
          geom_text(data=subset(text.table,SEX=="Maternal"),size=3.0, mapping = aes(x = 8, y = 0.400, label = EXPLABEL),hjust="middle", vjust="middle",show.legend=F) +
          geom_text(data=subset(text.table,SEX=="Maternal"),size=3.0, mapping = aes(x = 8, y = 0.365, label = VARLABEL),hjust="middle", vjust="middle",show.legend=F) +
          geom_text(data=subset(text.table,SEX=="Paternal"),size=3.0, mapping = aes(x = 8, y = 0.330, label = EXPLABEL),hjust="middle", vjust="middle",show.legend=F) +
          geom_text(data=subset(text.table,SEX=="Paternal"),size=3.0, mapping = aes(x = 8, y = 0.295, label = VARLABEL),hjust="middle", vjust="middle",show.legend=F) +
          facet_grid(COMPARISON ~ SPECIES + CHR) + 
          theme(strip.text.y=element_blank(),panel.grid=element_blank(),legend.title=element_blank(),legend.position="bottom") + 
          scale_x_continuous(breaks=c(1,15)) +
          scale_color_manual(values=c("tan2","steelblue2")) + scale_fill_manual(values=c("tan2","steelblue2")) +
          ylab("Density") + xlab ("Chromosome bin")
          ggsave("text.png",width=118,height=150,units="mm",dpi=1000)
          ggsave("text.pdf",width=118,height=150,units="mm")

     """
     
     )
     
     
     val plotBetweenIndividualVariance = REvaluate(
         table1=concatenateCrossovers3spn,
         table2=chromosomeMetadata3spn,
         table3=concatenateCrossovers9spn,
         table4=chromosomeMetadata9spn,
         script="""
         library(ggplot2)
         table.out=data.frame()
         rm('optOut1')
		  fig.dir <- get.output(cf, 'optOut1')
		  dir.create(fig.dir, recursive=TRUE)
		  setwd(fig.dir)
		  
		  ITER=100
		######
			h=9
			pascalTriangle <- function(h) {
			   lapply(0:h, function(i) choose(i, 0:i))
			}
			P=pascalTriangle(h)
		sample.co.count=function(CHR,PROB.TABLE){  
			   probs=as.numeric(PROB.TABLE[PROB.TABLE$CHR==CHR,2:11])
			   probs.variant=probs
			   probs.variant[3]=probs.variant[3]+sum(probs[4:10])
			   probs.variant[4:10]=0
			   co.count=sample(x=seq(0,9),size=1,prob=probs.variant)
			   observed.co.count=sample(x=seq(0,co.count),size=1,prob=P[[co.count+1]]*0.5^co.count)
			   return(observed.co.count)
			}
		#####
		##G. aculeatus
		PROB.TABLE.MALE=table2[,c(1,32:41)]
		PROB.TABLE.FEMALE=table2[,c(1,51:60)]
		
		d=table1 
		d1=subset(d,CHR!="chrXIX")
		PROB.TABLE.MALE=table2[,c("CHR",paste0("PATERNALn",0:5),paste0("PATERNALMLEp",0:9))]
		PROB.TABLE.FEMALE=table2[,c("CHR",paste0("MATERNALn",0:5),paste0("MATERNALMLEp",0:9))]
		male.means=aggregate(x=d1$PATERNALCOUNT,by=list(d1$MALE),mean)
		female.means=aggregate(x=d1$MATERNALCOUNT,by=list(d1$FEMALE),mean)
		n.offspring.male=as.numeric(table(subset(d1,CHR=="chrI")$MALE))
		n.offspring.female=as.numeric(table(subset(d1,CHR=="chrI")$FEMALE))
		n.parents.male=length(n.offspring.male)
		n.parents.female=length(n.offspring.female)
		variances.male=c()
		for(j in 1:ITER){
			out.male=matrix(ncol=n.parents.male)
			for(CHROM in paste0("chr",as.roman(c(1:18,20:21)))){
			    means=c()
			    i=1
			    while(i<=n.parents.male){
			        means=c(means,mean(replicate(n.offspring.male[i],sample.co.count(CHROM,PROB.TABLE.MALE))))
			        i=i+1
			    }
			    out.male=rbind(out.male,means)
			}
			out.male=out.male[-1,]
		    variances.male=c(variances.male,var(colMeans(out.male)))
		}
		var.test(colMeans(out.male),male.means[male.means$x>0.2,2])
		GRAND.MEAN=mean(colMeans(out.male))
		VAR=var(colMeans(out.male))
		####################
		variances.female=c()
		for(j in 1:ITER){
			out.female=matrix(ncol=n.parents.female)
			for(CHROM in paste0("chr",as.roman(c(1:18,20:21)))){
			    means=c()
			    i=1
			    while(i<=n.parents.female){
			        means=c(means,mean(replicate(n.offspring.female[i],sample.co.count(CHROM,PROB.TABLE.FEMALE))))
			        i=i+1
			    }
			    out.female=rbind(out.female,means)
			}
			out.female=out.female[-1,]
		    variances.female=c(variances.female,var(colMeans(out.female)))
		}
		GRAND.MEAN=mean(colMeans(out.female))
		VAR=var(colMeans(out.female))
		
		total=data.frame(SPECIES="G. aculeatus",SEX=rep(c("Paternal","Maternal"),times=c(length(male.means$x)+ITER,length(female.means$x)+ITER)),VALUE=c(male.means$x,variances.male,female.means$x,variances.female),GROUP=rep(c("Empirical","Simulated","Empirical","Simulated"),times=c(length(male.means$x),ITER,length(female.means$x),ITER)))
		variance.points=data.frame(SPECIES="G. aculeatus",SEX=c("Paternal","Maternal"),VALUE=c(var(male.means$x),var(female.means$x)),GROUP="Simulated")
		####
		# P. pungitius
		d=table3 
		d1=subset(d,CHR!="LG12" & !(MALE %in% c("88-m-2","72-m-2","69-m-1","66-m-2")))
		PROB.TABLE.MALE=table4[,c("CHR",paste0("PATERNALn",0:5),paste0("PATERNALMLEp",0:9))] 
		PROB.TABLE.FEMALE=table4[,c("CHR",paste0("MATERNALn",0:5),paste0("MATERNALMLEp",0:9))]
		male.means=aggregate(x=d1$PATERNALCOUNT,by=list(d1$MALE),mean)
		female.means=aggregate(x=d1$MATERNALCOUNT,by=list(d1$FEMALE),mean)
		n.offspring.male=as.numeric(table(subset(d1,CHR=="LG1")$MALE))
		n.offspring.female=as.numeric(table(subset(d1,CHR=="LG1")$FEMALE))
		n.parents.male=length(n.offspring.male)
		n.parents.female=length(n.offspring.female)
		variances.male=c()
		for(j in 1:ITER){
			out.male=matrix(ncol=n.parents.male)
			for(CHR in paste0("LG",c(seq(1,11),seq(13,21)))){
			    means=c()
			    i=1
			    while(i<=n.parents.male){
			        means=c(means,mean(replicate(n.offspring.male[i],sample.co.count(CHR,PROB.TABLE.MALE))))
			        i=i+1
			    }
			    out.male=rbind(out.male,means)
			}
			out.male=out.male[-1,]
		    variances.male=c(variances.male,var(colMeans(out.male)))
		}
		var.test(colMeans(out.male),male.means[male.means$x>0.2,2])
		GRAND.MEAN=mean(colMeans(out.male))
		VAR=var(colMeans(out.male))
		#####################
		variances.female=c()
		for(j in 1:ITER){
			out.female=matrix(ncol=n.parents.female)
			for(CHR in paste0("LG",c(seq(1,11),seq(13,21)))){
			    means=c()
			    i=1
			    while(i<=n.parents.female){
			        means=c(means,mean(replicate(n.offspring.female[i],sample.co.count(CHR,PROB.TABLE.FEMALE))))
			        i=i+1
			    }
			    out.female=rbind(out.female,means)
			}
			out.female=out.female[-1,]
		    variances.female=c(variances.female,var(colMeans(out.female)))
		}
		GRAND.MEAN=mean(colMeans(out.female))
		VAR=var(colMeans(out.female))
		
        ppun.total=data.frame(SPECIES="P. pungitius",SEX=rep(c("Paternal","Maternal"),times=c(length(male.means$x)+ITER,length(female.means$x)+ITER)),VALUE=c(male.means$x,variances.male,female.means$x,variances.female),GROUP=rep(c("Empirical","Simulated","Empirical","Simulated"),times=c(length(male.means$x),ITER,length(female.means$x),ITER)))
        ppun.variance.points=data.frame(SPECIES="P. pungitius",SEX=c("Paternal","Maternal"),VALUE=c(var(male.means$x),var(female.means$x)),GROUP="Simulated")
		both.total=rbind(ppun.total,total)
		variance.both=rbind(ppun.variance.points,variance.points)
		
		both.total$COL=paste(both.total$SPECIES,both.total$SEX,sep=" ")
		variance.both$COL=paste(variance.both$SPECIES,variance.both$SEX,sep=" ")
		print(head(both.total))
		ggplot(both.total,aes(VALUE,fill=COL,col=COL)) + geom_histogram() + geom_point(data=variance.both,aes(x=VALUE,y=0),shape=21,color="black") +
		scale_color_manual(values=c("tan3","steelblue3","tan1","steelblue1")) + scale_fill_manual(values=c("tan3","steelblue3","tan1","steelblue1")) +
		facet_wrap(GROUP~SPECIES,scales="free") +
		theme(panel.grid=element_blank(),legend.position="none")
        ggsave("both.pdf",width=110,height=110,units="mm")        
		"""		  
      
     
     
     )

}
// d$CHR=factor(d$CHR,levels=paste0("LG",1:21))
//p=ggplot(d,aes(x=BIN,y=DENSITY,fill=as.factor(CO),col=as.factor(CO)))+geom_line()+geom_ribbon(aes(x=BIN,ymax=DENSITY),ymin=0,alpha=0.3) + facet_grid(CHR ~ SEX + CO.GROUP,scales="free_y")
