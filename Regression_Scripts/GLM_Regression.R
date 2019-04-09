#Factorial experimental design: General linear mixed model
#Understanding the effects of rhizobia genotype on flowering time.

setwd("/Users/wendy/Dropbox/FloweringTimeConflict")
flo<-read.csv("FLOWERING_MASTER_FILE_JAN2014.csv")

no.nit<-flo[flo[,3]!="-N",]
no.nit<-no.nit[no.nit[,3]!="'+N",]

no.nit1<-no.nit[no.nit[,50]!=1,] # plants that died early in life

###FLOWERING TIME ANALYSIS:

library(car)

Anova(summary(lm(days_flower~soil.origin.geno*soil.origin.rhiz*as.factor(treat),family=poisson(link="log"),data=no.nit)->lm1)

library(lme4)
str(no.nit)

flo<-glmer(days_flower~soil.origin.geno*soil.origin.rhiz*as.factor(treat) + (1|rhiz) + (1|geno)+(1|id),family=poisson(link="log"), data=no.nit)
flo1<-glmer(days_flower~soil.origin.geno*soil.origin.rhiz*as.factor(treat) + (1|rhiz)+(1|id),family=poisson(link="log"), data=no.nit)
flo2<-glmer(days_flower~soil.origin.geno*soil.origin.rhiz*as.factor(treat) + (1|geno)+(1|id),family=poisson(link="log"), data=no.nit)

#determine significance of rhiz and geno random effects:

anova(flo1,flo) # Significant (1|geno) p-value = 2.2e-16
anova(flo2,flo) # Not Significant (1|rhiz) p-value = 0.3484


#Checking for overdispersion: http://glmm.wikidot.com/faq  This works for the lme4 and glmmADMB package

overdisp_fun <- function(model) {
  ## number of variance parameters in 
  ##   an n-by-n variance-covariance matrix
  vpars <- function(m) {
    nrow(m)*(nrow(m)+1)/2
  }
  model.df <- sum(sapply(VarCorr(model),vpars))+length(fixef(model))
  rdf <- nrow(model.frame(model))-model.df
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}


#Are there nodule effects on flowering time? YES!!!

Anova(lm(sqrt(days_flower)~soil.origin.geno*soil.origin.rhiz*as.factor(treat) + nod.T,data=no.nit)->lm1)

nod<-glmer(days_flower~soil.origin.geno*soil.origin.rhiz*as.factor(treat) +nod.T+ (1|rhiz) + (1|geno)+(1|id),family=poisson(link="log"), data=no.nit)

	# 1|geno:rhiz ==> geno:rhizobia interaction as a random effect in the model
	 
nod1<-glmer(days_flower~soil.origin.geno*soil.origin.rhiz*as.factor(treat) +nod.T + (1|geno:rhiz)+(1|id),family=poisson(link="log"), data=no.nit)

nod2<-glmer(days_flower~soil.origin.geno*soil.origin.rhiz*as.factor(treat) +nod.T +(1|id),family=poisson(link="log"), data=no.nit)

anova(nod2,nod1) # lr test: pval= 2.2e-16

cor.test(no.nit$days_flower,no.nit$nod.T,na.rm=TRUE)
plot(no.nit$days_flower,no.nit$nod.T)
abline(lm(no.nit$nod.T~no.nit$days_flower))

#Are there nodule effects on flowering time with respect to treatment?
par(mfcol=c(1,2))
cor.test(nod0$nod.T,nod0$days_flower,na.rm=TRUE)# NO CORR for 0mM salt: p-value= 0.4505
plot(nod0$nod.T,nod0$days_flower, ylab="days to flowering",xlab="nodule number", main="0mM Salt")
abline(lm(nod0$days_flower~nod0$nod.T))

cor.test(nod100$nod.T,nod100$days_flower,na.rm=TRUE) # There is a significant CORR: p-value=8.82e-05
plot(nod100$nod.T,nod100$days_flower,ylab="", xlab="nodule number", main="100mM Salt")
abline(lm(nod100$days_flower~nod100$nod.T))

#NODULE ANALYSIS:
#Are there differences in total nodule number between saline and nonsaline conditions? YES!!!!

nod0<-no.nit[no.nit[,6]==0,]
nod100<-no.nit[no.nit[,6]==100,]

wilcox.test(nod0[,42],nod100[,42],alternative = "greater") #p-value=0.0002224, mean.nod0=28.83, mean.nod100=23.62
	
	#BARPLOT of nodule number difference between treatment
nods<-cbind(nod0[,42],nod100[,42])
nod.mean<-apply(nods,2,mean,na.rm=TRUE)
nod.sd<-apply(nods,2,sd,na.rm=TRUE)
barx<-barplot(nod.mean,ylim=c(0,50),col=c("blue","red"),axis.lty=1,xlab="Treatment",ylab="nodule number")
error.bar(barx,nod.mean,1.96*nod.sd/10)

error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
stop("vectors must be same length")
arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

#Are there differences in total nodule number between matched and mismatched pairs? NOTHING interesting

nodNS.NS<-no.nit[no.nit[,4]=="NS" & no.nit[,5]=="NS",]
nodNS.S<-no.nit[no.nit[,4]=="NS" & no.nit[,5]=="S",]
nodS.S<-no.nit[no.nit[,4]=="S" & no.nit[,5]=="S",]
nodS.NS<-no.nit[no.nit[,4]=="S" & no.nit[,5]=="NS",]

wilcox.test(nodS.S[,42],nodS.NS[,42])

nodNS.S0<-nodNS.S[nodNS.S[,6]=="0",]
nodNS.NS0<-nodNS.NS[nodNS.NS[,6]=="0",]

wilcox.test(nodNS.NS0[,42],nodNS.S0[,42]) # no difference in nodule numbers pval=0.7707 mean.nodNS.S0

mean(nodNS.S0[,42],na.rm=TRUE)


###Are there flowering effects on Nodule numbers? YES!!!
Anova(lm(sqrt(nod.T)~soil.origin.geno*soil.origin.rhiz*as.factor(treat) + days_flower,data=no.nit)->lm1)
hist(lm1$residuals)


nod1<-glmer(nod.T~soil.origin.geno*soil.origin.rhiz*as.factor(treat) +days_flower + (1|geno:rhiz)+(1|id),family=poisson(link="log"), data=no.nit)

###Are there soil.origin effects with respect to treatment? YES!! 0mM Salt: sig. Soil.Orig.Geno*Soil.Orig.Rhiz*flowering : WRONG!!!! 

Anova(lm(sqrt(nod.T)~soil.origin.geno*soil.origin.rhiz*days_flower,data=nod0)->lm1)
Anova(lm(sqrt(nod.T)~soil.origin.geno*soil.origin.rhiz*days_flower,data=nod100)->lm2)

##Are there differences in pod/seed production between treatments?

pod0<-nod0[,33]
pod100<-nod100[,33]
mean(pod0,na.rm=TRUE)
mean(pod100,na.rm=TRUE)

seed0<-nod0[,31]
seed100<-nod100[,31]
mean(seed0,na.rm=TRUE)
mean(seed100,na.rm=TRUE)

wilcox.test(pod0,pod100)
wilcox.test(seed0,seed100)

	#BARPLOT of seed number difference between treatment
	
nods<-cbind(nod0[,31],nod100[,31])
nod.mean<-apply(nods,2,mean,na.rm=TRUE)
nod.sd<-apply(nods,2,sd,na.rm=TRUE)
barx<-barplot(nod.mean,ylim=c(0,20),col=c("blue","red"),axis.lty=1,xlab="Treatment",ylab="seed number")
error.bar(barx,nod.mean,1.96*nod.sd/10)

error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
stop("vectors must be same length")
arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

##Are there differences in leaf production between treatments?

leaf0<-nod0[,11]
leaf100<-nod100[,11]
mean(leaf0,na.rm=TRUE)
mean(leaf100,na.rm=TRUE)
wilcox.test(leaf0,leaf100)


	#BARPLOT of leaf number difference between treatment
	
nods<-cbind(nod0[,11],nod100[,11])
nod.mean<-apply(nods,2,mean,na.rm=TRUE)
nod.sd<-apply(nods,2,sd,na.rm=TRUE)
barx<-barplot(nod.mean,ylim=c(0,20),col=c("blue","red"),axis.lty=1,xlab="Treatment",ylab="leaf number")
error.bar(barx,nod.mean,1.96*nod.sd/10)

error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
stop("vectors must be same length")
arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}
