#The following code creates a QQplot of pvalues from logit regression analysis

ggd.qqplot = function(pvector, main=NULL, ...) {
  o = -log10(sort(pvector,decreasing=F))
  e = -log10( 1:length(o)/length(o) )
  plot(e,o,pch=19,cex=1, main=main, ...,
       xlab=expression(Expected~~-log[10](italic(p))),
       ylab=expression(Observed~~-log[10](italic(p))),
       xlim=c(0,max(e)), ylim=c(0,max(o)))
  abline(0,1,col="red")
}

set.seed(42)
pvalues=runif(10000)
pvalues[sample(10000,10)]=pvalues[sample(10000,10)]/5000

ggd.qqplot(pvalues)

##This QQplot is done on the all snps that have at least 55% genotype calls
pval <- read.table("/Users/wendy/Dropbox/Sol96_LogitRegression/LogitRegression_Files/FDR_NONFiltered_Logit_pval.txt",header=TRUE,sep=",")

pval.lat <-pval[which(pval[,1] != 0.000000e+00),1]
pval.lon <-pval[which(pval[,2] != 0),2]
pval.ele <-pval[which(pval[,3] != 0),3]
pval.ec <-pval[which(pval[,4] != 0),4]


par(mfrow=c(2,2))

hist(pval.lat,main="Latitude",xlab="Pvalue (FDR)")
hist(pval.lon,main="Longitude",xlab="Pvalue (FDR)")
hist(pval.ele,main="Elevation",xlab="Pvalue (FDR)")
hist(pval.ec,main="Salinity",xlab="Pvalue (FDR)")

par(mfrow=c(1,1))
ggd.qqplot(pval.lat, "Latitude")
ggd.qqplot(pval.lon, "Longitude")
ggd.qqplot(pval.ele, "Elevation")
ggd.qqplot(pval.ec, "Salinity")

snp.60<-read.table("/Users/wendy/Dropbox/Sol96_LogitRegression/LogitRegression_Files/Sol96_Freq_min2reads.60.Major_Minor_Allele.txt",header=T,sep=",")
snp.70<-read.table("/Users/wendy/Dropbox/Sol96_LogitRegression/LogitRegression_Files/Sol96_Freq_min2reads.70.Major_Minor_Allele.txt",header=T,sep=",")
snp.80<-read.table("/Users/wendy/Dropbox/Sol96_LogitRegression/LogitRegression_Files/Sol96_Freq_min2reads.80.Major_Minor_Allele.txt",header=T,sep=",")

snp60<-paste(snp.60$CHROM,snp.60$POS,snp.60$REF,sep="_")
snp70<-paste(snp.70$CHROM,snp.70$POS,snp.70$REF,sep="_")
snp80<-paste(snp.80$CHROM,snp.80$POS,snp.80$REF,sep="_")

pval60<-pval[pval[,1] %in% snp60,]
lat.60<-pval60[which(pval60[,2] != 0),1]
lon.60<-pval60[which(pval60[,3] != 0),2]
ele.60<-pval60[which(pval60[,4] != 0),3]
ec.60<-pval60[which(pval60[,5] != 0),4]

ggd.qqplot(lat.60, "Latitude")
ggd.qqplot(lon.60, "Longitude")
ggd.qqplot(ele.60, "Elevation")
ggd.qqplot(ec.60, "Salinity")

pval70<-pval[pval[,1] %in% snp70,]
lat.70<-pval70[which(pval70[,2] != 0),2]
lon.70<-pval70[which(pval70[,3] != 0),3]
ele.70<-pval70[which(pval70[,4] != 0),4]
ec.70<-pval70[which(pval70[,5] != 0),5]

ggd.qqplot(lat.70, "Latitude")
ggd.qqplot(lon.70, "Longitude")
ggd.qqplot(ele.70, "Elevation")
ggd.qqplot(ec.70, "Salinity")

pval80<-pval[pval[,1] %in% snp80,]
lat.80<-pval80[which(pval80[,2] != 0),2]
lon.80<-pval80[which(pval80[,3] != 0),3]
ele.80<-pval80[which(pval80[,4] != 0),4]
ec.80<-pval80[which(pval80[,5] != 0),5]

ggd.qqplot(lat.80, "Latitude")
ggd.qqplot(lon.80, "Longitude")
ggd.qqplot(ele.80, "Elevation")
ggd.qqplot(ec.80, "Salinity")

