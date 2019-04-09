#filter out snp positions that have less than 4 reads to eliminate the structure I found in the data
#######
#snp<-read.table("/Users/wendy/Dropbox/Sol96_LogitRegression/Minor_Major_SNP_data.txt",header=T,sep=",")
geno.reads<-read.table("/Users/wendy/Dropbox/Sol96_LogitRegression/vcf_geno_depth_filter_200KSNP.txt",header=T,sep="\t")
reads.mat<-as.matrix(geno.reads[,9:99])
#row.col.idx<-which(reads.mat < 4, arr.ind =TRUE)
pos.con <-paste(geno.reads$CHROM,geno.reads$POS,geno.reads$REF, sep="_")
row.names(reads.mat)<-pos.con
#######
#use file below for PCA: Genoytpes coded 0,0.5,1
#snp<-read.table("/Users/wendy/Dropbox/Sol96_LogitRegression/LogitRegression_Files/Sol96_Freq_min2reads.55.Major_Minor_Allele.txt",header=T,sep=",")

#use file below for Ordinal Log Regression: genotypes coded by genotype freq --> 1,2,3
snp <-read.table("/Users/wendy/Dropbox/Sol96_LogitRegression/LogitRegression_Files/LogisticReg_GenoFrq_Recode_Cat.txt",header=T,sep=",")

#file key for Log Reg analysis: identifies the genotypes (0,0.5,1) from genotype freq file above
snp.key <- read.table("/Users/wendy/Dropbox/Sol96_LogitRegression/LogitRegression_Files/LogisticReg_GenoFrq_Recode_Cat_Key.txt",header=T,sep=",")
snp.60<-read.table("/Users/wendy/Dropbox/Sol96_LogitRegression/LogitRegression_Files/Sol96_Freq_min2reads.60.Major_Minor_Allele.txt",header=T,sep=",")

het.55<-read.table("/Users/wendy/Dropbox/Sol96_LogitRegression/LogitRegression_Files/Sol96_Freq_min2reads.55.het",header=T, sep="\t" )
het.60<-read.table("/Users/wendy/Dropbox/Sol96_LogitRegression/LogitRegression_Files/Sol96_Freq_min2reads.60.het",header=T, sep="\t" )

#ref vcftools for genotype filter: https://github.com/jpuritz/dDocent/blob/master/tutorials/Filtering%20Tutorial.md


###Logit Regression#######################
#Preparing snp.pf dataframe for Logit regression

snp.frq<-t(snp.pos)
class(snp.frq) <- "numeric"

#preparing phenotype/predictor data for logit regression

base<-read.csv("/Users/wendy/Dropbox/Tails_Data_Soliman250/GPS_96tails/GPS_BaseCamp_Soliman_Points.csv")
sol96<-read.csv("/Users/wendy/Dropbox/Tails_Data_Soliman250/Sol96_Pheno_Analysis/Pheno_Sol_Seq96_Orig.csv")

geno.names<-colnames(snp[,4:length(colnames(snp))])

avgec96<-aggregate(sol96[,6],list(line=sol96$line),mean)
base96<-base[base[,1] %in% avgec96[,1],]


mean.seed<-aggregate(sol96[,9],list(line=sol96$line),mean,na.rm=T)
mean.leaf<-aggregate(sol96[,27],list(line=sol96$line),mean, na.rm=T)
mean.flo96<-aggregate(sol96[,34],list(line=sol96$line),mean,na.rm=T)
mean.pod<-aggregate(sol96[,37],list(line=sol96$line),mean, na.rm=T)

pheno<-cbind(base96[,c(1,3:5)],avgec96[,2],mean.seed[,2],mean.leaf[,2],mean.flo96[,2],mean.pod[,2])
pheno<-pheno[,-1]
colnames(pheno)<-c("Lat","Long","Ele","EC","SeedMass","LeafNum","FlowerTime","PodNum")
rownames(pheno)<-geno.names

#removing rows/genotypes with NA values from pheno and then corresponding snp.frq data

row.has.na <- apply(pheno, 1, function(x){any(is.na(x))})
remove.geno<-names(which(row.has.na=="TRUE"))

pheno1<-pheno[!(rownames(pheno) %in% remove.geno),]
snp.frq1<-snp.frq[!(rownames(snp.frq) %in% remove.geno),]


#STREAMLINING the OLR to analyze ~100K snps independently
#NOTE: This model will only work for response variable with 3 or more levels or categories
library(foreign)
library(ggplot2)
library(MASS)
library(Hmisc)
library(reshape2)

#####Components for function
my.snp<-snp.frq[,3]
polr(as.factor(my.snp) ~ lon + ele + ec, Hess=TRUE)->m
(t.val<-coef(summary(m))[,3])
pnorm(abs(t.val),lower.tail=FALSE)*2
#intercept == cutpoints: 

#####
#ERROR:OLR will only analyze 3 or more response levels/categories. Some of the snps have < 3 levels, so we must use tryCatch method to catch the error and continue to the next snp 
#http://stackoverflow.com/questions/8093914/skip-to-next-value-of-loop-upon-error-in-r-trycatch

#####FINAL SCRIPT to run Ordered logistic regression 
snp.frq1<-snp.frq[,1:4643]
snp.frq2<-snp.frq[,4645:12766]
snp.frq3<-snp.frq[,12768:16868]
snp.frq4<-snp.frq[,16870:length(snp.frq[1,])]
lat<-pheno[,1]
lon<-pheno[,2]
ele<-pheno[,3]
ec<-pheno[,4]

###Pvalue matrix for the additive model 
i=1
do.olr<-function(i){
  print(i)
  my.snp<-snp.frq[,i]
  error<-tryCatch(polr(as.factor(my.snp) ~ lon + log(ele) + log(ec), Hess=TRUE,method="logistic",na.action=na.omit),  error=function(e) e)
  # e= is any error messages that originate in the expression that you are tryCatch - ing.
  #print(error) This tells you the class of the error
  #class(error) This returns an object that can be used for a conditional statement below
  if(!inherits(error,"simpleError")){ #if the snp doesn't generate an error, run ols model and calculate pval
    model<-polr(as.factor(my.snp) ~ lat + lon + log(ele) + log(ec), Hess=TRUE, method="logistic",na.action=na.omit)
    t.val<-coef(summary(model))[,3] #extracting t.val from coefficient table
    pnorm(abs(t.val),lower.tail=FALSE)*2
  } else{c("NA","NA","NA","NA","NA","NA")} #returns NA's if there is an error (there are six variables in the output)
  
}

sapply(1:length(snp.frq[1,]),do.olr)->ols.pval


#ols.pval: columns correspond to snp position 
#rows correspond to lat, lon, ele, ec, intercept 0|0.5, intercept 0.5|1

gene_id<-colnames(snp.frq)
variable_id<-rownames(ols.pval)
colnames(ols.pval)<-gene_id #add gene id to column names
pval<-t(ols.pval) #transpose matrix column-->row

#Filter out the snps with less 3 levels before fdr correction: remove snps with NA's
class(pval)<-"numeric"
idx<-apply(pval, 1, function(x) all(is.na(x)))
pval<-pval[!idx,] #After filtering there are 5225 snps that were interrogated

##FDR correction for multiple testing
i=1
correct.multi.test<-function(i){
  p.adjust(as.numeric(pval[,i]),method="fdr")
}
sapply(1:length(pval[1,]),correct.multi.test)->fdr

gene_id<-rownames(pval)
variable_id<-colnames(pval)

colnames(fdr)<-variable_id
rownames(fdr)<-gene_id

fdr1<-fdr
fdr2<-fdr
fdr3<-fdr
fdr4<-fdr

fdr<-rbind(fdr1,fdr2,fdr3,fdr4)
fdr.fin <-as.data.frame(cbind(rownames(fdr),fdr))
colnames(fdr.fin)[1]<-"SNP_Pos_Ref"
write.table(fdr.fin,"/Users/wendy/Dropbox/Sol96_LogitRegression/LogitRegression_Files/FDR_NONFiltered_Logit_pval.txt",
            col.names=TRUE,row.names = FALSE,quote=FALSE,sep=",")
#FDR < 0.05

sig.lat<-fdr[which(fdr[,1] < 0.05),] 
sig.lon<-fdr[which(fdr[,2] < 0.05),] 
sig.ele<-fdr[which(fdr[,3] < 0.05),] 
sig.ec<-fdr[which(fdr[,4] < 0.05),] 

###SNPs filtered at 55% criteria
#filter out pval = 0 from data
  # snp.frq[,1:4643]
sig.lat1<-sig.lat[which(sig.lat[,1]!=0),1]
sig.lon1<-sig.lon[which(sig.lon[,2]!=0),2]
sig.ele1<-sig.ele[which(sig.ele[,3]!=0),3]
sig.ec1<-sig.ec[which(sig.ec[,4]!=0),4]

  #snp.frq[,4645:12766]
sig.lat2<-sig.lat[which(sig.lat[,1]!=0),1]
sig.lon2<-sig.lon[which(sig.lon[,2]!=0),2]
sig.ele2<-sig.ele[which(sig.ele[,3]!=0),3]
sig.ec2<-sig.ec[which(sig.ec[,4]!=0),4]

  #snp.frq[,12768:16868]
sig.lat3<-sig.lat[which(sig.lat[,1]!=0),1]
sig.lon3<-sig.lon[which(sig.lon[,2]!=0),2]
sig.ele3<-sig.ele[which(sig.ele[,3]!=0),3]
sig.ec3<-sig.ec[which(sig.ec[,4]!=0),4]

  #snp.frq[,16870:length(snp.frq[1,])]
sig.lat4<-sig.lat[which(sig.lat[,1]!=0),1]
sig.lon4<-sig.lon[which(sig.lon[,2]!=0),2]
sig.ele4<-sig.ele[which(sig.ele[,3]!=0),3]
sig.ec4<-sig.ec[which(sig.ec[,4]!=0),4]

##Combine 
#SNPs filtered at 55% criteria
sig.lat.fdr <- c(sig.lat1,sig.lat2,sig.lat3,sig.lat4) #6192 snps
sig.lon.fdr <- c(sig.lon1,sig.lon2,sig.lon3,sig.lon4) #300 snps
sig.ele.fdr <- c(sig.ele1,sig.ele2,sig.ele3,sig.ele4) #191 snps
sig.ec.fdr <- c(sig.ec1,sig.ec2,sig.ec3,sig.ec4) #60 snps

sig.lat.fdr.55 <-cbind(names(sig.lat.fdr),as.data.frame(sig.lat.fdr),rep("lat",times=length(sig.lat.fdr)))
sig.lon.fdr.55 <-cbind(names(sig.lon.fdr),as.data.frame(sig.lon.fdr),rep("lon",times=length(sig.lon.fdr)))
sig.ele.fdr.55 <-cbind(names(sig.ele.fdr),as.data.frame(sig.ele.fdr),rep("ele",times=length(sig.ele.fdr)))
sig.ec.fdr.55 <-cbind(names(sig.ec.fdr),as.data.frame(sig.ec.fdr),rep("ec",times=length(sig.ec.fdr)))

colnames(sig.ec.fdr.55) <- c("SNP_Pos","Pval_FDR","Factor")

sig.all.fdr.55 <-rbind(sig.lat.fdr.55, sig.lon.fdr.55, sig.ele.fdr.55, sig.ec.fdr.55)

write.table(sig.all.fdr.55,"/Users/wendy/Dropbox/Sol96_LogitRegression/LogitRegression_Files/Sig.Lat.Lon.Ele.EC.fdr.55.txt",quote=FALSE, row.names = FALSE,col.names = TRUE,sep="\t")


#What is the relationship between STRUCTURE and heterozygosity estimates?
plot(str.sort.geno[,6],het[,5],col=ifelse(str.sort.geno[,8]=="hap1","black","red"),pch=19,xlab="Proportion of Haplotype 1 (Dark blue)", ylab="Inbreeding coefficient (F)")

#What is the relationship between heterozygosity estimates and genetic diversity
idx<-c(1:96)
het<-cbind(idx,het)
het.sort<-het[order(het[,6],decreasing=FALSE),]
het.num<-c(rep("het1",times=length(which(het.sort[,6]<0))),rep("het2",times=length(which(het.sort[,6]>0))))
het1<-cbind(het.sort,het.num)
het.sort.geno<-het1[order(het1[,1]),]

plot(pc1, pc2, col=ifelse(het.sort.geno[,7]=="het1","black","red"),xlab="PC1 (71.7%)",ylab="PC2 (1.2%)",pch=19,main="PC1-PC2 by Het")





