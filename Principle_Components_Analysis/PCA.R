####################################################
#Principle components analysis for American oil palm samples. Input file (tab.frq) is a modified SNP file from genomic sequencing data: SNPs were recoded to correspond to major and minor allele frequencies (1/0.5/0). I used a custom python script to recode the SNPs.
####################################################

#####PRINCIPLE COMPONENTS ANALYSIS
tab.frq <- read.table("/Users/wvu/Dropbox/Oil_Palm/RAD_SequencingResults_Dec2017/American_diversity/PCA_Het/SOL_OILPALM.Eguineensis9.1.VF.SNP.DP8.95.Major_Minor_Allele_hcpc.tp.bra.CultWildG.hybrids.txt",header=T,sep=",")

####ALL SAMPLES 
#reformat data.frame columns=sampleID, rows=snps
pos.concat<-paste(tab.frq$CHROM, tab.frq$POS, tab.frq$REF,sep="_")
snp.pos<-cbind(pos.concat,tab.frq[,4:length(tab.frq[1,])])
row.names(snp.pos)<-snp.pos[,1]
snp.pos<-snp.pos[,-1]
dim(snp.pos)

#extract american samples
snp.pos <- snp.pos[,colnames(snp.pos) %in% as.character(amer.g.hyb[,1])]

#reorder columns
snp.pos = snp.pos[,order(match(colnames(snp.pos),het[,1]))]

#generate a covariance matrix
covmat<-cov(snp.pos,use="pairwise")

#estimate eigenvectors
ecov<-eigen(covmat)

#rename eigenvectors by country and cross
rownames(ecov$vectors)<-het$country
rownames(ecov$vectors)<-het$pop

pc1<-ecov$vectors[,1]
pc2<-ecov$vectors[,2]

#calculate Variation explained by PC1 and PC2
sum<-sum(ecov$values)
ecov$values[1]/sum #0.51
ecov$values[2]/sum #0.10

##############MAKING PCA PLOTS WITH SAMPLE LABELS (HCPC_cultG + wild hybrids)###################################
col = c(rep("lightseagreen",length(het[het[,2]=='HND',1])), rep("purple",length(het[het[,2]=='CRI',1])), 
        rep("gold",length(het[het[,2]=='PAN',1])), rep("deeppink",length(het[het[,2]=='COL',1])), 
        rep("black",length(het[het[,2]=='TAI',1])), rep("red",length(het[het[,2]=='PER',1])), 
        rep("darkgreen",length(het[het[,2]=='BRA',1])),rep("gold4",length(het[het[,2]=='Cultivated_G',1])))

col.legend = c("lightseagreen","purple", 
        "gold","deeppink", 
        "black","red", 
        "darkgreen","gold4")

plot(pc1,pc2,xlab="PC1 (0.40)",ylab="PC2 (0.23)",col=col,pch=19,cex=2,ylim=c(-0.30,0.05),xlim=c(-0.25,0.2))
legend(-0.265,-0.16,legend=unique(names(pc1)),col=col.legend,pch=19,cex=2)
text(pc1,pc2,labels=names(pc1),col=col)


#####################BOXPLOTS OF HETERZYGOSITY ESTIMATES
het <- read.csv("/Users/wvu/Dropbox/Oil_Palm/RAD_SequencingResults_Dec2017/American_diversity/PCA_Het/SOL_OILPALM.Eguineensis9.1.VF.SNP.DP8.90.vcf.hcpc.tp.bra.CultWildG.hybrids.het.csv",header=T)

#########Sort by country
sort1 <- c("HND","CRI","PAN","COL","TAI","PER","BRA","Cultivated_G","AFR")

het = het[order(match(het$country,sort1)),]
#sort factor levels
het$country <- factor(het$country, levels=sort1)

#######Sort by pop
sort1 <- c("HND","CRI","PAN","COL","TAI","PER","BRA","Cultivated_G","AGO","CMR","GAM","GHA","GUI","MDG","NGA","SEN","SLE","TZA","ZRE")
het = het[order(match(het$pop,sort1)),]
#sort factor levels
het$pop <- factor(het$pop, levels=sort1)


#####BOXPLOTS FOR HETEROZYGOSITY ESTIMATES
library(gplots)
boxplot2(F ~ country, data=het, las=3, ylab="Inbreeding Coefficient (F)", cex.axis=0.75, ylim=c(-1.5,1))
boxplot2(F ~ pop, data=het, las=3, ylab="Inbreeding Coefficient (F)", cex.axis=0.75, ylim=c(-1.5,1))






