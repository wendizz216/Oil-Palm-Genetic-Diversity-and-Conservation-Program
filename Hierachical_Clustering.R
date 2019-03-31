library(poppr) 
library(adegenet)

snp <-  read.table("/Users/wvu/Dropbox/Oil_Palm/RAD_SequencingResults_Dec2017/American_diversity/Nei_GeneticDist/SOL_OILPALM.Eguineensis9.1.VF.SNP.DP8.90.miA02.thin100.hcpc.tp.bra.cultG.vcf.tab",header=T,sep='\t')

##########COMPUTE SIMILARITY MEASURES: I will use calculate Nei's genetic similarity for this data set.
#Transform SNP Data frame: rows = sample ID and cols = snps
rownames(snp) <- paste(snp[,1],snp[,2],sep='_')
snp <- snp[,-c(1:3)]
snp <- as.data.frame(t(snp))

#convert snp data frame to genind object:  https://cran.r-project.org/web/packages/adegenet/adegenet.pdf
snp.gen <- df2genind(snp, sep='/',ploidy=2, NA.char="./.",ncode=1)

#decide how to deal with missing data: four options: https://cran.r-project.org/web/packages/poppr/poppr.pdf
#snp.gen <- missingno(snp.gen, type = 'mean')

#REPLACE missing data with zeros is best for this data set
snp.gen <- missingno(snp.gen, type = '0',cutoff=0.10)
nei.dist(snp.gen) -> snp.nei

###########IMPLEMENT LINKAGE FUNCTION: choose linkage method to generate cluster tree
#"complete" is the best method to cluster this data based correlation of the cophenetic distances and original distances as well as genetic relatedness expectations.

hc <- hclust(as.dist(snp.nei),method="ward.D2") #on the scale of distances
res.coph <- cophenetic(hc)
cor(as.dist(snp.nei),res.coph) #0.86

hc <- hclust(as.dist((snp.nei)^2),method="ward.D") #on the scale of distances squared
res.coph <- cophenetic(hc)
cor(as.dist(snp.nei),res.coph) #0.79

hc <- hclust(as.dist((snp.nei)^2),method="complete")
res.coph <- cophenetic(hc)
cor(as.dist(snp.nei),res.coph) #0.8909515

############ Cut the dendrogram into different groups

library("factoextra")
fviz_dend(hc, cex=0.5)

#based on dengram, there are 4 distinct genetic groups
grp <- cutree(hc, k=4)
table(grp) #number of members in each clustered group

fviz_dend(hc, cex=0.5, k=4, color_labels_by_K=TRUE, rect=T, k_colors=c("#1B9E77", "#D95F02", "#7570B3", "#E7298A"))

fviz_dend(hc, k=4, k_colors = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A"), type="phylogenic", repel = T)

############# K-MEANS CLUSTERING 

snp.mat <- snp.gen$tab

#Scale data
snp.scaled <- scale(snp.mat)

#K-means clustering
km.snps <- kmeans(snp.scaled, 4, nstart=30)

#VISUALIZE K-MEANS CLUSTERING

fviz_cluster(km.snps, snp.mat, ellipse.type="norm") + scale_colour_manual(values =c("#1B9E77", "#D95F02", "#7570B3", "#E7298A")) + scale_fill_manual(values=c("#1B9E77", "#D95F02", "#7570B3", "#E7298A"))

fviz_cluster(km.snps, snp.mat, ggtheme = theme_minimal()) + scale_colour_manual(values =c("#1B9E77", "#D95F02", "#E7298A", "#7570B3")) + scale_fill_manual(values=c("#1B9E77", "#D95F02", "#E7298A","#7570B3"))


