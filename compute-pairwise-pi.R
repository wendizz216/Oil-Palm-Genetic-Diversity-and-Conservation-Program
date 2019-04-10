##this was done to run a PCA to understand genetic diversity among the 96 individuals with a total of 64K SNPs and potentially explain why there are two distinct clusters.


x <- matrix( rbinom(40,size=2,prob=0.5), ncol=4 )/2
x=snp.pos

pimat <- matrix( NA, nrow=ncol(x), ncol=ncol(x) )

for (i in 1:ncol(x)) {
    for (j in 1:i) {
        if (i==j) {
            pimat[i,i] <- 0.5 * sum( x[,i] == 0.5 ) / sum( !is.na(x[,i]) )
        } else {
            pimat[i,j] <- pimat[j,i] <- ( sum( abs( x[,i] - x[,j] ), na.rm=TRUE ) + 0.5 * sum( x[,i]==x[,j] & x[,i]==0.5, na.rm=TRUE ) ) / sum( !is.na(x[,i]) & !is.na(x[,j]) )
        }
    }
}

pimat<-as.matrix(as.dist(pimat))

pi.eig<-eigen(pimat)

plot(pi.eig$vectors[,1],pi.eig$vectors[,2],xlab="PC1",ylab="PC2",col=ifelse(str.sort.geno[,8]=="hap1","black","red"),pch=19)