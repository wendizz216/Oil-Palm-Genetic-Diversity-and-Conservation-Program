##########################################
#Making Boxplot figures in ggplot for heterozygosity estimates for Latin American and Oil palm species
#Heterozygosity estimates were calculated using VCFtools
##########################################

library(ggplot2)

het <- read.csv("/Users/wvu/Dropbox/Github/Hetero/HETEROZYGOSITY_OILPALM.csv")

head(het)

p <- ggplot(het, aes(x=POP, y=F_ADJUST, fill=SPECIES)) + geom_boxplot() + labs(titles="Genetic Diversity Estimates",x="Oil Palm Populations",y="1 - F (Inbreeding Coefficient)") + scale_fill_brewer(palette="Dark2")

#Resizing axis titles and text
p + theme(
		plot.title = element_text(size=30, face="bold"), 
		axis.title.x=element_text(size=20, face="bold"),
		axis.title.y=element_text(size=20, face="bold"),
		axis.text.x=element_text(size=12, face="bold"),
		axis.text.y=element_text(size=12, face="bold"),
		legend.title=element_text(size=20, face="bold"),
		legend.text=element_text(size=15,face="bold")) 



