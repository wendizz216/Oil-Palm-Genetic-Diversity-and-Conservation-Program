#######IDENTIFY MINIMUM SAMPLE SIZE TO CAPTURE AT LEAST 90% OF GENETIC VARIATION: polymorphic sites sampling 

p <- read.table("/Users/wvu/Dropbox/Oil_Palm/RAD_SequencingResults_Dec2017/Diversity/GeneticVariation_Sampling.txt",header=T,sep='\t')

#Genetic Diversity Index: %polymorphism @ N=x / Total %polymorphism N=154
mean <- aggregate(p$prop_SegSNPs,list(N_size=p$N_size),mean)
diversity_index <- data.frame(mean,div_index=mean[,2]/mean[31,2])
s.var <- aggregate(p$prop_SegSNPs,list(N_size=p$N_size),var)

#Looks like the data is non-linear and follows a Michaelis-Menten function
y <- diversity_index[,3]
x <- diversity_index[,1]
plot(x,y)

#Fitting a Michaelis-Menten function to data

fit <- nls(y ~ SSmicmen(x, Vm, K), data = data.frame(x,y))
summary(fit)

#Goodness of fit 
cor(y,predict(fit)) #0.9970813

#######RESULT Model Function: y = Vm*x/(K+x)  ---> y = 1.12x/(19.04 + x)

#Plot nonlinear regression curve to sampling data
plot(y~x, axes = F, ylab = "Genetic Diversity Index", xlab="Sample Size",pch=19)
axis(side = 1, at=unique(x))
axis(side = 2, at=c(0.1,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0))
lines(seq(0, 150, length.out = 100), predict(fit, newdata= data.frame(x=seq(0, 150, length.out=100))))
abline(a=0.91,b=0.80,col="red")
#arrows(x, y-sdev[,2], x, y+sdev[,2], length=0.05, angle=90, code=3)

#Useful NOTE: List of self-starting functions in R to estimate the start parameters for 10 models
#1	SSasymp       asymptotic regression models
#2	SSasympOff    asymptotic regression models with an offset
#3	SSasympOrig   asymptotic regression models through the origin
#4	SSbiexp       biexponential models
#5	SSfol         first-order compartment models
#6	SSfpl         four-parameter logistic models
#7	SSgompertz    Gompertz growth models
#8	SSlogis       logistic models
#9 	SSmicmen      Michaelis–Menten models
#10 SSweibull     Weibull growth curve models

https://www.zoology.ubc.ca/~schluter/R/fit-model/
https://dataconomy.com/2017/08/nonlinear-least-square-nonlinear-regression-r/
