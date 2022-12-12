rm(list=ls())

library(tidyverse)
library(tictoc)
# library(geosphere)
library(mousetrap)
# library(diptest)
library(moments)
# library(raster)
# library(BimodalIndex)
library(LaplacesDemon)
library(MASS)
library(ggpubr)
# library(rgdal)
library(data.table)
# library("rnaturalearth")
# library("rnaturalearthdata")
# library(rgeos)
library(matrixStats)
library(factoextra)
library(ggfortify)
# library(caret)
# library(disttree)
# library(ranger)
library(gamlss)
library(gamlss.add)

####----FunDivEUROPE_EU: GAMLSS zero adjusted Poisson (ZAP) distribution----####
AzapMx <- readRDS("**Continental_Analysis/FunDivEurope_EU/Data/EU_zapMx_EV_FSD.rds")
BzapMx <- readRDS("**Continental_Analysis/FunDivEurope_EU/Data/EU_zapMx_DE_FSD.rds")
iterations=1000
GebnMx_Bi <- numeric(iterations)
KebnMx_Bi <- numeric(iterations)
binwidthx = 0.05
AS_Bi<-list()



for(i in 1:iterations){
  tic()
  #create simulated counts data frame
  Gx_BI<-data.frame(A=AzapMx[i,],B=BzapMx[i,])
  
  Gx_BI$D<-Gx_BI$A/rowSums(Gx_BI[,1:2]) 
  
  
  #filter to have only plots with min tree ==10
  G1x_BI<-Gx_BI[rowSums(Gx_BI[,1:2])>=10,]
  
  # spearman correlation of counts
  GebnMx_Bi[i]<-cor(G1x_BI$A,G1x_BI$B,method="spearman")
  
  
  
  #here we create the bimodality histogram based on our null model draws in each iteration 
  a2x_BI<-ggplot(G1x_BI,aes(x=D,fill=..x..))+geom_histogram(binwidth = binwidthx)
  aDx_BI<-ggplot_build(a2x_BI)$data[[1]]
  AS_Bi[[i]]<-data.frame(x=aDx_BI$x,count=aDx_BI$count,relf=aDx_BI$count/sum(aDx_BI$count),id=i)
  toc()
}

#here we combine the simulated bimodality histograms to generate the confidence interval under the null for each bin
ASx_Bi<-rbindlist(AS_Bi)
ASx_Bi<-ASx_Bi %>% group_by(x)%>%summarise(low=quantile(relf,0.025),high=quantile(relf,0.975),med=quantile(relf,0.5))
# we compare to the observed bimodality histogram
aL_Bi<-ggplot(AX2,aes(x=ev.density/(de.density+ev.density),fill=..x..))+geom_histogram(binwidth = binwidthx)+xlab("relEV")+theme_classic()+
  scale_fill_gradient(low='red', high='blue')+ggtitle(paste0("model"))+ theme(legend.title = element_blank())
#+ylim(0,15000)
aDest_Bi<-ggplot_build(aL_Bi)$data[[1]]
ASest_Bi<-data.frame(x=aDest_Bi$x,count=aDest_Bi$count,relf=aDest_Bi$count/sum(aDest_Bi$count),id=i)

# Figure S3B
#this plot gives us the comparison. black points and confidence intervals give the binlength expected under the nullmodels
ggplot(ASx_Bi,aes(x=x,y=med))+geom_point()+geom_bar(data=ASest_Bi,aes(y=relf,fill=..x..),show.legend = FALSE,stat="identity",alpha=0.6)+geom_errorbar(aes(ymin=low,ymax=high))+theme_classic()+
  scale_fill_gradient("",low='red', high='blue')+
  ylab("Relative Frequency")+xlab("Relative Abundance Evergreen Trees")+
  theme(axis.text.x=element_text(size=12, face="bold"),
        axis.text.y=element_text(size=12, face="bold"),
        axis.title.x=element_text(size=16,face="bold"),
        axis.title.y=element_text(size=16,face="bold"),
        legend.text = element_text(colour = "black", size = 10, face = "bold"),
        legend.title = element_text(colour = "black", size = 12, face = "bold"))

# Figure S3C
#this plot gives us the comparison of spearman's coefficient between observation and null model predictions.
corsample<-cor(AX2$ev.density,AX2$de.density,method="spearman")
cordistribution<-data.frame(cor=GebnMx_Bi)
ggplot()+geom_histogram(data=cordistribution,aes(x=cor),bins=200, fill="black")+
  geom_vline(xintercept = corsample,color="red", size=2)+
  xlim(c(-0.75,0.1))+
  xlab("Spearman's rank correlation coefficient")+
  ylab("Count")+
  theme_classic()+
  theme(axis.text.x=element_text(size=12, face="bold"),
        axis.text.y=element_text(size=12, face="bold"),
        axis.title.x=element_text(size=16,face="bold"),
        axis.title.y=element_text(size=16,face="bold"),
        legend.text = element_text(colour = "black", size = 10, face = "bold"),
        legend.title = element_text(colour = "black", size = 12, face = "bold"))

####----FunDivEUROPE_EU: GAMLSS zero negative binomial (ZANBI) distribution----####
AbnMx <- readRDS("**Continental_Analysis/FunDivEurope_EU/Data/EU_bnMx_EV_FSD.rds")
BbnMx <- readRDS("**Continental_Analysis/FunDivEurope_EU/Data/EU_bnMx_DE_FSD.rds")
iterations=1000
GebnMx_Bi <- numeric(iterations)
KebnMx_Bi <- numeric(iterations)
binwidthx = 0.05
AS_Bi<-list()



for(i in 1:iterations){
  tic()
  #create simulated counts data frame
  Gx_BI<-data.frame(A=AbnMx[i,],B=BbnMx[i,])
  
  Gx_BI$D<-Gx_BI$A/rowSums(Gx_BI[,1:2]) 
  
  
  #filter to have only plots with min tree ==10
  G1x_BI<-Gx_BI[rowSums(Gx_BI[,1:2])>=10,]
  
  # spearman correlation of counts
  GebnMx_Bi[i]<-cor(G1x_BI$A,G1x_BI$B,method="spearman")
  
  
  
  #here we create the bimodality histogram based on our null model draws in each iteration 
  a2x_BI<-ggplot(G1x_BI,aes(x=D,fill=..x..))+geom_histogram(binwidth = binwidthx)
  aDx_BI<-ggplot_build(a2x_BI)$data[[1]]
  AS_Bi[[i]]<-data.frame(x=aDx_BI$x,count=aDx_BI$count,relf=aDx_BI$count/sum(aDx_BI$count),id=i)
  toc()
}

#here we combine the simulated bimodality histograms to generate the confidence interval under the null for each bin
ASx_Bi<-rbindlist(AS_Bi)
ASx_Bi<-ASx_Bi %>% group_by(x)%>%summarise(low=quantile(relf,0.025),high=quantile(relf,0.975),med=quantile(relf,0.5))
# we compare to the observed bimodality histogram
aL_Bi<-ggplot(AX2,aes(x=ev.density/(de.density+ev.density),fill=..x..))+geom_histogram(binwidth = binwidthx)+xlab("relEV")+theme_classic()+
  scale_fill_gradient(low='red', high='blue')+ggtitle(paste0("model"))+ theme(legend.title = element_blank())
#+ylim(0,15000)
aDest_Bi<-ggplot_build(aL_Bi)$data[[1]]
ASest_Bi<-data.frame(x=aDest_Bi$x,count=aDest_Bi$count,relf=aDest_Bi$count/sum(aDest_Bi$count),id=i)

# Figure S3D
#this plot gives us the comparison of spearman's coefficient between observation and null model predictions.
corsample<-cor(AX2$ev.density,AX2$de.density,method="spearman")
cordistribution<-data.frame(cor=GebnMx_Bi)
ggplot()+geom_histogram(data=cordistribution,aes(x=cor),bins=200, fill="black")+
  geom_vline(xintercept = corsample,color="red", size=2)+
  xlim(c(-0.75,0.1))+
  xlab("Spearman's rank correlation coefficient")+
  ylab("Count")+
  theme_classic()+
  theme(axis.text.x=element_text(size=12, face="bold"),
        axis.text.y=element_text(size=12, face="bold"),
        axis.title.x=element_text(size=16,face="bold"),
        axis.title.y=element_text(size=16,face="bold"),
        legend.text = element_text(colour = "black", size = 10, face = "bold"),
        legend.title = element_text(colour = "black", size = 12, face = "bold"))

####----FIA_US: GAMLSS zero negative binomial (ZANBI) distribution----####
AbnMx <- readRDS("**/Continental_Analysis/FIA_US/Data/FIA_bnMx_EV_FSD.rds")
BbnMx <- readRDS("**/Continental_Analysis/FIA_US/Data/FIA_bnMx_DE_FSD.rds")

iterations=1000

GebnMx_Bi <- numeric(iterations)
KebnMx_Bi <- numeric(iterations)
binwidthx = 0.05
AS_Bi<-list()

for(i in 1:iterations){
  tic()
  #create simulated counts data frame
  Gx_BI<-data.frame(A=AbnMx[i,],B=BbnMx[i,])
  
  Gx_BI$D<-Gx_BI$A/rowSums(Gx_BI[,1:2]) 
  
  
  #filter to have only plots with min tree ==10 
  G1x_BI<-Gx_BI[rowSums(Gx_BI[,1:2])>=10,]
  
  # spearman correlation of counts
  GebnMx_Bi[i]<-cor(G1x_BI$A,G1x_BI$B,method="spearman")
  
  
  
  #here we create the bimodality histogram based on our null model draws in each iteration 
  a2x_BI<-ggplot(G1x_BI,aes(x=D,fill=..x..))+geom_histogram(binwidth = binwidthx)
  aDx_BI<-ggplot_build(a2x_BI)$data[[1]]
  AS_Bi[[i]]<-data.frame(x=aDx_BI$x,count=aDx_BI$count,relf=aDx_BI$count/sum(aDx_BI$count),id=i)
  toc()
}

#here we combine the simulated bimodality histograms to generate the confidence interval under the null for each bin
ASx_Bi<-rbindlist(AS_Bi)
ASx_Bi<-ASx_Bi %>% group_by(x)%>%summarise(low=quantile(relf,0.025),high=quantile(relf,0.975),med=quantile(relf,0.5))

# we compare to the observed bimodality histogram
aL_Bi<-ggplot(AX2,aes(x=ev.density/(de.density+ev.density),fill=..x..))+geom_histogram(binwidth = binwidthx)+xlab("relEV")+theme_classic()+
  scale_fill_gradient(low='red', high='blue')+ggtitle(paste0("model"))+ theme(legend.title = element_blank())

aDest_Bi<-ggplot_build(aL_Bi)$data[[1]]
ASest_Bi<-data.frame(x=aDest_Bi$x,count=aDest_Bi$count,relf=aDest_Bi$count/sum(aDest_Bi$count),id=i)

## Figure S3E
#this plot gives us the comparison of spearman's coefficient between observation and null model predictions.
corsample<-cor(AX2$ev.density,AX2$de.density,method="spearman")
cordistribution<-data.frame(cor=GebnMx_Bi)
ggplot()+geom_histogram(data=cordistribution,aes(x=cor),bins=200, fill="black")+
  geom_vline(xintercept = corsample,color="red", size=2)+
  xlim(c(-0.7,-0.1))+
  xlab("Spearman's rank correlation coefficient")+
  ylab("Count")+
  theme_classic()+
  theme(axis.text.x=element_text(size=12, face="bold"),
        axis.text.y=element_text(size=12, face="bold"),
        axis.title.x=element_text(size=16,face="bold"),
        axis.title.y=element_text(size=16,face="bold"),
        legend.text = element_text(colour = "black", size = 10, face = "bold"),
        legend.title = element_text(colour = "black", size = 12, face = "bold"))


####----GFBi: GAMLSS zero negative binomial (ZANBI) distribution----####
AbnMx <- readRDS("**/Global_Analysis/Data/GFBI_bnMx_EV_sub.rds")
BbnMx <- readRDS("**/Global_Analysis/Data/GFBI_bnMx_DE_sub.rds")
iterations=1000
GebnMx_Bi <- numeric(iterations)
KebnMx_Bi <- numeric(iterations)
binwidthx = 0.05
AS_Bi<-list()

for(i in 1:iterations){
  tic()
  #create simulated counts data frame
  Gx_BI<-data.frame(A=AbnMx[i,],B=BbnMx[i,])
  
  Gx_BI$D<-Gx_BI$A/rowSums(Gx_BI[,1:2]) 
  
  
  #filter to have only plots with min tree ==10 
  G1x_BI<-Gx_BI[rowSums(Gx_BI[,1:2])>=10,]
  
  # spearman correlation of counts
  GebnMx_Bi[i]<-cor(G1x_BI$A,G1x_BI$B,method="spearman")
  
  
  
  #here we create the bimodality histogram based on our null model draws in each iteration 
  a2x_BI<-ggplot(G1x_BI,aes(x=D,fill=..x..))+geom_histogram(binwidth = binwidthx)
  aDx_BI<-ggplot_build(a2x_BI)$data[[1]]
  AS_Bi[[i]]<-data.frame(x=aDx_BI$x,count=aDx_BI$count,relf=aDx_BI$count/sum(aDx_BI$count, na.rm=T),id=i)
  toc()
}

#here we combine the simulated bimodality histograms to generate the confidence interval under the null for each bin
ASx_Bi<-rbindlist(AS_Bi)
ASx_Bi<-ASx_Bi %>% group_by(x)%>%summarise(low=quantile(relf,0.025),high=quantile(relf,0.975),med=quantile(relf,0.5))
# we compare to the observed bimodality histogram
aL_Bi<-ggplot(GFBI_Df2_Aggregated_Full_2_filtered_stratSample_Fig1,aes(x=EV_Density/(DE_Density+EV_Density),fill=..x..))+geom_histogram(binwidth = binwidthx)+xlab("relEV")+theme_classic()+
  scale_fill_gradient(low='red', high='blue')+ggtitle(paste0("model"))+ theme(legend.title = element_blank())
#+ylim(0,15000)
aDest_Bi<-ggplot_build(aL_Bi)$data[[1]]
ASest_Bi<-data.frame(x=aDest_Bi$x,count=aDest_Bi$count,relf=aDest_Bi$count/sum(aDest_Bi$count),id=i)

## Figure S3F
#this plot gives us the comparison of spearman's coefficient between observation and null model predictions.
corsample<-cor(GFBI_Df2_Aggregated_Full_2_filtered_stratSample_Fig1$EV_Density,GFBI_Df2_Aggregated_Full_2_filtered_stratSample_Fig1$DE_Density,method="spearman")
cordistribution<-data.frame(cor=GebnMx_Bi)
ggplot()+geom_histogram(data=cordistribution,aes(x=cor),bins=200, fill="black")+
  geom_vline(xintercept = corsample,color="red", size=2)+
  xlim(c(-0.75, 0.1))+
  xlab("Spearman's rank correlation coefficient")+
  ylab("Count")+
  theme_classic()+
  theme(axis.text.x=element_text(size=12, face="bold"),
        axis.text.y=element_text(size=12, face="bold"),
        axis.title.x=element_text(size=16,face="bold"),
        axis.title.y=element_text(size=16,face="bold"),
        legend.text = element_text(colour = "black", size = 10, face = "bold"),
        legend.title = element_text(colour = "black", size = 12, face = "bold"))