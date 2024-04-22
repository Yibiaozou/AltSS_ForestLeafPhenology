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

####----FIA_US: GAMLSS zero adjusted Poisson (ZAP) distribution----####
storage.dir <- '**/Continental_Analysis/FIA_US/'
setwd(storage.dir)
source('Project_Functions/paths.r')

AX <- read.csv("**/Source_Data/SourceData_Fig1A-C.csv")

# For FIA we already remove managed plot, therefore here we only remove small plots with less than 10 trees
AX <- AX[AX$count>=10,]
AX <- AX %>% mutate(relEV_Density=ev.density/(ev.density+de.density))
hist(AX$relEV_Density,breaks=20)
hist(AX$relEV,breaks=20)


# the data
AX2<-AX[,c("latitude","longitude","county.id","ev.density","de.density","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10", "relEV")]
AX2<-AX2[complete.cases(AX2),]

### Get spatial cluster of each plot. we partitioned the global forest zones using a ‘fishing net’ with 10 arc-min (~20km) grid size.
fishNet <- raster() # create an empty raster template 
res(fishNet) <- 1/6 # set the resolution of the template as 1/6 degree. 
values(fishNet) <- 1:ncell(fishNet) # allocate values to the 'fishing net' raster
# extract spatial cluster for each plot
AX2$Cluster <- raster::extract(fishNet, AX2[,c("longitude", "latitude")])

AzapMx <- readRDS("Data/FIA_zapMx_EV_FSD.rds")
BzapMx <- readRDS("Data/FIA_zapMx_DE_FSD.rds")

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

aDest_Bi<-ggplot_build(aL_Bi)$data[[1]]
ASest_Bi<-data.frame(x=aDest_Bi$x,count=aDest_Bi$count,relf=aDest_Bi$count/sum(aDest_Bi$count),id=i)

## Figure 1B
#this plot gives us the comparison. black points and confidence intervals give the binlength expected under the nullmodels
ggplot(ASx_Bi,aes(x=x,y=med))+geom_point()+geom_bar(data=ASest_Bi,aes(y=relf,fill=..x..),show.legend = FALSE,stat="identity",alpha=0.8)+geom_errorbar(aes(ymin=low,ymax=high))+theme_classic()+
  scale_fill_gradient("",low='red', high='blue')+
  ylab("Relative Frequency")+xlab("Relative Abundance Evergreen Trees")+
  theme(axis.text.x=element_text(size=12, face="bold"),
        axis.text.y=element_text(size=12, face="bold"),
        axis.title.x=element_text(size=16,face="bold"),
        axis.title.y=element_text(size=16,face="bold"),
        legend.text = element_text(colour = "black", size = 10, face = "bold"),
        legend.title = element_text(colour = "black", size = 12, face = "bold"))

## Figure 1C
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


####----GFBI: GAMLSS zero adjusted Poisson (ZAP) distribution----####
# GFBI_Df2_Aggregated_Full_2_filtered <- readRDS("**/Global_Analysis/Data/GFBI_Df2_Aggregated_Full_2_FSD_Fig1.rds")
GFBI_Df2_Aggregated_Full_2_filtered <- read.csv("**/Source_Data/SourceData_Fig1D-F.csv")

AzapMx <- readRDS("**/Global_Analysis/Data/GFBI_zapMx_EV_FSD.rds")
BzapMx <- readRDS("**/Global_Analysis/Data/GFBI_zapMx_DE_FSD.rds")

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
  G1x_BI <- G1x_BI %>% drop_na()
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
aL_Bi<-ggplot(GFBI_Df2_Aggregated_Full_2_filtered,aes(x=EV_Density/(DE_Density+EV_Density),fill=..x..))+geom_histogram(binwidth = binwidthx)+xlab("relEV")+theme_classic()+
  scale_fill_gradient(low='red', high='blue')+ggtitle(paste0("model"))+ theme(legend.title = element_blank())

aDest_Bi<-ggplot_build(aL_Bi)$data[[1]]
ASest_Bi<-data.frame(x=aDest_Bi$x,count=aDest_Bi$count,relf=aDest_Bi$count/sum(aDest_Bi$count),id=i)

## Figure 1E
#this plot gives us the comparison. black points and confidence intervals give the binlength expected under the nullmodels
ggplot(ASx_Bi,aes(x=x,y=med))+geom_point()+geom_bar(data=ASest_Bi,aes(y=relf,fill=..x..),show.legend = FALSE,stat="identity",alpha=0.8)+geom_errorbar(aes(ymin=low,ymax=high))+theme_classic()+
  scale_fill_gradient("",low='red', high='blue')+
  ylab("Relative Frequency")+xlab("Relative Abundance Evergreen Trees")+
  theme(axis.text.x=element_text(size=12, face="bold"),
        axis.text.y=element_text(size=12, face="bold"),
        axis.title.x=element_text(size=16,face="bold"),
        axis.title.y=element_text(size=16,face="bold"),
        legend.text = element_text(colour = "black", size = 10, face = "bold"),
        legend.title = element_text(colour = "black", size = 12, face = "bold"))

## Figure 1F
#this plot gives us the comparison of spearman's coefficient between observation and null model predictions.
corsample<-cor(GFBI_Df2_Aggregated_Full_2_filtered$EV_Density,GFBI_Df2_Aggregated_Full_2_filtered$DE_Density,method="spearman")
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
