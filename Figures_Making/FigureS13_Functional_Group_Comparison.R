rm(list=ls())
library(tidyverse)


# load the tree-level GFBi data
load("**/Global_Analysis/Data/GFBI_Lab_Full.R")

# read in species level mycorrhizal association information from TRY
TRYdata_EM_AM_Combined <- readRDS("**/Global_Analysis/Data/TRYdata_EM_AM.rds")

####---Extract Mycorrhizal type data for GFBI----####

SppLst1 <- unique(GFBI_Df1$SPCD)
Len1 <- length(SppLst1)
Myco_explore_Df <- as.data.frame(cbind(SppName = SppLst1,  
                                       Myco_status = rep(NA, Len1),
                                       GenName = rep(NA, Len1)))
Myco_explore_Df$GenName <- sapply(1:Len1, function(x) strsplit(Myco_explore_Df$SppName[x], split=" ")[[1]][1])
Myco_Name <- c("AM", "EM")

# Species level search
for(i in 1:Len1){
  tic()
  MC_ts_Vector <- TRYdata_EM_AM_Combined$Myco_status[TRYdata_EM_AM_Combined$SppName==Myco_explore_Df$SppName[i]]
  MC_ts_Table <- table(MC_ts_Vector) 
  if(length(MC_ts_Table)==2){  
    Myco_explore_Df$Myco_status[i] <- Myco_Name[which(MC_ts_Table==max(MC_ts_Table))[1]]
  }else if(length(MC_ts_Table)==1){
    Myco_explore_Df$Myco_status[i] <- MC_ts_Vector[1]
  }
  toc()
}

# Genus level search
tic()
for(i in 1:Len1){
  if(is.na(Myco_explore_Df$Myco_status[i])){
    MC_ts_Vector <- Myco_explore_Df$Myco_status[Myco_explore_Df$GenName==Myco_explore_Df$GenName[i]]
    MC_ts_Table <- table(MC_ts_Vector) 
    if(length(MC_ts_Table)==2){  
      Myco_explore_Df$Myco_status[i] <- Myco_Name[which(MC_ts_Table==max(MC_ts_Table))[1]]
    }else if(length(MC_ts_Table)==1){
      Myco_explore_Df$Myco_status[i] <- MC_ts_Vector[1]
    }
  }
}
toc()


## add all species/genus level information back to GFBi dataset
GFBI_Df1$Myco_status <- NA
for(i in 1:Len1){
  tic()
  if(!is.na(Myco_explore_Df$Myco_status[i])){
    GFBI_Df1$Myco_status[GFBI_Df1$SPCD==Myco_explore_Df$SppName[i]] <- Myco_explore_Df$Myco_status[i]
  }
  toc()
}

GFBI_Df1 <- GFBI_Df1%>%mutate(Myco_known=ifelse(is.na(Myco_status),0,1), 
                              EM_Status=ifelse(Myco_status=="EM", 1, 0),
                              AM_Status=ifelse(Myco_status=="AM", 1, 0))
GFBI_Df1 <- GFBI_Df1%>%mutate(Myco_BA=BA*Myco_known, EM_BA=BA*EM_Status, AM_BA=BA*AM_Status)
save(GFBI_Df1, file="**/Global_Analysis/Data/GFBI_Lab_Full.R")



spp_dominance <- function(x){
  TT <- table(x)
  maxProp <- max(TT)/sum(TT)
  return(maxProp)
}


####---Aggregate the tree-level data to plot-level data----####
tic()
GFBI_Df2_LP_LT_Myco <- GFBI_Df1%>%group_by(PLT)%>%summarise(LAT = mean(LAT, na.rm=T),
                                                       LON = mean(LON, na.rm=T),
                                                       #BA=sum(BA, na.rm = T), 
                                                       #LP_BA=sum(LP_BA, na.rm = T),
                                                       #LT_BA=sum(LP_BA, na.rm = T),
                                                       Myco_BA=sum(Myco_BA, na.rm = T),
                                                       #EV_BA=sum(EV_BA, na.rm = T), # plot basal area of evergreen trees
                                                       #BL_BA=sum(BL_BA, na.rm = T), # plot basal area of broadleaf trees
                                                       EM_BA=sum(EM_BA, na.rm = T), # plot basal area of Ectomycorrhiza-associated trees
                                                       #spp_maxprop=spp_dominance(SPCD),
                                                       #spp.count=length(unique(SPCD)),
                                                       #count=sum(count, na.rm=T)
                                                       )
toc()


GFBI_Df2_Aggregated_Full_2 <- readRDS("**/Global_Analysis/Data/GFBI_Df2_aggregated_New.rds")

tic()
GFBI_Df2_Aggregated_Full_2 <- merge(GFBI_Df2_Aggregated_Full_2, GFBI_Df2_LP_LT_Myco[,c("PLT", "Myco_BA", "EM_BA")], by="PLT")
toc()

GFBI_Df2_Aggregated_Full_2 <- GFBI_Df2_Aggregated_Full_2 %>% mutate(relEM=EM_BA/BA, relBL=BL_BA/BA)
GFBI_Df2_Aggregated_Full_2 <- GFBI_Df2_Aggregated_Full_2 %>% mutate(relNL=1-relBL)

GFBI_Df2_Aggregated_Full_2_filtered <- GFBI_Df2_Aggregated_Full_2[GFBI_Df2_Aggregated_Full_2$LP_BA/GFBI_Df2_Aggregated_Full_2$BA>=0.9 & GFBI_Df2_Aggregated_Full_2$LT_BA/GFBI_Df2_Aggregated_Full_2$BA>=0.9 & GFBI_Df2_Aggregated_Full_2$Myco_BA/GFBI_Df2_Aggregated_Full_2$BA>=0.9, ]


# we accounted for the potential impact of human management on the 
# establishment of monocultures by only keeping plots 
# 1) with at least ten trees, 
# 2) with at least two species, 
# 3) in which the relative basal area of the tree species with the largest cumulative basal area in the plot was smaller than 75%
GFBI_Df2_Aggregated_Full_2_filtered <- GFBI_Df2_Aggregated_Full_2_filtered[GFBI_Df2_Aggregated_Full_2_filtered$spp_maxprop<=0.75 & GFBI_Df2_Aggregated_Full_2_filtered$count>=10 & GFBI_Df2_Aggregated_Full_2_filtered$spp.count>=2,]



### density plot
# Figure S13C
ggplot(GFBI_Df2_Aggregated_Full_2_filtered, aes(x=relEV, y=relEM))+
  geom_hex()+
  # geom_smooth(method="lm",se=F)+
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),label.x.npc = 0.5, label.y.npc = "bottom", method = "pearson",label.sep = "/n", size = 5,color="black")+
  # # geom_smooth(method='lm',size=2, color="red", alpha=0.5)+
  stat_smooth(geom="line", alpha=0.7, size=2, color="black", span=0.5, method = "lm", linetype="dashed") +
  # scale_color_viridis(name="Density")+
  scale_fill_gradientn(colors=rev(rainbow(9)[1:6]))+
  # scale_fill_continuous(type = "viridis") +
  # scale_fill_manual(values=rainbow(9)) +
  labs(x = "Relative Evergreen Abundance", y = "Relative EM Trees Abundance") +
  theme_classic()+
  # scale_x_log10()+
  theme(axis.text.x=element_text(size=12, face="bold"),
        axis.text.y=element_text(size=12, face="bold"),
        axis.title.x=element_text(size=16,face="bold"),
        axis.title.y=element_text(size=16,face="bold"),
        legend.text = element_text(colour = "black", size = 10, face = "bold"),
        legend.title = element_text(colour = "black", size = 12, face = "bold"))

# Figure S13D
ggplot(GFBI_Df2_Aggregated_Full_2_filtered, aes(x=relEV, y=relNL))+
  geom_hex()+
  # geom_smooth(method="lm",se=F)+
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),label.x.npc = 0.5, label.y.npc = "bottom", method = "pearson",label.sep = "/n", size = 5,color="black")+
  # # geom_smooth(method='lm',size=2, color="red", alpha=0.5)+
  stat_smooth(geom="line", alpha=0.7, size=2, color="black", span=0.5, method = "lm", linetype="dashed") +
  # scale_color_viridis(name="Density")+
  scale_fill_gradientn(colors=rev(rainbow(9)[1:6]))+
  # scale_fill_continuous(type = "viridis") +
  # scale_fill_manual(values=rainbow(9)) +
  labs(x = "Relative Evergreen Abundance", y = "Relative Needleleaf Abundance") +
  theme_classic()+
  # scale_x_log10()+
  theme(axis.text.x=element_text(size=12, face="bold"),
        axis.text.y=element_text(size=12, face="bold"),
        axis.title.x=element_text(size=16,face="bold"),
        axis.title.y=element_text(size=16,face="bold"),
        legend.text = element_text(colour = "black", size = 10, face = "bold"),
        legend.title = element_text(colour = "black", size = 12, face = "bold"))


####---Aggregate plot-level information to get grid-level index for comparison----####
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}


tic()
GFBI_Df3_Scaled_10min <- GFBI_Df2_Aggregated_Full_2_filtered%>%
  group_by(Cluster_10min)%>%summarise(latitude=mean(latitude),longitude=mean(longitude), 
                                      dip_Stat=dip.test(relEV)$statistic,
                                      skewness=skewness(relEV), 
                                      dip_Stat_Myco=dip.test(relEM)$statistic,
                                      skewness_Myco=skewness(relEM), 
                                      dip_Stat_LT=dip.test(relNL)$statistic,
                                      skewness_LT=skewness(relNL), 
                                      SampleSize=sum(SampleSize)
  )
toc()


GFBI_Df3_Scaled_10min <- GFBI_Df3_Scaled_10min %>% mutate(dip_Stat_AD=sqrt(SampleSize)*dip_Stat,
                                                          dip_Stat_AD_Myco=sqrt(SampleSize)*dip_Stat_Myco,
                                                          dip_Stat_AD_LT=sqrt(SampleSize)*dip_Stat_LT)

GFBI_Df3_Scaled_10min <- GFBI_Df3_Scaled_10min %>% mutate(BI=-exp(-(dip_Stat_AD)^2*6)*sign(skewness),
                                                          BI_Myco=-exp(-(dip_Stat_AD_Myco)^2*6)*sign(skewness_Myco),
                                                          BI_LT=-exp(-(dip_Stat_AD_LT)^2*6)*sign(skewness_LT))

GFBI_Df3_Scaled_10min <- GFBI_Df3_Scaled_10min%>%drop_na()
Df_Map_Comparison <- as.data.frame(cbind(values_LP=GFBI_Df3_Scaled_10min$BI, values_Myco=GFBI_Df3_Scaled_10min$BI_Myco,
                                         values_LT=GFBI_Df3_Scaled_10min$BI_LT))


# Figure S13A
Df_Map_Comparison <- Df_Map_Comparison%>%mutate(density=get_density(values_LP, values_Myco, n=100))
ggplot(Df_Map_Comparison, aes(x=values_LP, y=values_Myco, color=density))+
  geom_point()+
  # geom_smooth(method="lm",se=F)+
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),label.x.npc = 0.5, label.y.npc = "bottom", method = "spearman",label.sep = "/n", size = 5)+
  # # geom_smooth(method='lm',size=2, color="red", alpha=0.5)+
  stat_smooth (geom="line", alpha=0.7, size=2, color="black", span=0.5, method = "lm", linetype="dashed") +
  # scale_color_viridis(name="Density")+
  scale_color_gradientn(name="Density", colors=rev(rainbow(9)[1:6]))+
  labs(x = "BI for forest leaf phenology types", y = "BI for forest mycorrhizal types") +
  theme_classic()+
  # scale_x_log10()+
  theme(axis.text.x=element_text(size=12, face="bold"),
        axis.text.y=element_text(size=12, face="bold"),
        axis.title.x=element_text(size=16,face="bold"),
        axis.title.y=element_text(size=16,face="bold"),
        legend.text = element_text(colour = "black", size = 10, face = "bold"),
        legend.title = element_text(colour = "black", size = 12, face = "bold"))

# Figure S13B
Df_Map_Comparison <- Df_Map_Comparison%>%mutate(density=get_density(values_LP, values_LT, n=100))
ggplot(Df_Map_Comparison, aes(x=values_LP, y=values_LT, color=density))+
  geom_point()+
  # geom_smooth(method="lm",se=F)+
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),label.x.npc = 0.5, label.y.npc = "bottom", method = "spearman",label.sep = "/n", size = 5)+
  # # geom_smooth(method='lm',size=2, color="red", alpha=0.5)+
  stat_smooth (geom="line", alpha=0.7, size=2, color="black", span=0.5, method = "lm", linetype="dashed") +
  # scale_color_viridis(name="Density")+
  scale_color_gradientn(name="Density", colors=rev(rainbow(9)[1:6]))+
  labs(x = "BI for forest leaf phenology types", y = "BI for forest leaf forms") +
  theme_classic()+
  # scale_x_log10()+
  theme(axis.text.x=element_text(size=12, face="bold"),
        axis.text.y=element_text(size=12, face="bold"),
        axis.title.x=element_text(size=16,face="bold"),
        axis.title.y=element_text(size=16,face="bold"),
        legend.text = element_text(colour = "black", size = 10, face = "bold"),
        legend.title = element_text(colour = "black", size = 12, face = "bold"))