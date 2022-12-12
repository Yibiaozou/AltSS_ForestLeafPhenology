rm(list=ls())

library(tidyverse)
library(ggpubr)

## load plot-level aggregated forest data
GFBI_Df2_Aggregated_Full_2 <- readRDS("**/Global_Analysis/Data/GFBI_Df2_aggregated_New.rds")

# function to get density of one variable along another variable, for visualization
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

# we accounted for the potential impact of human management on the 
# establishment of monocultures by only keeping plots 
# 1) with at least ten trees, 
# 2) with at least two species, 
# 3) in which the relative basal area of the tree species with the largest cumulative basal area in the plot was smaller than 75%
GFBI_Df2_Aggregated_Full_2_filtered <- GFBI_Df2_Aggregated_Full_2[GFBI_Df2_Aggregated_Full_2$spp_maxprop<=0.75 & GFBI_Df2_Aggregated_Full_2$count>=10 & GFBI_Df2_Aggregated_Full_2$spp.count>=2,]

# get density between the two types of relEV
GFBI_Df2_Aggregated_Full_2_filtered <- GFBI_Df2_Aggregated_Full_2_filtered%>%drop_na()%>%mutate(density=get_density(relEV, relEV_Density, n=100))

# Figure S7
ggplot(GFBI_Df2_Aggregated_Full_2_filtered, aes(x=relEV_Density, y=relEV, color=density))+
  geom_point()+
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),label.x.npc = 0.5, label.y.npc = "top", method = "spearman",label.sep = "/n", size = 5)+
  stat_smooth (geom="line", alpha=0.7, size=2, color="black", span=0.5, method = "lm") +
  scale_color_gradientn(name="Density", colors=rev(rainbow(9)[1:6]))+
  labs(x = "Individual-based relEV", y = "Area-based relEV") +
  theme_classic()+
  theme(axis.text.x=element_text(size=12, face="bold"),
        axis.text.y=element_text(size=12, face="bold"),
        axis.title.x=element_text(size=20,face="bold"),
        axis.title.y=element_text(size=20,face="bold"),
        legend.text = element_text(colour = "black", size = 10, face = "bold"),
        legend.title = element_text(colour = "black", size = 12, face = "bold"))