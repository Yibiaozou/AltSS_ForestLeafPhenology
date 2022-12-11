#Subset p1 and p2 to get a subset of data for eastern US, which later we used for ecoregion-level demographic analysis
rm(list=ls())

setwd(storage.dir)
source('Project_Functions/paths.r')
library(ggplot2)
library(ggalt)
library(data.table)

#set output path.----
output.path <- Product_2.subset.path

#load data.----
p2 <- data.table(readRDS(Product_2.path))
states <- read.csv('required_products_utilities/FIA_state_codes_regions.csv')
states.east <- states[states$east_plus == 1,]
states$northeast <- ifelse(states$state_name %in% c('Connecticut','New York','Massachusetts','Rhode Island','New Hampshire','Vermont','Maine'),1,0)

#Some final data cleaning (should be moved?).-----
#Would be nice to move these to full data filtering script (#2).
p2 <- p2[p2$REMPER >=4.9 & p2$REMPER <= 5.1,]
p2 <- p2[p2$n.trees >  5,]
p2 <- p2[DIA.cm > 0,]

#Complete cases of environmental covariates.
env.drop <- p2[,c('PLT_CN','mat','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10')]
env.drop <- env.drop[!complete.cases(env.drop),]
env.drop <- unique(env.drop$PLT_CN)
p2 <- p2[!(p2$PLT_CN %in% env.drop),]

#Grab a plot table with PLT_CN, lat-lon and STATECD.----
d <- p2[,.(PLT_CN,LAT,LON,STATECD,relEM,REMPER)]
setkey(d, 'PLT_CN')
d <- unique(d)

#Subset.----
#subset to eastern US.
d <- d[d$STATECD %in% states.east$STATECD,]
#NORTHEAST FOR TESTING.
p2 <- p2[PLT_CN %in% d$PLT_CN,]

#Save output.----
saveRDS(p2, output.path)
