rm(list=ls())
setwd("D:/Zeus/ETH_zurich_MSc/ETHz_S2/Eco_Evo_ETH-D_assistantship/Leaf_Phenology_Project/Code/YibiaoZ_202x_LeafHabit_ColinA/altSS_forest_mycorrhizas-master/")
library(mgcv)
library(data.table)
source('paths.r')
source('project_functions/predict_gam_well.r')
library(tidyverse)
library(raster)

#function to draw randomly from multivariate normal distribution.-----
rmvn <- function(n,mu,sig) { ## MVN random deviates
  L <- mroot(sig);m <- ncol(L);
  t(mu + L%*%matrix(rnorm(m*n),m,n)) 
}
#define number of posterior draws for plotting.----
N <- 300

#Load fits
fit <- readRDS("D:/Zeus/ETH_zurich_MSc/ETHz_S4/MatserThesis_Crowther_Lab/Results/Analysis_Results/Demo_GAMs/demographic_fits_gam_separate_FunDivEurope_FSD.rds")

# load the plot-level and tree-level FunDivEUROPE dataset
FunDivEUROPE_plots_EU_Full <- readRDS("D:/Zeus/ETH_zurich_MSc/ETHz_S4/MatserThesis_Crowther_Lab/Data/FunDivEurope/FunDivEurope_Inventory_DE_SW_ES_2022-09/FunDivEUROPE_plots_EU_Full.rds")
FunDivEUROPE_trees_EU_Full <- readRDS("D:/Zeus/ETH_zurich_MSc/ETHz_S4/MatserThesis_Crowther_Lab/Data/FunDivEurope/FunDivEurope_Inventory_DE_SW_ES_2022-09/FunDivEUROPE_trees_EU_Full.rds")

# we accounted for the potential impact of human management on the 
# establishment of monocultures by only keeping plots 
# 1) with at least ten trees, 
# 2) with at least two species, 
# 3) in which the relative basal area of the tree species with the largest cumulative basal area in the plot was smaller than 75%
d1 <- FunDivEUROPE_plots_EU_Full[FunDivEUROPE_plots_EU_Full$spp_maxprop<=0.75 & FunDivEUROPE_plots_EU_Full$spp.count>=2 & FunDivEUROPE_plots_EU_Full$stem.density>=10,]

# subset the tree-level data to have the same plot-id as the plot-level dataset
d2 <- FunDivEUROPE_trees_EU_Full[FunDivEUROPE_trees_EU_Full$plotcode%in%d1$plotcode,]

# compute deciduous basal area of each plot
d1$BASAL.de <- d1$BA-d1$BA_EV
# change names of some columns
setnames(d1,'BA','BASAL.plot')
setnames(d1,'BA_EV','BASAL.ev')
setnames(d1,'BA_DE' ,'BASAL.de')
setnames(d2,'BA_EV','BASAL.ev')
setnames(d2,'BA_DE','BASAL.de')
setnames(d2,'dbh1' ,'PREVDIA.mm')
setnames(d2,'dbh2' ,'DIA.mm')

# add plot-level properties to tree-level data
d2 <- left_join(d2,d1[,c("plotcode","stem.density","relEV","BASAL.ev", "BASAL.plot")],by="plotcode")

### Get spatial cluster of each plot. we partitioned the global forest zones using a ‘fishing net’ with 10 arc-min (~20km) grid size.
fishNet <- raster() # create an empty raster template 
res(fishNet) <- 1/6 # set the resolution of the template as 1/6 degree. 
values(fishNet) <- 1:ncell(fishNet) # allocate values to the 'fishing net' raster
d1$Cluster <- raster::extract(fishNet, d1[,c("longitude", "latitude")])
d2$Cluster <- raster::extract(fishNet, d2[,c("longitude", "latitude")])

# set plot code and cluster as factor
d1$plotcode <- as.factor(d1$plotcode)
d2$plotcode <- as.factor(d2$plotcode)
d2$Cluster <- as.factor(d2$Cluster)
d1$Cluster <- as.factor(d1$Cluster)

# remove NA
d1 <- d1%>%drop_na()
d2 <- d2%>%drop_na()

# change units from mm to cm
d2$PREVDIA.cm <- d2$PREVDIA.mm/10
d2$DIA.cm <- d2$DIA.mm/10
d2$BASAL.plot <- d2$BASAL.plot/100
d1$BASAL.plot <- d1$BASAL.plot/100
d1$BASAL.ev <- d1$BASAL.ev/100

# extract GAMs and environmental conditions from fits
env <- fit$env.cov
fit <- fit$y.feedback
# extract mortality GAMs
fit.de <- fit$M.mod.de
fit.ev <- fit$M.mod.ev

#Estimate impact of relEV on DE vs. EV mortality.
de.de.dat <- c(0,0,25000,30,22,env)
names(de.de.dat)[1:5] <- c('ev','relEV','BASAL.plot','stem.density','PREVDIA.cm')
#add cluster level random effect that will be ignored.
check1 <- fit.de$model$Cluster
check2 <- fit.ev$model$Cluster
check <- check1[check1 %in% check2]
de.de.dat <- data.frame(t(de.de.dat))
de.de.dat$Cluster <- as.factor(check[10])

#make the rest of the data sets.
de.ev.dat <- de.de.dat
de.ev.dat$relEV <- 1
ev.ev.dat <- de.ev.dat
ev.ev.dat$ev <- 1
ev.de.dat <- ev.ev.dat
ev.de.dat$relEV <- 0
de.de <- predict_gam_well(fit.de, newdata = de.de.dat, ranef.lab = "Cluster")
de.ev <- predict_gam_well(fit.de, newdata = de.ev.dat, ranef.lab = "Cluster")
ev.de <- predict_gam_well(fit.ev, newdata = ev.de.dat, ranef.lab = "Cluster")
ev.ev <- predict_gam_well(fit.ev, newdata = ev.ev.dat, ranef.lab = "Cluster")

y.de.de <- rnorm(N, de.de$fit, de.de$se.fit)
y.de.ev <- rnorm(N, de.ev$fit, de.ev$se.fit)
y.ev.de <- rnorm(N, ev.de$fit, ev.de$se.fit)
y.ev.ev <- rnorm(N, ev.ev$fit, ev.ev$se.fit)
y.de <- boot::inv.logit(c(y.de.de, y.ev.de)) #undo log link to get to scale of mortality.
y.ev <- boot::inv.logit(c(y.de.ev, y.ev.ev))

#convert to survival, rather than mortality.
y.de.surv <- 1 - y.de
y.ev.surv <- 1 - y.ev

# extract recruitment GAMs
fit.de <- fit$R.mod.de
fit.ev <- fit$R.mod.ev

#get mean and se.
ref.dat <- c(mean(d2$PREVDIA.cm, na.rm = T),mean(d1$BASAL.de, na.rm = T), mean(d1$BASAL.ev, na.rm = T), mean(d1$relEV, na.rm = T), mean(d1$BASAL.plot, na.rm = T),mean(d1$stem.density, na.rm = T))
names(ref.dat) <- c('PREVDIA.cm','BASAL.de','BASAL.ev','relEV','BASAL.plot','stem.density')
ref.dat['relEV'] <- 0
de.de.dat <- c(ref.dat,env)
check1 <- fit.de$model$Cluster
check2 <- fit.ev$model$Cluster
check <- check1[check1 %in% check2]
de.de.dat <- data.frame(t(de.de.dat))
de.de.dat$Cluster <- as.factor(check[1])

de.ev.dat <- de.de.dat
de.ev.dat$relEV <- 1
ev.ev.dat <- de.ev.dat
ev.ev.dat$ev <- 1
ev.de.dat <- ev.ev.dat
ev.de.dat$relEV <- 0
de.de <- predict_gam_well(fit.de, newdata = de.de.dat, ranef.lab = "Cluster")
de.ev <- predict_gam_well(fit.de, newdata = de.ev.dat, ranef.lab = "Cluster")
ev.de <- predict_gam_well(fit.ev, newdata = ev.de.dat, ranef.lab = "Cluster")
ev.ev <- predict_gam_well(fit.ev, newdata = ev.ev.dat, ranef.lab = "Cluster")

N <- 300
y.de.de <- rnorm(N, de.de$fit, de.de$se.fit)
y.de.ev <- rnorm(N, de.ev$fit, de.ev$se.fit)
y.ev.de <- rnorm(N, ev.de$fit, ev.de$se.fit)
y.ev.ev <- rnorm(N, ev.ev$fit, ev.ev$se.fit)

y.de.recr <- exp(c(y.de.de, y.ev.de)) #undo log link to get to scale of recruitment.
y.ev.recr <- exp(c(y.de.ev, y.ev.ev))

#load, workup growth data.----
fit.de <- fit$G.mod.de
fit.ev <- fit$G.mod.ev

#get mean and se.
ref.dat <- c(mean(d2$PREVDIA.cm, na.rm = T),mean(d1$BASAL.de, na.rm = T), mean(d1$BASAL.ev, na.rm = T), mean(d1$relEV, na.rm = T), mean(d1$BASAL.plot, na.rm = T),mean(d1$stem.density, na.rm = T))
names(ref.dat) <- c('PREVDIA.cm','BASAL.de','BASAL.ev','relEV','BASAL.plot','stem.density')
ref.dat['relEV'] <- 0
de.de.dat <- c(ref.dat,env)
check1 <- fit.de$model$Cluster
check2 <- fit.ev$model$Cluster
check <- check1[check1 %in% check2]
de.de.dat <- data.frame(t(de.de.dat))
de.de.dat$Cluster <- as.factor(check[1])
de.ev.dat <- de.de.dat
de.ev.dat$relEV <- 1
ev.ev.dat <- de.ev.dat
ev.ev.dat$ev <- 1
ev.de.dat <- ev.ev.dat
ev.de.dat$relEV <- 0
de.de <- predict(fit.de, newdata = de.de.dat, se.fit = T, exclude = c("s(Cluster)"), newdata.guaranteed = T)
de.ev <- predict(fit.de, newdata = de.ev.dat, se.fit = T, exclude = c("s(Cluster)"), newdata.guaranteed = T)
ev.de <- predict(fit.ev, newdata = ev.de.dat, se.fit = T, exclude = c("s(Cluster)"), newdata.guaranteed = T)
ev.ev <- predict(fit.ev, newdata = ev.ev.dat, se.fit = T, exclude = c("s(Cluster)"), newdata.guaranteed = T)
N <- 300
y.de.de <- rnorm(N, de.de$fit, de.de$se.fit)
y.de.ev <- rnorm(N, de.ev$fit, de.ev$se.fit)
y.ev.de <- rnorm(N, ev.de$fit, ev.de$se.fit)
y.ev.ev <- rnorm(N, ev.ev$fit, ev.ev$se.fit)
y.de.grow <- (c(y.de.de, y.ev.de))
y.ev.grow <- (c(y.de.ev, y.ev.ev))

#subtract initial diameter to get increment.
y.de.grow <- y.de.grow - ref.dat['PREVDIA.cm']
y.ev.grow <- y.ev.grow - ref.dat['PREVDIA.cm']
y.de.grow[y.de.grow<0] <- 0.01
y.ev.grow[y.ev.grow<0] <- 0.01



###---ggplot for demogrophic comparison----####
## need to run growth ananlysis in the supp._Fig._3 code
demoVisualization_df <- as.data.frame(cbind(values=c(y.de.surv, y.ev.surv, y.de.grow, 
                                                     y.ev.grow, y.de.recr, y.ev.recr),
                                            forests = rep(rep(c("DE Forest", "EV Forest"), 3), each=600), 
                                            trees = rep(rep(c("All DE trees", "All EV trees"), 6), each=300),
                                            demo = rep(c("survival", "growth", "recruitment"), each=1200)
))
demoVisualization_df$values <- as.numeric(demoVisualization_df$values)

# surival
ggplot(demoVisualization_df[demoVisualization_df$demo=="survival",], aes(x = forests, y = values, color = trees, fill=trees)) +  # ggplot function
  geom_boxplot(alpha=0.2, width=0.6)+
  scale_color_manual(values=c("red", "blue"))+
  scale_fill_manual(values=c("red", "blue"))+
  ggtitle("(A)")+
  ylim(c(0.7,1))+
  theme_classic()+
  labs(
    y = "Survival Probability")+
  # scale_x_log10()+
  theme(title = element_text(size=14,face="bold"),
        axis.text.x=element_text(size=12, face="bold"),
        axis.text.y=element_text(size=12, face="bold"),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=14,face="bold"),
        legend.text = element_text(colour = "black", size = 10, face = "bold"),
        legend.title = element_blank())

# growth
ggplot(demoVisualization_df[demoVisualization_df$demo=="growth",], aes(x = forests, y = values/10, color = trees, fill=trees)) +  # ggplot function
  geom_boxplot(alpha=0.2, width=0.6)+
  scale_color_manual(values=c("red", "blue"))+
  scale_fill_manual(values=c("red", "blue"))+
  ggtitle("(B)")+
  # ylim(c(0.7,1))+
  theme_classic()+
  labs(
    y = "Diameter Increment (cm/yr)")+
  # scale_x_log10()+
  theme(title = element_text(size=14,face="bold"),
        axis.text.x=element_text(size=12, face="bold"),
        axis.text.y=element_text(size=12, face="bold"),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=14,face="bold"),
        legend.text = element_text(colour = "black", size = 10, face = "bold"),
        legend.title = element_blank())

# recruitment
ggplot(demoVisualization_df[demoVisualization_df$demo=="recruitment",], aes(x = forests, y = values, color = trees, fill=trees)) +  # ggplot function
  geom_boxplot(alpha=0.2, width=0.6)+
  scale_color_manual(values=c("red", "blue"))+
  scale_fill_manual(values=c("red", "blue"))+
  ggtitle("(C)")+
  # ylim(c(0.7,1))+
  theme_classic()+
  labs(
    y = 'Recruitment per year')+
  # scale_x_log10()+
  theme(title = element_text(size=14,face="bold"),
        axis.text.x=element_text(size=12, face="bold"),
        axis.text.y=element_text(size=12, face="bold"),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=14,face="bold"),
        legend.text = element_text(colour = "black", size = 10, face = "bold"),
        legend.title = element_blank())

### get percentage demographic differences between the two groups
surv_DE <- demoVisualization_df$values[1:300]/demoVisualization_df$values[301:600]-1
mean(surv_DE)
sd(surv_DE)

surv_EV <- demoVisualization_df$values[901:1200]/demoVisualization_df$values[601:900]-1
mean(surv_EV)
sd(surv_EV)

grow_DE <- demoVisualization_df$values[1201:1500]/demoVisualization_df$values[1501:1800]-1
mean(grow_DE)
sd(grow_DE)

grow_EV <- demoVisualization_df$values[2101:2400]/demoVisualization_df$values[1801:2100]-1
mean(grow_EV)
sd(grow_EV)

