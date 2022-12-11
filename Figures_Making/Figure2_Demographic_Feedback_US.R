#Plot growth, recruitment and survival of all species together as DE vs. EV.
#Convert mortality to survival.
rm(list=ls())
storage.dir <- '**/Continental_Analysis/FIA_US/'
setwd(storage.dir)
source('Project_Functions/paths.r')
source('Project_Functions/predict_gam_well.r')
library(mgcv)
library(data.table)
library(tidyverse)

#function to draw randomly from multivariate normal distribution.-----
rmvn <- function(n,mu,sig) { ## MVN random deviates
  L <- mroot(sig);m <- ncol(L);
  t(mu + L%*%matrix(rnorm(m*n),m,n)) 
}
#define number of posterior draws for plotting.----
N <- 300

#load, workup mortality (survival) data.----
#Load fits and data.
fit <- readRDS(demographic_fits_gam_separate.path)

env <- fit$env.cov
fit.de <- fit$y.feedback$M.mod.de
fit.ev <- fit$y.feedback$M.mod.ev

#Estimate impact of relEV on DE vs. EV mortality.
de.de.dat <- c(0,0,25000,30,22,env)
names(de.de.dat)[1:5] <- c('ev','relEV','BASAL.plot','stem.density','PREVDIA.cm')
#add county level randome ffect that will be ignored.
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

#load, workup recruitment data.----
#Load fits and data.
env <- fit$env.cov
fit <- fit$y.feedback
fit.de <- fit$R.mod.de
fit.ev <- fit$R.mod.ev

d1 <- readRDS(Product_1.path)
d2 <- readRDS(Product_2.path)

# For FIA we already remove managed plot, therefore here we only remove small plots with less than 10 trees
d1 <- d1[d1$count>=10,]
d1 <- d1[d1$PLT_CN %in% d2$PLT_CN,]

#get mean and se.
ref.dat <- c(mean(d2$PREVDIA.cm, na.rm = T),mean(d1$BASAL.Deciduous), mean(d1$BASAL.Evergreen), mean(d1$relEV), mean(d1$plot.BASAL),mean(d1$stem.density))
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
#Note- none of this is plotted.
fit.de <- fit$G.mod.de
fit.ev <- fit$G.mod.ev


#get mean and se.
ref.dat <- c(mean(d2$PREVDIA.cm, na.rm = T),mean(d1$BASAL.Deciduous, na.rm = T), mean(d1$BASAL.Evergreen, na.rm = T), mean(d1$relEV, na.rm = T), mean(d1$plot.BASAL, na.rm = T),mean(d1$stem.density, na.rm = T))
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


###---ggplot for demogrophic comparison----####
## need to run growth ananlysis in the supp._Fig._3 code
demoVisualization_df <- as.data.frame(cbind(values=c(y.de.surv, y.ev.surv, y.de.grow, 
                                        y.ev.grow, y.de.recr, y.ev.recr),
                               forests = rep(rep(c("DE Forest", "EV Forest"), 3), each=600), 
                               trees = rep(rep(c("All DE trees", "All EV trees"), 6), each=300),
                               demo = rep(c("survival", "growth", "recruitment"), each=1200)
                               ))
demoVisualization_df$values <- as.numeric(demoVisualization_df$values)
saveRDS(demoVisualization_df, file="D:/Zeus/ETH_zurich_MSc/ETHz_S4/MatserThesis_Crowther_Lab/Results/Analysis_Results/demoVisualization_df_US_FSD.rds")

demoVisualization_df <- readRDS("D:/Zeus/ETH_zurich_MSc/ETHz_S4/MatserThesis_Crowther_Lab/Results/Analysis_Results/demoVisualization_df_US_FSD.rds")

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

# Figure 2A
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

# Figure 2B
ggplot(demoVisualization_df[demoVisualization_df$demo=="growth",], aes(x = forests, y = values/5, color = trees, fill=trees)) +  # ggplot function
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

# Figure 2C
ggplot(demoVisualization_df[demoVisualization_df$demo=="recruitment",], aes(x = forests, y = values, color = trees, fill=trees)) +  # ggplot function
  geom_boxplot(alpha=0.2, width=0.6)+
  scale_color_manual(values=c("red", "blue"))+
  scale_fill_manual(values=c("red", "blue"))+
  ggtitle("(C)")+
  # ylim(c(0.7,1))+
  theme_classic()+
  labs(
    y = expression(paste('Recruitment per five years')))+
  # scale_x_log10()+
  theme(title = element_text(size=14,face="bold"),
        axis.text.x=element_text(size=12, face="bold"),
        axis.text.y=element_text(size=12, face="bold"),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=14,face="bold"),
        legend.text = element_text(colour = "black", size = 10, face = "bold"),
        legend.title = element_blank())

