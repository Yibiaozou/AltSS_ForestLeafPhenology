#Fitting gams for separate DE-EV models and PC scores.

rm(list = ls())

storage.dir <- '**/Continental_Analysis/FIA_US/'
setwd(storage.dir)
source('Project_Functions/paths.r')
source('Project_Functions/tic_toc.r')

#set output path.----

memory.limit(5e5)

output.path <- demographic_fits_gam_separate.path

d1 <- readRDS(Product_1.path)
d2 <- readRDS(Product_2.path)

d1$plot_id <- 1:nrow(d1)
# FIA_US_Coordinates <- d1[,c("plot_id","latitude", "longitude")]
# FIA_US_Coordinates <- FIA_US_Coordinates %>% drop_na()

# For FIA we already remove managed plot, therefore here we only remove small plots with less than 10 trees
d1 <- d1[d1$count>=10,]

#Subset and rename some things.-----
d1 <- d1[d1$PLT_CN %in% d2$PLT_CN,]
setnames(d1,'plot.BASAL','BASAL.plot')
setnames(d2,'plot.BASAL','BASAL.plot')
setnames(d1,'BASAL.Evergreen','BASAL.ev')
setnames(d1,'BASAL.Deciduous' ,'BASAL.de')


### Get spatial cluster of each plot. we partitioned the global forest zones using a ‘fishing net’ with 10 arc-min (~20km) grid size.
fishNet <- raster() # create an empty raster template 
res(fishNet) <- 1/6 # set the resolution of the template as 1/6 degree. 
values(fishNet) <- 1:ncell(fishNet) # allocate values to the 'fishing net' raster
# extract spatial cluster for each plot
d1$Cluster <- raster::extract(fishNet, d1[,c("longitude", "latitude")])
d2$Cluster <- raster::extract(fishNet, d2[,c("longitude", "latitude")])

#Set global k.----
kk <- 5    #global max number of 'knots' to use in smoothing functions.
bs <- 'tp' #thin plate regression splines are default as well, but you can change things here if you like.
nt <- 8    #number of threads.

#Fit growth, recruitment and mortality models.----
#Environmental models without feedbacks.
cat('Fitting null models...\n');tic()
# evergreen recruitment
R.mod.ev <- bam(recruit.ev ~        s(BASAL.ev, k=kk, bs=bs) + s(mat, k=kk, bs=bs)  + s(BASAL.plot, k=kk, bs=bs) + s(stem.density, k=kk, bs=bs) +
                  s(Cluster, bs = 're') +
                  s(PC1, k=kk, bs=bs) + s(PC2, k=kk, bs=bs) + s(PC3, k=kk, bs=bs) + s(PC4, k=kk, bs=bs) + s(PC5, k=kk, bs=bs) + s(PC6, k=kk, bs=bs) + s(PC7, k=kk, bs=bs) + s(PC8, k=kk, bs=bs) + s(PC9, k=kk, bs=bs) + s(PC10, k =kk, bs=bs), 
                data = d1, family = 'poisson', nthreads=nt, discrete=T) 
toc()
tic()
# deciduous recruitment
R.mod.de <- bam(recruit.de ~        s(BASAL.de, k=kk, bs=bs) + s(mat, k=kk, bs=bs)  + s(BASAL.plot, k=kk, bs=bs) + s(stem.density, k=kk, bs=bs) +
                  s(Cluster, bs = 're') +
                  s(PC1, k=kk, bs=bs) + s(PC2, k=kk, bs=bs) + s(PC3, k=kk, bs=bs) + s(PC4, k=kk, bs=bs) + s(PC5, k=kk, bs=bs) + s(PC6, k=kk, bs=bs) + s(PC7, k=kk, bs=bs) + s(PC8, k=kk, bs=bs) + s(PC9, k=kk, bs=bs) + s(PC10, k =kk, bs=bs), 
                data = d1, family = 'poisson', nthreads=nt, discrete=T) 
toc()
tic()
# evergreen mortality
M.mod.ev <- bam(mortality ~          s(Soil_PH, k=kk, bs=bs) + s(mat, k=kk, bs=bs) + s(BASAL.plot, k=kk, bs=bs) + s(stem.density, k=kk, bs=bs) + s(PREVDIA.cm, k=kk, bs=bs) + 
                  s(Cluster, bs = 're') +
                  s(PC1, k=kk, bs=bs) + s(PC2, k=kk, bs=bs) + s(PC3, k=kk, bs=bs) + s(PC4, k=kk, bs=bs) + s(PC5, k=kk, bs=bs) + s(PC6, k=kk, bs=bs) + s(PC7, k=kk, bs=bs) + s(PC8, k=kk, bs=bs) + s(PC9, k=kk, bs=bs) + s(PC10, k =kk, bs=bs), 
                data = d2[d2$ev == 1,], family = 'binomial', nthreads=nt, discrete=T) 
toc()

tic()
# deciduous mortality
M.mod.de <- bam(mortality ~          s(Soil_PH, k=kk, bs=bs) +s(mat, k=kk, bs=bs) +  s(BASAL.plot, k=kk, bs=bs) + s(stem.density, k=kk, bs=bs) + s(PREVDIA.cm, k=kk, bs=bs) + 
                  s(Cluster, bs = 're') +
                  s(PC1, k=kk, bs=bs) + s(PC2, k=kk, bs=bs) + s(PC3, k=kk, bs=bs) + s(PC4, k=kk, bs=bs) + s(PC5, k=kk, bs=bs) + s(PC6, k=kk, bs=bs) + s(PC7, k=kk, bs=bs) + s(PC8, k=kk, bs=bs) + s(PC9, k=kk, bs=bs) + s(PC10, k =kk, bs=bs), 
                data = d2[d2$ev == 0,], family = 'binomial', nthreads=nt, discrete=T) 
toc()
tic()
# evergreen growth
G.mod.ev <- bam(DIA.cm   ~            s(Soil_PH, k=kk, bs=bs) + s(mat, k=kk, bs=bs) + s(BASAL.plot, k=kk, bs=bs) + s(stem.density, k=kk, bs=bs) + s(PREVDIA.cm, k=kk, bs=bs) + 
                  s(Cluster, bs = 're') +
                  s(PC1, k=kk, bs=bs) + s(PC2, k=kk, bs=bs) + s(PC3, k=kk, bs=bs) + s(PC4, k=kk, bs=bs) + s(PC5, k=kk, bs=bs) + s(PC6, k=kk, bs=bs) + s(PC7, k=kk, bs=bs) + s(PC8, k=kk, bs=bs) + s(PC9, k=kk, bs=bs) + s(PC10, k =kk, bs=bs), 
                data = d2[d2$ev == 1,], nthreads=nt, discrete=T)
toc()
tic()
# deciduous growth
G.mod.de <- bam(DIA.cm   ~            s(Soil_PH, k=kk, bs=bs) + s(mat, k=kk, bs=bs) + s(BASAL.plot, k=kk, bs=bs) + s(stem.density, k=kk, bs=bs) + s(PREVDIA.cm, k=kk, bs=bs) + 
                  s(Cluster, bs = 're') +
                  s(PC1, k=kk, bs=bs) + s(PC2, k=kk, bs=bs) + s(PC3, k=kk, bs=bs) + s(PC4, k=kk, bs=bs) + s(PC5, k=kk, bs=bs) + s(PC6, k=kk, bs=bs) + s(PC7, k=kk, bs=bs) + s(PC8, k=kk, bs=bs) + s(PC9, k=kk, bs=bs) + s(PC10, k =kk, bs=bs), 
                data = d2[d2$ev == 0,], nthreads=nt, discrete=T)
cat('Null models fit.\n');toc()
#wrap output and name.
n.feedback <- list(G.mod.de, G.mod.ev, M.mod.de, M.mod.ev, R.mod.de, R.mod.ev)
names(n.feedback) <- c('G.mod.de','G.mod.ev','M.mod.de','M.mod.ev','R.mod.de','R.mod.ev')

#Environmental models with feedbacks.
#gam models.
cat('Fitting feedback models...\n');tic()
# evergreen recruitment
R.mod.ev <- bam(recruit.ev ~ s(relEV, k=kk, bs=bs) + s(BASAL.ev, k=kk, bs=bs) + s(mat, k=kk, bs=bs) + s(Soil_PH, k=kk, bs=bs) + s(BASAL.plot, k=kk, bs=bs) + s(stem.density, k=kk, bs=bs) +
                  s(county.ID, bs = 're') +
                  s(PC1, k=kk, bs=bs) + s(PC2, k=kk, bs=bs) + s(PC3, k=kk, bs=bs) + s(PC4, k=kk, bs=bs) + s(PC5, k=kk, bs=bs) + s(PC6, k=kk, bs=bs) + s(PC7, k=kk, bs=bs) + s(PC8, k=kk, bs=bs) + s(PC9, k=kk, bs=bs) + s(PC10, k =kk, bs=bs), 
                data = d1, family = 'poisson', nthreads=nt, discrete=T)
toc()
tic()
# deciduous recruitment
R.mod.de <- bam(recruit.de ~ s(relEV, k=kk, bs=bs) + s(BASAL.de, k=kk, bs=bs) + s(mat, k=kk, bs=bs) + s(Soil_PH, k=kk, bs=bs) + s(BASAL.plot, k=kk, bs=bs) + s(stem.density, k=kk, bs=bs) +
                  s(county.ID, bs = 're') +
                  s(PC1, k=kk, bs=bs) + s(PC2, k=kk, bs=bs) + s(PC3, k=kk, bs=bs) + s(PC4, k=kk, bs=bs) + s(PC5, k=kk, bs=bs) + s(PC6, k=kk, bs=bs) + s(PC7, k=kk, bs=bs) + s(PC8, k=kk, bs=bs) + s(PC9, k=kk, bs=bs) + s(PC10, k =kk, bs=bs), 
                data = d1, family = 'poisson', nthreads=nt, discrete=T)
toc()
tic()
# evergreen mortality
M.mod.ev <- bam(mortality ~ s(relEV, k=kk, bs=bs) + s(mat, k=kk, bs=bs) + s(Soil_PH, k=kk, bs=bs) + s(BASAL.plot, k=kk, bs=bs) + s(stem.density, k=kk, bs=bs) + s(PREVDIA.cm, k=kk, bs=bs) + 
                  s(county.ID, bs = 're') +
                  s(PC1, k=kk, bs=bs) + s(PC2, k=kk, bs=bs) + s(PC3, k=kk, bs=bs) + s(PC4, k=kk, bs=bs) + s(PC5, k=kk, bs=bs) + s(PC6, k=kk, bs=bs) + s(PC7, k=kk, bs=bs) + s(PC8, k=kk, bs=bs) + s(PC9, k=kk, bs=bs) + s(PC10, k =kk, bs=bs), 
                data = d2[d2$ev == 1,], family = 'binomial', nthreads=nt, discrete=T)
toc()
tic()
# deciduous mortality
M.mod.de <- bam(mortality ~ s(relEV, k=kk, bs=bs) + s(mat, k=kk, bs=bs) + s(Soil_PH, k=kk, bs=bs) + s(BASAL.plot, k=kk, bs=bs) + s(stem.density, k=kk, bs=bs) + s(PREVDIA.cm, k=kk, bs=bs) + 
                  s(county.ID, bs = 're') +
                  s(PC1, k=kk, bs=bs) + s(PC2, k=kk, bs=bs) + s(PC3, k=kk, bs=bs) + s(PC4, k=kk, bs=bs) + s(PC5, k=kk, bs=bs) + s(PC6, k=kk, bs=bs) + s(PC7, k=kk, bs=bs) + s(PC8, k=kk, bs=bs) + s(PC9, k=kk, bs=bs) + s(PC10, k =kk, bs=bs), 
                data = d2[d2$ev == 0,], family = 'binomial', nthreads=nt, discrete=T)
toc()
tic()
# evergreen growth
G.mod.ev <- bam(DIA.cm   ~ s(relEV, k=kk, bs=bs) + s(mat, k=kk, bs=bs) + s(Soil_PH, k=kk, bs=bs) + s(BASAL.plot, k=kk, bs=bs) + s(stem.density, k=kk, bs=bs) + s(PREVDIA.cm, k=kk, bs=bs) + 
                  s(county.ID, bs = 're') +
                  s(PC1, k=kk, bs=bs) + s(PC2, k=kk, bs=bs) + s(PC3, k=kk, bs=bs) + s(PC4, k=kk, bs=bs) + s(PC5, k=kk, bs=bs) + s(PC6, k=kk, bs=bs) + s(PC7, k=kk, bs=bs) + s(PC8, k=kk, bs=bs) + s(PC9, k=kk, bs=bs) + s(PC10, k =kk, bs=bs), 
                data = d2[d2$ev == 1,], nthreads=nt, discrete=T)
toc()
tic()
# deciduous growth
G.mod.de <- bam(DIA.cm   ~ s(relEV, k=kk, bs=bs) + s(mat, k=kk, bs=bs) + s(Soil_PH, k=kk, bs=bs) + s(BASAL.plot, k=kk, bs=bs) + s(stem.density, k=kk, bs=bs) + s(PREVDIA.cm, k=kk, bs=bs) + 
                  s(county.ID, bs = 're') +
                  s(PC1, k=kk, bs=bs) + s(PC2, k=kk, bs=bs) + s(PC3, k=kk, bs=bs) + s(PC4, k=kk, bs=bs) + s(PC5, k=kk, bs=bs) + s(PC6, k=kk, bs=bs) + s(PC7, k=kk, bs=bs) + s(PC8, k=kk, bs=bs) + s(PC9, k=kk, bs=bs) + s(PC10, k =kk, bs=bs), 
                data = d2[d2$ev == 0,], nthreads=nt, discrete=T)
cat('Feedback models fit.\n');toc()

#wrap output and name.
y.feedback <- list(G.mod.de, G.mod.ev, M.mod.de, M.mod.ev, R.mod.de, R.mod.ev)
names(y.feedback) <- c('G.mod.de','G.mod.ev','M.mod.de','M.mod.ev','R.mod.de','R.mod.ev')

#Grab plot environmental covariates for reference.----
#all plot-level environmental covariates. we will sample this later when simulating forests.
all.cov <- d1[,c('mat','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10')]

#mean plot-level environmental covariates.
cov <- colMeans(all.cov, na.rm=T)
names(cov) <- colnames(all.cov)

#Save models and size categories.----
cat('Saving output...\n');tic()
output <- list(n.feedback, y.feedback, cov, all.cov)
names(output) <- c('n.feedback', 'y.feedback','env.cov','all.cov')
saveRDS(output, output.path, version = 2)
cat('Output saved, script complete.\n');toc()
