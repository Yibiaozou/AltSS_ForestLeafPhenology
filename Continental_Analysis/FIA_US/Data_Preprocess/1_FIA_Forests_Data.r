#Querying the FIA database.
#This code is adapted from code used in Averill et al. (2022)
#Yibiao downloaded the most recent PLOT, COND, TREE_GRM_EST, and TREE tables from the FIA on 24 Aug, 2020
#This uses RSQlite to query the respective tables, rather than pSQL as done previously by R. Kelly.
#This script takes a while to run (>20 minutes).
#clear environment, load packages.
rm(list=ls())
library(data.table)
library(RSQLite)
# library(doParallel)
library(DBI)
# library(foreach)
storage.dir <- '**/Continental_Analysis/FIA_US/'
setwd(storage.dir)
source('Project_Functions/paths.r')
source('Project_Functions/tic_toc.r')

#Connect to FIA9 database.----
con <- dbConnect(SQLite(), dbname = FIAdb.path)

#File paths to reference data.----
file.pft = "Data/Other_Data_Product/gcbPFT.csv"
file.myc = "Data/Other_Data_Product/mycorrhizal_SPCD_data.csv"
file.gym = "Data/Other_Data_Product/gymnosperm_family_genera.csv"

#state codes.----
states <- read.csv('Data/Other_Data_Product/FIA_state_codes_regions.csv')
#Drop Hawaii and Alaska, thus to focus on the mainland US. Still 49 'states' because DC counts.
states <- states[!(states$STATECD %in% c(2,15)),]$STATECD
states <- paste(states, collapse = ', ')
states <- paste0('(',states,')')

###---Query PLOT table.----
cat("Query PLOT.../n")
tic()
query <- paste("select 
                CN, STATECD, PREV_PLT_CN, REMPER, LAT, LON, ELEV, COUNTYCD 
                from 
                PLOT 
                where 
                STATECD IN",states)
PLOT = dbGetQuery(con, query)
PLOT = data.table(PLOT)
toc()
setnames(PLOT,"CN","PLT_CN")
states  = sort(unique(PLOT$STATECD))
n.state = length(states)

# Remove this one miscellaneous plot, per Trevor Andrews
PLOT = PLOT[ PLT_CN!= 134680578010854 ]


###---Query COND table.----
cat("Query COND.../n")
tic()
query <- paste("select 
                        PLT_CN, CONDID, STDORGCD, STDAGE, CONDPROP_UNADJ, AFFORESTATION_CD
                        from 
                        COND")
COND = dbGetQuery(con, query)
COND = data.table(COND)
toc()

###---Query SUBP_COND table.----
cat("Query SUBP_COND.../n")
tic()
SUBP_COND = dbGetQuery(con, "select 
                       CN, PLT_CN, CONDID, SUBP, SUBPCOND_PROP
                       from 
                       SUBP_COND")
SUBP_COND = data.table(SUBP_COND)
toc()

#Subset to subplots in forested condition (CONDID == 1)
forest_prop <- SUBP_COND[CONDID == 1,]
#Calculate forested condition proportion.
forest_prop <- data.table(aggregate(SUBPCOND_PROP ~ PLT_CN, data = forest_prop, FUN = 'sum'))
forest_prop[,SUBPCOND_PROP := SUBPCOND_PROP / 4]
colnames(forest_prop)[2] <- 'forest_proportion'
#merge into COND table
COND <- merge(COND,forest_prop, all.x = T)

#need to have at least one plot with condition =1 (forested)
COND <- COND[CONDID == 1,]

#only keep plots that are 100% forested.
COND <- COND[COND$forest_proportion == 1,]

###---merge PLOT and COND tables.----
PC = merge(COND, PLOT, by="PLT_CN")
PC <- data.table(PC)

#remove PLOT, COND and SUBP_COND tables from memory.
rm(PLOT,COND, SUBP_COND)

cl<-makeCluster(8)
registerDoParallel(cl)
clusterEvalQ(cl , c(library(data.table),library(foreach), library(RSQLite)))
###---Query TREE table.----
cat("Query TREE.../n") 
tic()

TREE <- list()
for(i in 1:length(states)){
  query = paste('select CN, PREV_TRE_CN, PLT_CN, INVYR, CONDID, DIA, TPA_UNADJ, CARBON_AG, CARBON_BG,
              SPCD, STOCKING, STATUSCD, PREVDIA, PREV_STATUS_CD, P2A_GRM_FLG, RECONCILECD, AGENTCD, TPAMORT_UNADJ,
              DIAHTCD, HT, HTCD, ACTUALHT, CCLCD
                from TREE 
                WHERE (PREVDIA>5 OR DIA>5) AND (STATUSCD=1 OR PREV_STATUS_CD=1) AND 
                STATECD IN (', paste(states[i],collapse=','), ')')
  pre.tree = as.data.table(dbGetQuery(con, query))
  TREE[[i]] <- pre.tree
}

toc()

#collapse dataframe.
TREE <- do.call(rbind, TREE)

###--- Filter TREE table.----
cat("Filter TREE .../n")
# By plot/cond criteria
TREE = TREE[PLT_CN %in% PC$PLT_CN,]

# CONDID ("Remove edge effects" --TA)
#Colin- this is calculated but never used to filter. Can probably drop.
TREE[, CONmax := max(CONDID, na.rm = T), by=PLT_CN]

# STATUSCD
#Get largest statuscd observe within a plot.
# *** RK: Next line looks wrong. It's a sum, not max, despite the name. I did rewrite the line but this is equivalent to what Travis had so keeping for now.
TREE[, STATUSCDmax := sum(3*as.integer(STATUSCD==3), na.rm = T), by=PLT_CN]

# RECONCILECD. This is just signaling that a tree isn't a new tree to the plot.
TREE[is.na(RECONCILECD), RECONCILECD :=0] # Set NA values to 0.

# Filter
#STATUSCD=0 are remeasured trees that shouldn't be there.
#STATUSCD=3 means a tree was cut down by humans.
#RECONCILECD<=4 means these are acceptable reasons to have missed counting trees in the previous inventory.
TREE = TREE[STATUSCDmax!=3 & STATUSCD!=0 & RECONCILECD<=4 ]

#Count the number of trees in a plot. 
TREE[,n.trees  := length(TPA_UNADJ), by=PLT_CN]

#Gotta keep CN numbers straight between trees and plots.
setnames(TREE,'CN','TRE_CN')

###--- Merge in PFTs and mycorrhizal associations.----
cat("Merge in PFTs and mycorrhizal associations.../n")
MCDPFT     = as.data.table(read.csv(file.pft, header = T)) #load in PFT assignments. Mike Dietze (MCD) made this file long ago.
CA_myctype = as.data.table(read.csv(file.myc, header = T)) #load in mycorrhizal associations.
CA_myctype = CA_myctype[,c("SPCD","MYCO_ASSO","GENUS","SPECIES"),with=F]
TREE = merge(TREE, MCDPFT    , all.x=T, by = "SPCD")
TREE = merge(TREE, CA_myctype, all.x=T, by = "SPCD")

#Assign whether or not a tree is a gymnosperm.----
gymnosperm <- read.csv(file.gym)
TREE$gymno <- ifelse(TREE$GENUS %in% gymnosperm$genus, 1, 0)

#Link together temporal sequences of trees and PC plot codes.----
PC.present   <-   PC[!(PLT_CN %in% PREV_PLT_CN),] #newest observations are not a PREV_PLT_CN of anything.
PC.past1     <-   PC[PLT_CN %in% PC.present$PREV_PLT_CN,]
PC.past2     <-   PC[PLT_CN %in% PC.past1$PREV_PLT_CN,]
PC.past3     <-   PC[PLT_CN %in% PC.past2$PREV_PLT_CN,]
TREE.present <- TREE[PLT_CN %in% PC.present$PLT_CN]
TREE.past1   <- TREE[PLT_CN %in% PC.past1$PLT_CN]
TREE.past2   <- TREE[PLT_CN %in% PC.past2$PLT_CN]
TREE.past3   <- TREE[PLT_CN %in% PC.past3$PLT_CN]
#merge PC and TREE tables.
all.present <- merge(TREE.present, PC.present, by = 'PLT_CN')
all.past1   <- merge(TREE.past1  , PC.past1  , by = 'PLT_CN')
all.past2   <- merge(TREE.past2  , PC.past2  , by = 'PLT_CN')
all.past3   <- merge(TREE.past3  , PC.past3  , by = 'PLT_CN')
#grab data for composite.
for.composite <- all.present[,.(PLT_CN,LAT,LON)]
for.composite <- for.composite[!duplicated(for.composite),]
for.composite$PLT_CN <- as.numeric(gsub('"', "",for.composite$PLT_CN))
colnames(for.composite) <- c('PLT_CN','latitude','longitude')

#save outputs.----
cat("Saving output.../n")
tic()
saveRDS(all.present,file = all.present.path)
saveRDS(all.past1  ,file = all.past1.path  )
saveRDS(all.past2  ,file = all.past2.path  )
saveRDS(all.past3  ,file = all.past3.path  )
write.csv(for.composite, data_for_composite.path)
toc()
###end script.
cat("Script complete./n")
