#paths.r
#Main file dump directories.----
host <- system('hostname', intern = T)
storage.dir <- '**/Continental_Analysis/FIA_US/'

#FIA input paths.----
FIA9.dir.path <- paste0(storage.dir,'Data/FIA_DataMart/FIADB_USA_2020_08_24_SQLite/')

FIAdb.path <- paste0(FIA7.dir.path,'FIADB_USA.db')

#Other raster data outside of main storage directory.----
EPA_L2_ecoregions_raster.path <- paste0(storage.dir,'Data/EPA_US/na_cec_eco_l2/NA_CEC_Eco_Level2.shp')


#FIA filtered output paths.----
fia.dir <- paste0(storage.dir,'FIA_Output/')

#All FIA data broken up by remeasurement.
all.present.path <- paste0(fia.dir,'FIA.all.present.rds')
all.past1.path   <- paste0(fia.dir,'FIA.all.past1.rds')
all.past2.path   <- paste0(fia.dir,'FIA.all.past2.rds')
all.past3.path   <- paste0(fia.dir,'FIA.all.past3.rds')

#paths for grabbing and returning lab composite environmental data.----
data_for_composite.path <- paste0(storage.dir, 'data_for_composite_FIA.csv')

#FIA formatted analysis products.----
Product_1.path         <- paste0(fia.dir,"Product_1.rds")
Product_2.path         <- paste0(fia.dir,"Product_2.rds")

Product_1.path         <- paste0(fia.dir,"Product_1.subset.rds")
Product_2.subset.path  <- paste0(fia.dir,"Product_2.subset.rds")

#GAM and simulation model output paths.----
model.dir <- paste0(storage.dir,'Model_Output/')
cmd <- paste0('mkdir -p ',model.dir)
system(cmd)

#GAM fits.
demographic_fits_gam_separate.path <- paste0(model.dir,'demographic_fits_gam_separate.rds')
demographic_fits_gam_ecoregion.path <- paste0(model.dir,'demographic_fits_gam_ecoregion.rds')
demographic_fits_gam_species.path <- paste0(model.dir,'demographic_fits_gam_species.rds')

#Demographic simulation output paths.
null_vs_feedback_simulation_output.path <- paste0(model.dir,'null_vs_feedback_simulation_output.rds')
null_vs_feedback_simulation_output_RE_uniform.path <- paste0(model.dir,'null_vs_feedback_simulation_output_RE_uniform.rds')

#hysteresis simulations.
initial_condition_hysteresis_simulation.path <- paste0(model.dir,'initial_condition_hysteresis_simulation.rds')


#spatial variogram output.----
variogram_data.path <- paste0(model.dir,'variogram_data.rds')
