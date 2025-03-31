
#####################################
### niche modeling #############
######################################

## install and load the required packages
library(biomod2)
library(ggplot2)
library(gridExtra)
library(raster)
library(rasterVis)
library(sf)
library(ncdf4)
library(spThin)
library(usdm)
library(caret)
library(rJava)
library(dismo)

## read data (bolivica.txt or camura.txt)
ProLau_occ <- read.table('bolivica.txt', header = T, sep = ",")
summary(ProLau_occ)

#thinning data points
occ_thin <- thin(ProLau_occ, verbose=T, 
               long.col= "longitude", 
               lat.col="latitude", 
               spec.col= "A_bolivica",  #or A_camura
               thin.par= 3, 
               reps=1,
               locs.thinned.list.return=T,
               write.files=F)

occ_thin1<-data.frame(occ_thin)
occ_thin<-data.frame(occ_thin)
occ_thin$sp1<-1

#creating a bbox
points_sf <- st_as_sf(occ_thin, coords = c("Longitude", "Latitude"), crs = 4326) # WGS84
points_utm <- st_transform(points_sf, crs = 32723) 
buffer_utm <- st_buffer(points_utm, dist = 400000)  # create a calibration area with 400 km buffer around occurrence points
buffer_wgs84 <- st_transform(buffer_utm, crs = 4326)
plot(buffer_wgs84, col = alpha("blue", 0.5))
calibration <- st_union(buffer_wgs84)
plot(calibration)
bbox <- st_bbox(calibration)


# read files from  Bioclim
files <- list.files('G:/Meu Drive/Modelagem/script/wc5', pattern='.bil$', full.names=TRUE)
predictors <- stack(files)

# read slope file
slope <- raster("./slope_res.tif")
slope

# read soil files
clay <- raster("./SoilGrids/clay_0_30cm.tif")
sand <- raster("./SoilGrids/sand_0_30cm.tif")
carbon <- raster("./SoilGrids/soc_0_30cm.tif")
nitro <- raster("./SoilGrids/nitrogen_0_30cm.tif")

# changing resolution
slope <- resample(slope, predictors, method = "ngb") 
clay <- resample(clay, predictors, method = "ngb")
sand <- resample(sand, predictors, method = "ngb")
carbon <- resample(carbon, predictors, method = "ngb")
nitro <- resample(nitro, predictors, method = "ngb")

# stacking
r1 <- stack(slope, predictors, clay, sand, carbon, nitro)

# crop predictors
layers <- crop(r1, bbox)
plot(layers$bio1)

#plotting layers
plot(layers$sand_0_30cm) # example with bio1
points(occ_thin[,1:2])

# remove variables that are highly correlated
spx <- extract(layers, occ_thin1)
spx <- data.frame(spx)
v <- vifstep(spx, th=5) 
v
bioclim <- exclude(layers, v)
bioclim <- stack(bioclim)


abr = "bolivica" #definindo abreviacao (or camura)

## format the data 
niche_data <- 
  BIOMOD_FormatingData(
    resp.var = occ_thin['sp1'],
    resp.xy = occ_thin[, c('Longitude','Latitude')],
    expl.var = bioclim,
    resp.name = abr,
    PA.nb.rep = 1,
    PA.nb.absences = 10000,
    PA.strategy = 'random')

# # k-fold selection
cv.k <- bm_CrossValidation(bm.format = niche_data,
                           strategy = "kfold",
                           nb.rep = 2,
                           k = 3)


opt.df <- bm_ModelingOptions(data.type = 'binary',
                             models = "MAXENT",
                             strategy = 'bigboss',
                             bm.format = niche_data,
                             calib.lines = cv.k)

niche_models <- 
  BIOMOD_Modeling(
    bm.format = niche_data,
    #CV.strategy = 'kfold',
    CV.nb.rep = 2,
    CV.perc = 0.8,
    models = c("MAXENT"),
    OPT.strategy = 'bigboss',
    bm.options = opt.df,
    var.import = 1,
    modeling.id = "demo1")


## formatted object summary
niche_data

bm_PlotVarImpBoxplot(bm.out = niche_models, group.by = c('expl.var', 'algo', 'algo'))

# Get evaluation scores & variables importance
get_evaluations(niche_models)
get_variables_importance(niche_models)

# Represent evaluation scores & variables importance
# see tables generated
bm_PlotEvalMean(bm.out = niche_models, metric.eval = c("TSS", "KAPPA")) #TSS (mean 2)
bm_PlotEvalBoxplot(bm.out = niche_models, group.by = c('algo', 'algo'))
bm_PlotEvalBoxplot(bm.out = niche_models, group.by = c('algo', 'run'))
bm_PlotVarImpBoxplot(bm.out = niche_models, group.by = c('expl.var', 'algo', 'algo'))
bm_PlotVarImpBoxplot(bm.out = niche_models, group.by = c('expl.var', 'algo', 'run'))
bm_PlotVarImpBoxplot(bm.out = niche_models, group.by = c('algo', 'expl.var', 'run'))


# Represent response curves
bm_PlotResponseCurves(
  bm.out = niche_models,
  models.chosen = get_built_models(niche_models),
  fixed.var = "median")


## do models projections ----
## current projections
niche_models_proj_current <- 
  BIOMOD_Projection(
    bm.mod = niche_models,
    new.env = bioclim,
    models.chosen = 'all',
    metric.binary = 'all',
    metric.filter = 'all',
    proj.name = "current",
    build.clamping.mask = FALSE,
    compress = FALSE)

#final files
current <- raster(niche_models_proj_current@proj.out@link[1])
plot(current) 

current_bin <- raster(niche_models_proj_current@proj.out@link[3])
plot(current_bin) 


