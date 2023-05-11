
source("Functions_sdm.R")


####Petit script pour travailler avec blockCV jusqu'a l'arrivée sur biomod2

require(blockCV)
require(biomod2)
library(tidyverse)
library(sp)

#######Data : Je te laisse charger tes données
###Ici la colonne de Pres/absences s'appelle "Presence", les coordonnées "lon" et "lat"
###le stack de raster pour chaque variable s'appelle "predictors"

###On va rester sur la med cette fois
data <- read.csv("species_datasets/withICES/PENNRUB.CSV")
data_var <- extractVals(data)
predictors <- loadMedstack()

xmin <- -6
xmax <- 39
ymax <- 46
ymin <- 10
ext_study <- extent(xmin,xmax,ymin,ymax)
predictors <- crop(predictors,ext_study)
data_var <- dplyr::filter(data_var, lon > xmin & lon < xmax & lat > ymin &lat < ymax)

plot(predictors$bathymetry)
sp::coordinates(data_var) <- ~lon+lat
sp::proj4string(data_var) <- "+proj=longlat +datum=WGS84 +no_defs "


#######Methode simple : si je veut faire des blocs de taille précises (ex : 300 km)
sb <- cv_spatial(x = data_var, # sf or SpatialPoints of sample data (e.g. species data)
                 column = "Presence", # the response column (binary or multi-class)
                 r = predictors$bathymetry, # a raster for background (optional)
                 size = 300000, # size of the blocks in metres
                 k = 5, # number of folds
                 hexagon = TRUE, # use hexagonal blocks - defualt
                 selection = "random", # random blocks-to-fold
                 iteration = 100, # to find evenly dispersed folds
                 biomod2 = TRUE) # also create folds for biomod2
?cv_spatial
###le dernier


# On va voir ce que ça donne
cv_plot(cv = sb, # a blockCV object
        x = data_var, # sample points
        r =  predictors$bathymetry, # optionally add a raster background
        points_alpha = 0.5,
        nrow = 2)



#####Mais : Block CV peut faire l'analyse d'autocor spatiale
# exploring the effective range of spatial autocorrelation in raster covariates or sample data

range <- cv_spatial_autocor(
    x = data_var, # species data
    column = "Presence",
    predictors# column storing presence-absence records (0s and 1s)
    plot = TRUE
)
range$range/1000
###range$range est la taille de grille minimale por 
#####
# environmental clustering : on peut choisir les blocs selon les dissimilarité 
set.seed(6)

scv2 <- cv_nndm(
    x = data_var,
    column = "Presence",
    r = predictors,
    size = range$range, # range of spatial autocorrelation
    num_sample = 100, # number of samples of prediction points
    sampling = "regular", # sampling methods
    min_train = 0.1, # minimum portion to keep in each train fold
    plot = TRUE
)


########SI je veut aller sur biomod,
#les outputs de ces fonctions sont directement implémentable pour la CV
#
data_var <- as.data.frame(data_var)
BM_data <- BIOMOD_FormatingData(resp.name = "Penn",
                     resp.var = data_var$Presence,
                     expl.var = data_var[,14:17], 
                     resp.xy = data_var[,2:3])


# 4. Model fitting
biomod_model_out <- BIOMOD_Modeling(BM_data,
                                    models = c('GLM','MARS','GBM'),
                                    data.split.table = sb$biomod_table,##### Doing the CV
                                    var.import = 0,
                                    metric.eval = c('ROC'),
                                    do.full.models = TRUE)

biomod_model_eval <- get_evaluations(biomod_model_out)


biomod_model_eval[c("run", "algo", "metric.eval", "calibration", "validation")]
###On a des modele pas tres performant ... overfit discutable.