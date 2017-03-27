##### Model Results
cal.F <- read.csv("ModOut/resM.californicusChaturvediFALSE.csv") #No Pd
cal.T <- read.csv("ModOut/resM.californicusChaturvediTRUE.csv") #With Pd

##### Environmental Data
library(raster)
Mtemp <- raster("parameterFiles/wxMeanTempUSv1.asc")
RH <- raster("parameterFiles/wxRHUS.asc")
nights <- raster("parameterFiles/wxnightsUS.asc")

##### Read in Species Distributions 
library(rgdal)
map_MyYu <- readOGR("parameterFiles/ShapeFiles", "myotis_yumanensis")
map_MyCa <- readOGR("parameterFiles/ShapeFiles", "moyotis_californicus")

##### Read in World Map
library(ggplot2)
world_map <- map_data("world")

calRas.f <- survivalRast(cal.F)
calRas.t <- survivalRast(cal.T)

calplot.f <- surv.ploter(calRas.f, map_MyCa)
calplot.t <- surv.ploter(calRas.t, map_MyCa)

calhist <- diff.hist(Comp1.rast = calRas.f,Comp2PD.rast = calRas.t, dist.map =  map_MyCa,SpeciesName = "M. cal",key = "A",keylocX = 24, keylocY = 6)
