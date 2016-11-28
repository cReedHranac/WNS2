# packagesHaymanWNS.R

### Checks for required packages and installs if not present ###


packs <- rownames(installed.packages()) # names of all installed packages
req <- c("plyr","rgl","akima","tgp","shape","nlme","lattice","fields","lhs",
				 "RNCEP", "adehabitat", "rgdal", "raster", "maps", "maptools", "sp",
				 "ggplot2", "rgeos", "ggmap", "mapdata", "mapproj", "grid", "gridExtra",
				 "PBSmapping", "data.table", "deSolve", "knitr", "tibble", "rmarkdown") # packages required to run code
toInstall <- req[!is.element(req,packs)] # packages needing installation

if(length(toInstall)>0){
	if(verbose){
		print("The following packages are required and will be installed:")
		print(toInstall)
	}
	install.packages(toInstall)
}

# clean up workspace:
rm(packs,req,toInstall)
