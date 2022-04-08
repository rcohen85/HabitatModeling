
# Some informative webpages:
# https://rspatial.org/raster/analysis/4-interpolation.html
# http://132.72.155.230:3838/r/spatial-interpolation-of-point-data.html
# https://mgimond.github.io/Spatial/interpolation-in-r.html

library(sp)
library(oce)
library(rgdal)
library(rspatial)
library(gstat)
library(automap)

method = "idw"
newx = seq(278,297,by=0.08)
newy = seq(24,44,by=0.08)

# Get data
load('J:/Chpt_3/Chla/0.0466deg/Chlorophyll_0_20170415.Rdata')
thisData = data.frame(data=stack(data.frame(data))[,1],
                      lat=rep(lats,length.out=length(lons)*length(lats)),
                      lon=rep(lons-360,each=length(lats)))
ggplot(thisData,aes(x=lon,y=lat))+geom_tile(aes(fill=data))

# Get rid of NAs
goodDat = which(!is.na(thisData$data))
thisData = thisData[goodDat,]

# Convert to a SpatialPointsDataFrame
coordinates(thisData) = ~lat+lon
proj4string(thisData) = CRS('+proj=longlat +datum=NAD83')

# Transform coordinate system
x = spTransform(thisData,CRS(paste(oceCRS('North Atlantic')," +datum=WGS84")))

if (method=="idw"){
  # Create raster to interpolate
  newCoords = data.frame(lat=rep(newy,length.out=length(newx)*length(newy)),
                         lon=rep(newx-360,each=length(newy)))
  coordinates(newCoords) = ~lat+lon
  proj4string(newCoords) = CRS('+proj=longlat +datum=NAD83')
  newGrid = spTransform(newCoords,CRS(paste(oceCRS('North Atlantic')," +datum=WGS84")))
  
  # Inverse Distance Weighting interpolation
  gs = gstat(formula=data~1,locations=x,maxdist=500)
  z = predict(gs,newGrid)
  
  newDat = data.frame(data=z@data[["var1.pred"]],
                      y=rep(newy,length.out=length(newx)*length(newy)),
                      x=rep(newx-360,each=length(newy)))
  ggplot(newDat,aes(x=x,y=y))+geom_tile(aes(fill=data))+labs(title="IDW, maxdist=500")

  
} else if (method=="krige"){
  # Create raster to interpolate
  newCoords = data.frame(lat=rep(newy,length.out=length(newx)*length(newy)),
                         lon=rep(newx-360,each=length(newy)))
  newGrid = SpatialGrid(points2grid(SpatialPoints(newCoords)),
                         proj4string=CRS(paste(oceCRS('North Atlantic')," +datum=WGS84")))
  
  # # Calculate empirical variogram       
  # gs = gstat(formula=data~1,locations=x)
  # v = variogram(gs,width=10)
  # fve <- fit.variogram(v, vgm(85, "Exp", 75, 20))
  # 
  # # Ordinary Kriging interpolation
  # k <- gstat(formula=data~1,locations=x, model=fve)
  # kp <- predict(k, template)
  # spplot(kp)
  
  # Carry out Kriging interpolation
  autoKrige(formula=data~1,input_data=x,new_data=newGrid)
  
}

