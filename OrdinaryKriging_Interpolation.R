library(sp)
library(oce)
library(rgdal)
library(rspatial)
library(gstat)

newx = seq(278,297,by=0.08)
newy = seq(24,44,by=0.08)

# Get data
load('J:/Chpt_3/Chla/0.0466deg/Chlorophyll_0_20160201.Rdata')
thisData = data.frame(data=stack(data.frame(data))[,1],
                      lat=rep(lats,length.out=length(lons)*length(lats)),
                      lon=rep(lons-360,each=length(lats)))

# Get rid of NAs
goodDat = which(!is.na(thisData$data))
thisData = thisData[goodDat,]

# Convert to a SpatialPointsDataFrame
coordinates(thisData) = ~lat+lon
proj4string(thisData) = CRS('+proj=longlat +datum=NAD83')

# Transform coordinate system
x = spTransform(thisData,CRS(paste(oceCRS('North Atlantic')," +datum=WGS84 +units=km")))

# Create raster to interpolate
newCoords = data.frame(lat=rep(newy,length.out=length(newx)*length(newy)),
                     lon=rep(newx-360,each=length(newy)))
template = SpatialGrid(points2grid(SpatialPoints(newCoords)),
                       proj4string=CRS(paste(oceCRS('North Atlantic')," +datum=WGS84 +units=km")))

# Calculate empirical variogram       
gs = gstat(formula=data~1,locations=x)
v = variogram(gs,width=10)
fve <- fit.variogram(v, vgm(85, "Exp", 75, 20))

# Ordinary Kriging interpolation
k <- gstat(formula=data~1,locations=x, model=fve)
kp <- predict(k, template)
spplot(kp)
