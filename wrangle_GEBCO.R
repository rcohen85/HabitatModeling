library(ncdf4)
library(EFDR)
library(ggplot2)
library(viridis)
library(grid)
library(pracma)
library(geodist)


# If necessary, load original high-res file and downsample to desired resolution
# This step is slow (few hours)
# file = 'J:/Chpt_3/GEBCO/gebco_2021_n44.0_s24.0_w-82.0_e-63.0.nc'
# ncid = nc_open(file)
# 
# elev = t(ncvar_get(ncid,'elevation'))
# lat = ncvar_get(ncid,'lat')
# lon = ncvar_get(ncid,'lon')
# 
# # downsample data to coarser lat/lon grid
# highRes = data.frame(z=stack(data.frame(elev))[,1],
#                       y=rep(lat,length.out=length(lon)*length(lat)),
#                       x=rep(lon,each=length(lat)))
# 
# newx = seq(278,297,by=0.02)
# newy = seq(24,44,by=0.02)
# newGrid = regrid(highRes,n1=length(newx),n2=length(newy),method="idw")
# colnames(newGrid) = c("Lon","Lat","Depth")
# depthData = matrix(newGrid$Depth,ncol=length(newx),byrow=TRUE)
# save(newGrid,depthData,newx,newy,file=('J:/Chpt_3/GEBCO/DownsampledGrid_02deg.Rdata'))

# Load depth grid of desired resolution
load('J:/Chpt_3/GEBCO/DownsampledGrid_08deg.Rdata')

# Plot bathymetry and save
ggplot(newGrid,aes(x=Lon,y=Lat))+geom_tile(aes(fill=Depth))+ggtitle("0.08deg Grid")+scale_fill_viridis()
ggsave('J:/Chpt_3/GEBCO/08deg_BathMap.png',device="png")

# Calculate distances between lat grid points and lon grid points in meters
dy=as.numeric()
dx=as.numeric()

for (i in 1:length(newx)){
  latMat = data.frame(cbind(newy,newx[i]))
  colnames(latMat) = c('latitude','longitude')
  latDist = geodist_vec(latMat$longitude,latMat$latitude,sequential=TRUE,measure="geodesic")
  dy = cbind(dy,latDist)
}
colnames(dy) = as.character(newx)
rownames(dy) = as.character(newy[-length(newy)])


for (j in 1:length(newy)){
  lonMat = data.frame(cbind(newy[j],newx))
  colnames(lonMat) = c('latitude','longitude')
  lonDist = geodist_vec(lonMat$longitude,lonMat$latitude,sequential=TRUE,measure="geodesic")
  dx = cbind(dx,lonDist)
}
dx = t(dx)
colnames(dx) = as.character(newx[-length(newx)])
rownames(dx) = as.character(newy)

# Function to calculate components of unit normal vector at each depth grid point
# taken from insol package
cgrad_RC <-
  function(dem, dlx, dly, cArea=FALSE){
    if (nargs() < 1) {
      cat("USAGE: cgrad(dem, dx, dly=dlx, cArea=FALSE) \n")
      return()
    }
    if ("RasterLayer" %in% class(dem)) {
      dlx = raster::res(dem)[1]
      dly = raster::res(dem)[2]
      dem = raster::as.matrix(dem)
    }
    if (dlx == 0){
      cat("Input data is not a RasterLayer, then I need the DEM resolution dlx \n")
      return()
    }
    mm=as.matrix(dem)
    rows=nrow(mm)
    cols=ncol(mm)
    cellgr=array(dim=c(rows,cols,3))
    md=mm[-rows,-1]
    mr=mm[-1,-cols]
    mrd=mm[-1,-1]
    dlx=dlx[-rows,]
    dly=dly[,-cols]
    cellgr[-rows,-cols,2]=.5*dlx*(mm[-rows,-cols]+md-mr-mrd)
    cellgr[-rows,-cols,1]=.5*dly*(mm[-rows,-cols]-md+mr-mrd)
    cellgr[-rows,-cols,3]=dlx*dly
    #last row and col are undefined. Replicate last value from previous row/col
    cellgr[rows,,]=cellgr[(rows-1),,]
    cellgr[,cols,]=cellgr[,(cols-1),]
    cellArea=sqrt(cellgr[,,1]^2+cellgr[,,2]^2+cellgr[,,3]^2)
    if (cArea) return(cellArea) else {
      for (i in 1:3) cellgr[,,i] = cellgr[,,i]/cellArea
      return(cellgr)
    }
  }

# get x/y/z components of unit vector normal to each grid point
g = cgrad_RC(depthData,dlx=dx,dly=dy)

# calculate slope
floorSlope = rad2deg(acos(g[,,3]))
slopeDF = data.frame(z=stack(data.frame(floorSlope))[,1],
                       y=rep(newy,length.out=length(newx)*length(newy)),
                       x=rep(newx,each=length(newy)))
ggplot(slopeDF,aes(x=x,y=y))+geom_tile(aes(fill=z))+scale_fill_viridis()

# calculate aspect
floorAspect = (pi/2) - atan2(g[,,2],g[,,1])
floorAspect[floorAspect<0] = floorAspect[floorAspect<0]+(2*pi)
floorAspect = rad2deg(floorAspect)
aspectDF = data.frame(z=stack(data.frame(floorAspect))[,1],
                      y=rep(newy,length.out=length(newx)*length(newy)),
                      x=rep(newx,each=length(newy)))
ggplot(aspectDF,aes(x=x,y=y))+geom_tile(aes(fill=z))+scale_fill_viridis()

# Calculate locations of key isobaths



# Grab data points nearest to each HARP site, save time series
sites = c('HZ','OC','NC','BC','WC','NFC','HAT','GS','BP','BS','JAX')
HARPs = data.frame(t(data.frame(c(41.06165, -66.35155), # WAT_HZ
                                c(40.26333,-67.9861,40.22999, -67.97798),  # WAT_OC
                                c(39.83295, -69.98194),  # WAT_NC
                                c(39.19192, -72.22735),  # WAT_BC
                                c(38.37337, -73.36985),  # WAT_WC
                                c(37.16452, -74.46585),  # NFC
                                c(35.30183,-74.8789,35.5841,-74.7499),  # HAT_A & HAT_B
                                c(33.66992, -75.9977),   # WAT_GS
                                c(32.10527, -77.09067),  # WAT_BP
                                c(30.58295, -77.39002),  # WAT_BS
                                c(30.27818, -80.22085))))  # JAX_D
rownames(HARPs) = sites
colnames(HARPs) = c("Lat1", "Lon1", "Lat2", "Lon2")
HARPs$Lon1 = HARPs$Lon1+360
HARPs$Lon2 = HARPs$Lon2+360

slopeMat = matrix(nrow=11,ncol=2)
aspectMat = matrix(nrow=11,ncol=2)

for (m in 1:nrow(HARPs)){ # for each HARP site
  for (j in 1:2){
    # find data points nearest this HARP site
    eval(parse(text=paste('sitelat = which.min(abs(HARPs$Lat',j,'[m]-newy))',sep="")))
    eval(parse(text=paste('sitelon = which.min(abs(HARPs$Lon',j,'[m]-newx))',sep="")))
    
    # grab slope and aspect values at this HARP site
    slopeMat[m,j] = floorSlope[sitelat,sitelon]
    aspectMat[m,j] = floorAspect[sitelat,sitelon]
    
  }
}

save(slopeDF,aspectDF,slopeMat,aspectMat,file='J:/Chpt_3/GEBCO/SlopeAspect_08deg.Rdata')