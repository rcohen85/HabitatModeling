

thisData = data.frame(data=stack(data.frame(data))[,1],
                      lat=rep(lats,length.out=length(lons)*length(lats)),
                      lon=rep(lons-360,each=length(lats)))
