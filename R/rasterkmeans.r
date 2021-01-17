
raster.kmeans=function(x,k=12,iter.max=100,nstart=10,geo=T,geo.weight=1){
  
  if("RasterStack" %in% class(x)){
    stk=x
  }else{
    stl=list.files(x,full.names = TRUE,include.dirs = FALSE,pattern = "tif")
    stk=raster::stack(stl)
  }
  
  
  oDF=as.data.frame(stk)
  oo=raster::xyFromCell(stk,1:(stk@ncols * stk@nrows))
  if(geo==T){
    oDF$x=oo[,1]
    oDF$y=oo[,2]
  }
  
  for(idx in 1:length(oDF)){
    oDF[,idx]=normalize(oDF[,idx])
  }
  
  if(geo==T){
    oDF$x = oDF$x * geo.weight
    oDF$y = oDF$y * geo.weight
  }
  
  oDFk=oDF
  oDFk$idx=1:length(oDFk[,1])
  oDFi=subset(oDFk,complete.cases(oDFk))
  oDF=subset(oDF,complete.cases(oDF))
  
  E <- kmeans(oDF, k, iter.max = iter.max, nstart = nstart)
  
  oDFi$cluster=E$cluster
  
  of=plyr::join(oDFk,oDFi,by="idx",type="left")
  
  EM <- matrix(of$cluster, nrow=stk@nrows,ncol=stk@ncols, byrow=TRUE)
  E.raster <- raster::raster(EM,crs=stk@crs, xmn=stk@extent@xmin, ymn=stk@extent@ymin, xmx=stk@extent@xmax, ymx=stk@extent@ymax)
  E.raster
}
