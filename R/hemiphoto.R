
darkpix=function(file,cutoff=100,clip.circle=T)
{
  img<-raster::stack(file) # Load JPG as a raster stack (3 layers)
  k1=img@layers[1] # Isolate each layer into a separate raster
  k2=img@layers[2]
  k3=img@layers[3]
  
  ratioB = k3[[1]] / (k1[[1]] + k2[[1]])
 # bw = (k1[[1]] + k2[[1]] + k3[[1]]) / 3
  
  img=(k1[[1]] + k2[[1]] + k3[[1]]) / 3  # Calculate the mean of the three layers, put in new raster (img)
  #m <- c(0, 1, 255,  2, cutoff-1, 0, cutoff-1, 255, 255)  # Define reclassification matrix
  #img=reclass(img,m) # Reclassify raster
  ratioC=ratioB < .7
  h=img<cutoff # Set values less that cutoff to TRUE
  oo = h==1 & ratioC == 1
  
  h=oo
  
  if(clip.circle){
  # Mask circle
  
  # First determine the dimensions of the image
    vdim=h@nrows
    hdim=h@ncols
  
  # Now work out the centre
    midpx=hdim/2
    midpy=vdim/2
  
  # Create two matrices, i and j, with row/column coordinates stored in the cells.
    i=as.matrix(h)
    for(x in 1:hdim){
      i[,x]=x
    }
    j=as.matrix(h)
    for(y in 1:vdim){
      j[y,]=y
    }
  
  # Adjust matrices so the values are relative to the midpoint
    i=i-midpx
    j=j-midpy
  
    # Pythagoras - ci matirx now contains distance from midpoint
    ci=sqrt(i^2+j^2)
  
    # All points closer than image height, set to 1, otherwise set to NA
    cc=ci<vdim/(2)
    cc[cc==T]=1
    cc[cc==F]=NA
  
    # A lot of the raster code below is so we can plot nice images if we want to. 
    # We could just deal with raw matrices from this point, but it's not that much faster really.
    cc=raster::raster(cc)
    extent(cc)=raster::extent(h)
    # Clip circle
    h=h*cc
  }
  hc=h@data@values # Pull the raster contents out as a vector
  pres=sum(hc==1,na.rm=T)
  abs=sum(hc==0,na.rm=T)
  
  (pres/(pres+abs))*100  # TRUE / Total pixels * 100 
}

darkpix.dir=function(directory,pattern="jpg",cutoff=100,clip.circle=T){
  file_list=list.files(directory,pattern=pattern)
  output=data.frame(ID=seq_along(file_list),files=file_list,cover=NA)
  for(idx in 1:length(file_list)){
    print(idx/length(file_list))
    output$cover[idx]=darkpix(paste(directory,file_list[idx],sep=""),cutoff,clip.circle)
  }
  output
}

