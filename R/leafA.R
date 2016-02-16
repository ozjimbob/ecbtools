
leafa=function(file,cutoff=100,dpi=300)
{
  dots_per_cm = dpi / 2.54
  dot_side = 1/dots_per_cm
  dot_area = dot_side**2
  
  img<-raster::stack(file) # Load JPG as a raster stack (3 layers)
  k1=img@layers[1] # Isolate each layer into a separate raster
 
  
  img2=k1[[1]] 
  #m <- c(0, 1, 255,  2, cutoff-1, 0, cutoff-1, 255, 255)  # Define reclassification matrix
  #img=reclass(img,m) # Reclassify raster
 
  h=img2>cutoff # Set values less that cutoff to TRUE

  
  hc=h@data@values # Pull the raster contents out as a vector
  pres=sum(hc==0,na.rm=T) * dot_area
  return(pres)
}

leafa.dir=function(directory,pattern="jpg",cutoff=100,dpi=300){
  file_list=list.files(directory,pattern=pattern)
  output=data.frame(ID=seq_along(file_list),files=file_list,cover=NA)
  for(idx in 1:length(file_list)){
    print(idx/length(file_list))
    output$cover[idx]=darkpix(paste(directory,file_list[idx],sep=""),cutoff,dpi)
  }
  output
}


