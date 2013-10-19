ascGridResample=function(x,template,filename=""){
  if(class(template)[1]=="RasterLayer"){  
    template=template
  }else{
    template=raster(x)
  }
  
  proj4string(template)=CRS("+proj=longlat +datum=WGS84")
  if(class(x)[1]=="RasterLayer"){  
    the_raster_file=x
  }else{
    the_raster_file=raster(x)
  }
  proj4string(the_raster_file)=CRS("+proj=longlat +datum=WGS84")
  h=resample(the_raster_file,template)
  extent(h)=extent(template)
  tf=tempfile(fileext=".asc")
  writeRaster(h,tf,format="ascii",overwrite=T)
  
  header=file(template@file@name,"w+")
  head_lines=readLines(header,6)
  close(header)

  input_file=file(tf)
  data_lines=readLines(input_file)
  patched=c(head_lines,data_lines[7:(length(data_lines))])
  close(input_file)
  output_file=file(filename)
  writeLines(patched,output_file)
  close(output_file)
  unlink(tf)
  tf
}

