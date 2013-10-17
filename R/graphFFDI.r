#FFDI Plot

graphFFDI=function(state="N",station="94757",DF=9){
  options(digits=14,scipen=20)
  URL=paste("http://www.bom.gov.au/fwo/ID",state,"60901/ID",state,"60901.",station,".axf",sep="")
  data=read.csv(URL,skip=19)
  
  TEMP=data$air_temp
  RH=data$rel_hum
  WS=data$wind_spd_kmh
  
  FFDI=ecbtools::index_FFDI_DF(TEMP,rep(DF,length(RH)),RH,WS) 
  FFDI=rev(FFDI)
  
  #dt=as.character(data$local_date_time_full.80.)
  dt=data$local_date_time_full.80.
  dt2=paste(substr(dt,1,4),"-",substr(dt,5,6),"-",substr(dt,7,8),"-",substr(dt,9,10),"-",substr(dt,11,12),"-",substr(dt,13,14),sep="")
  dt3=as.POSIXct(dt2,format="%Y-%m-%d-%H-%M-%S")
  dt3=rev(dt3)
  df=data.frame(Time=dt3,FFDI=FFDI)
  plot(df$FFDI ~ df$Time,xlab="Time",ylab="FFDI",type="l")
  rect(min(df$Time,na.rm=T),min(df$FFDI,na.rm=T),max(df$Time,na.rm=T),5,col="#00FF0044")
  rect(min(df$Time,na.rm=T),5,max(df$Time,na.rm=T),12,col="#66FF0044")
  rect(min(df$Time,na.rm=T),12,max(df$Time,na.rm=T),25,col="#FFFF0044")
  rect(min(df$Time,na.rm=T),25,max(df$Time,na.rm=T),50,col="#FF660044")
  rect(min(df$Time,na.rm=T),50,max(df$Time,na.rm=T),75,col="#FF000044")
  rect(min(df$Time,na.rm=T),75,max(df$Time,na.rm=T),100,col="#FF006644")
  rect(min(df$Time,na.rm=T),100,max(df$Time,na.rm=T),200,col="#FF00FF44")
  df
}