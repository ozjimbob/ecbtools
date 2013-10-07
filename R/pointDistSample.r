pointDistSample=function(points,sample.size){
  cdist=crossdist(X=coordinates(points)[,1],Y=coordinates(points)[,2],x2=coordinates(points)[,1],y2=coordinates(points)[,2])
  cmn=apply(cdist,1,mean)
  mdist=cmn
  mdist=(ecbtools::normalize(mdist)/dim(ath)[1])
  selvec=1:(dim(ath)[1])
  grb=sample(selvec,size=sample.size,replace=F,prob=mdist)
  new=points[grb,]
  new
}

