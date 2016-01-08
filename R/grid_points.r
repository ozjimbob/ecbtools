grid_points=function(points,spacing=25,gridsize=3){
  points$PointID=1:nrow(points)
  griddim=spacing*(gridsize-1)
  bias_st=-((gridsize-1)*.5)*spacing
  bias_end=bias_st+griddim
  
  xmat=matrix(rep(seq(bias_st,bias_end,spacing),gridsize),nrow=gridsize,ncol=gridsize)
  ymat=t(xmat)
  
  for(idx in seq_along(xmat)){
  
    this_points=points@data
    j=nrow(points)
    if(dim(sp::coordinates(points))[2]==3){
      sp::coordinates(this_points)=sp::coordinates(points)+c(rep(xmat[idx],j),rep(ymat[idx],j),rep(0,j))
    }else{
      sp::coordinates(this_points)=sp::coordinates(points)+c(rep(xmat[idx],j),rep(ymat[idx],j))
    }
    this_points$SubPointID=idx
  
    if(idx==1){
      new_points=this_points
    }else{
      new_points=rbind(new_points,this_points)
    }
  }
  new_points
}