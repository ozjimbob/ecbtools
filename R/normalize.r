# Normalize

normalize=function(x,low=0,high=1){
  low+(x-min(x,na.rm=T))*(high-low)/(max(x,na.rm=T)-min(x,na.rm=T))
}

