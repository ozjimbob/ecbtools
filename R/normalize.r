# Normalize

normalize=function(x,low=0,high=1){
  low+(x-min(x))*(high-low)/(max(x)-min(x))
}

