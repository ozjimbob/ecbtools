BISE=function(x,slide_period=20,slope_threshold=0.2){
  slope_thres_value=0
  days=length(x)
  x=c(x,x,x)
  cor_x=rep(0,length(x))
 
  for(i in 2:(3*days)){
    if (x[i] >= x[i-1]){
      cor_x[i]=x[i]
    }else{
      if(i + slide_period > 3 * days){
        period=3*days-i-1
      }else{
        period=slide_period
      }
      slope_thres_value = x[i] + slope_threshold * (x[i-1]-x[i])
      bypassed_elems=0
      ndvi_chosen=0
      for(j in (i+1):(i+period-1)){
        if((x[j]>slope_thres_value) & (x[j] > ndvi_chosen)){
          ndvi_chosen=x[j]
          bypassed_elems=j-i
        }
        if(ndvi_chosen >= x[i-1]){break}
      }
      if(ndvi_chosen==0){
        cor_x[i]=x[i]
      }else{
        for(j in 1:bypassed_elems){
          cor_x[i-1+j]=-1
        }
        i=i+bypassed_elems
        cor_x[i]=ndvi_chosen
      }
    }
    
  }
  cor_x=cor_x[days:(days*2)]
  
  for(i in 2:length(cor_x)){
    if(cor_x[i]==-1)
      cor_x[i]=cor_x[i-1]
  }
  cor_x
}
