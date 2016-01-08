div_pal=function (n, h = c(260, 0), c = 80, l = c(22, 90), power = 0.866666666666667, 
                  fixup = TRUE, gamma = NULL, ...) 
{
  if (!is.null(gamma)) 
    warning("'gamma' is deprecated and has no effect")
  if (n < 1) 
    return(character(0))
  h <- rep(h, length.out = 2)
  c <- c[1]
  l <- rep(l, length.out = 2)
  power <- rep(power, length.out = 2)
  rval <- seq(1, -1, length = n)
  rval <- hex(polarLUV(L = l[2] - diff(l) * abs(rval)^power[2], 
                       C = c * abs(rval)^power[1], H = ifelse(rval > 0, h[1], 
                                                              h[2])), fixup = fixup, ...)
  return(rval)
}

rad = function(x){x * (pi/180)}
deg = function(x){x * (180/pi)}

timespline=function(data,datevar,plotvar,k=40,hinge=F,col=div_pal(256)){
  year_list=unique(format(data[,datevar],"%Y"))
  x_len=365
  y_len=length(year_list)
  this_mat=matrix(NA,nrow=y_len,ncol=x_len)
  for(the_year_idx in seq_along(year_list)){
    the_year=year_list[the_year_idx]
    year_sub=subset(data,format(data[,datevar],"%Y")==the_year & format(data[,datevar],"%j")<366)
    year_sub[,datevar]=as.numeric(format(year_sub[,datevar],"%j"))
    this_mat[the_year_idx,year_sub[,datevar]]=year_sub[,plotvar]
  }
  rownames(this_mat)=year_list
  colnames(this_mat)=1:365

  as_table = melt(this_mat) # Now turn it into a data-frame
  if(hinge){
    as_table$Var2 = cos(rad((as_table$Var2/365)*360))
  }
  as_table=subset(as_table,!is.na(value))

  par(mfrow=c(2,1))

  quilt.plot(as_table,nx=length(year_list),ny=365,col=col(256))  # Plot of raw data

  fit=gam(value ~ s(Var1,Var2,k=k,bs="tp"),data=as_table)  # GAM is the quickest, nicest way to smooth. Couldn't get Tps to work.
  pred=predict(fit)
  quilt.plot(x=as_table$Var1,y=as_table$Var2,z=pred,nx=length(year_list),ny=365,col=col(256)) # Plot of smoothed model
}

