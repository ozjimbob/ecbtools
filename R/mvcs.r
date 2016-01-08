subsamp <- function(de,ss){
  de = de[rownames(de) %in% sample(rownames(de),ss),]
  de
}

class_count=function(data,f_list){
  all_factors=list()
  for(col in seq_along(f_list)){
    all_factors[[col]] = levels(data[,f_list[col]])
  }
  factor_frame=expand.grid(all_factors,stringsAsFactors=F)
  names(factor_frame)=f_list
  
  factor_frame$count=0
  for(row in seq_along(factor_frame$count)){
    row_ss=factor_frame[row,1:3]
    paste_vec=c()
    for(col in seq_along(names(row_ss))){
      paste_vec[col]=paste(names(row_ss)[col],"==\"",row_ss[col],"\"",sep="")
    }
    search=paste(paste_vec,collapse=" & ")
    search=paste("frame_ss = subset(data,",search,")",sep="")
    eval(parse(text=search))
    n_rows=dim(frame_ss)[1]
    factor_frame$count[row]=n_rows
  }
  p_sum=sum(factor_frame$count)
  factor_frame$p = factor_frame$count / p_sum
  factor_frame
}

mvs=function(data,number,variables,iter=20){
  tally=0
  for(idx in 1:iter){
    df2=subsamp(data,number)
    a=class_count(data,variables)
    b=class_count(df2,variables)
    a_p = which(!a$count < 1)
    a_ct = a$count[a_p]
    b_ct = b$count[a_p]
    cs=chisq.test(a_ct, p = b_ct, rescale.p=T)
    cs_p=cs$p.value
    print(tally)
    if(cs_p>tally){
      out_frame=df2
      tally=cs_p
    }
  }
  out_frame
}



mvcs<-function(data,number,variables,iter=200){
  # Main Cramer test loop.
  max=0
  for(test in 1:iter){
    #cat(paste("Iteration ",test,": ",sep=""))
    #Generate a subsample of the data to test
    samp <- subsamp(data,number)
    
    d2 <- data
    
    # Create a matrix with all the variables in the full dataset in "a".  Here we use elevation and rainfall.
    a <- as.matrix(cbind(d2[variables]))
    
    # Create a matrix with the same variables from the sampled dataset.
    b <- as.matrix(cbind(samp[variables]))
    
    # Perform Cramer test on the two matrices
    c <- cramer::cramer.test(b,a,sim="ordinary")
    
    # See if the results of the Cramer test are the "best yet" - if so, keep this sample
    cat(paste(c$p.value,"\n",sep=""))
    if(c$p.value > max){
      print("New optimum.")
      max = c$p.value
      bestdata = samp
    }
  }
  
  # Write the table of the "best" sample to disk.
  bestdata
}

perlin_noise <- function( 
  vx = 7,   vy = 8,   
  ix = 100, iy = 100  
) {

  vector_field <- apply(
    array( rnorm( 2 * vx * vy ), dim = c(2,vx,vy) ),
    2:3,
    function(u) u / sqrt(sum(u^2))
  )
  f <- function(x,y) {

    i <- floor(x)
    j <- floor(y)
    stopifnot( i >= 1 || j >= 1 || i < vx || j < vy )

    v1 <- vector_field[,i,j]
    v2 <- vector_field[,i+1,j]
    v3 <- vector_field[,i,j+1]
    v4 <- vector_field[,i+1,j+1]

    u1 <- c(x,y) - c(i,j)
    u2 <- c(x,y) - c(i+1,j)
    u3 <- c(x,y) - c(i,j+1)
    u4 <- c(x,y) - c(i+1,j+1)
 
    a1 <- sum( v1 * u1 )
    a2 <- sum( v2 * u2 )
    a3 <- sum( v3 * u3 )
    a4 <- sum( v4 * u4 )

    s <- function(p) 3 * p^2 - 2 * p^3
    p <- s( x - i )
    q <- s( y - j )
    b1 <- (1-p)*a1 + p*a2
    b2 <- (1-p)*a3 + p*a4
    (1-q) * b1 + q * b2
  }
  xs <- seq(from = 1, to = vx, length = ix+1)[-(ix+1)]
  ys <- seq(from = 1, to = vy, length = iy+1)[-(iy+1)]
  outer( xs, ys, Vectorize(f) )
}
