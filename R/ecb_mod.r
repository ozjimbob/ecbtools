ecb_mod=function(data,
                 term.names,
                 resp,
                 random.def,
                 model.family=gaussian(),
                 model.type="lm"){
  require(combinat)
  require(lme4)
  require(MASS)
  
  
  n.terms <- length(term.names)
  resp.mod <- paste(resp,"~",sep="")
  n.resp <- rep(resp.mod,n.terms)
  
  term.seq <- seq(1,n.terms,1)
  n.mod <- 2^n.terms
  
  l.dim <- 0
  for (m in 1:(n.terms)) {
    l.dim[m] <- length(combn(term.names,m,simplify=F))
  }
  
  mod_list=c()
  
  for (y in 1:n.terms) {
    lll  <- combn(term.names,y,simplify=F)
    col_list=sapply(lll,paste,collapse="+")
    mod_list=c(mod_list,col_list)
    
  }
  
  mod.vec= mod_list
  mod.vec <- paste0(resp.mod,mod.vec)
  mod.vec=mod.vec[!is.na(mod.vec)]
  mod.vec=c(paste0(resp.mod,"1"),mod.vec)
  
  # Set functions
  AICc <- function(...) {
    models <- list(...)
    num.mod <- length(models)
    AICcs <- numeric(num.mod)
    ns <- numeric(num.mod)
    ks <- numeric(num.mod)
    AICc.vec <- rep(0,num.mod)
    for (i in 1:num.mod) {
      if (length(models[[i]]$df.residual) == 0) n <- models[[i]]$dims$N else n <- length(models[[i]]$residuals)
      if (length(models[[i]]$df.residual) == 0) k <- sum(models[[i]]$dims$ncol) else k <- (length(models[[i]]$coeff))+1
      AICcs[i] <- (-2*logLik(models[[i]])) + ((2*k*n)/(n-k-1))
      ns[i] <- n
      ks[i] <- k
      AICc.vec[i] <- AICcs[i]
    }
    return(AICc.vec)
  }
  
  delta.AIC <- function(x) x - min(x) ## where x is a vector of AIC
  weight.AIC <- function(x) (exp(-0.5*x))/sum(exp(-0.5*x)) ## Where x is a vector of dAIC
  #ch.dev <- function(x) ((( as.numeric(x[12]) - as.numeric(x[10]) )/ as.numeric(x[12]))*100) ## % change in deviance, where x is glm object
  
  
  #Load ModelData
  Modnum <- length(mod.vec)
  n <- dim(data)[1]
  
  
  ###################################################################################
  # Model fitting and logLik output loop
  
  
  FitList <- LikList <- SaveCount <- AICc.vec <- 0
  results <- data.frame(mod.vec)
  dummy.mat <- matrix(0,ncol=(log(Modnum))/(log(2)),nrow=Modnum)
  
  for(i in 1:Modnum) {
    
    if (model.type == "lm") fit <- lm(as.formula(mod.vec[i]), data=data)
    if (model.type == "glm") fit <- glm(as.formula(mod.vec[i]), data=data,family=model.family)
    if (model.type == "lmer") fit <- lmer(as.formula(paste0(mod.vec[i]," + ",random.def)),data=data,REML= FALSE)
    if (model.type == "glmer") fit <- glmer(as.formula(paste0(mod.vec[i]," + ",random.def)),data=data,family=model.family)
    if (model.type == "glmer.nb") fit <- glmer.nb(as.formula(paste0(mod.vec[i]," + ",random.def)),data=data)
    
    AICc.mod <- AIC(fit)
    
    AICc.vec[i] <- AIC(fit)
    
    for (k in 1:(length(term.names))) {
      if (length(grep(term.names[k],mod.vec[i])) == 1) dummy.mat[i,k] <- 1
    }
    
    FitList[i] <- list(fit)
    LikList[i] <- logLik(fit)[1]
    #print(i)
  }
  
  
  # AIC weights
  dAICc <- delta.AIC(AICc.vec)
  wAICc <- weight.AIC(dAICc)
  
  # w+
  w.plus <- dummy.mat[,]*wAICc	
  w.plus.vec <- 0
  for (t in 1:(length(term.names))) {
    w.plus.vec[t] <- sum(w.plus[,t])
  }
  
  
  #mod.vec
  
  output_table=data.frame(Model=mod.vec,LL=LikList,AIC=AICc.vec,dAIC=dAICc,wAIC=wAICc)
  output_table=output_table[order(-output_table$wAIC),]
  w_table=data.frame(Term=term.names,Wplus=w.plus.vec)
  
  print(output_table)
  print(w_table)
  
  list(AICtable=output_table,Wtable=w_table)
}
