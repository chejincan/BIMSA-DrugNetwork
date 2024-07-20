rm(list = ls()) 
library(mvtnorm)
library(pbapply)
library(parallel)
library(orthopolynom)
library(glmnet)
library(ggplot2)
library(reshape2)

setwd("E:/")
load("WGpre.Rdata")


power_fit <- function(x, par){t(sapply(1:nrow(par),function(c)par[c,1]*x^par[c,2]))}
SS<-power_fit(times1,par.mu1)
data<- SS

times= times1

power_equation <- function(x, power_par){t(sapply(1:nrow(power_par),
                                                  function(c) power_par[c,1]*x^power_par[c,2] ) )}

power_equation_base <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  min_value = min(y[y != 0])
  
  lmFit <- lm(log(y + runif(1, min = 0, max = min_value)) ~ log(x))
  coefs <- coef(lmFit)
  a <- exp(coefs[1])
  b <- coefs[2]
  #, algorithm = "port"
  model <- try(nls(y ~ a * x^b, start = list(a = a, b = b), control=nls.control(maxiter= 5000)))
  
  if ('try-error' %in% class(model)) {
    result = NULL
  } else {
    result = model
  }
  return(result)
}

power_equation_all <- function(x, y, maxit = 1e2) {
  result <- power_equation_base(x, y)
  iter <- 1
  while (is.null(result) && iter <= maxit) {
    iter <- iter + 1
    try(result <- power_equation_base(x, y), silent = TRUE)
  }
  return(result)
}

get_init_par <- function(data, k, times) {
  n1 = length(times)
  
  # get initial pars based on k-means
  init_cluster <- kmeans(data, centers = k, algorithm="Lloyd",iter.max = 10000)
  
  cuM <- init_cluster$centers
  
  fit1 <- lapply(1:k, function(c) power_equation_all(times1, init_cluster$centers[c, 1:n1]))
  init_curve_par <- t(sapply(fit1, coef))
 
  init_SAD_par <- c(0.5,0.5)
  init_pro <- table(init_cluster$cluster) / nrow(data)
  
  return_object <- list(init_SAD_par, init_curve_par, init_pro)
  names(return_object) <- c("init_SAD_par", "init_curve_par", "init_pro")
  return(return_object)
}

#biFunClu_intial_pars <- get_init_par(data=data, k=3, times=times)
get_cluster <- function(data,k,input){ 
  Delta <- 100; iter <- 1; itermax <- 100;
  get_SAD1_covmatrix <- function(par,d){
    phi <- par[1]; gamma <- par[2]; 
    sigma <- array(dim=c(d,d))
    #formula 1, diag element
    diag(sigma) <- sapply(1:d, function(c)(1-phi^(2*c))/(1-phi^2) )
    #formula 2, non-diag element
    sigma[lower.tri(sigma)] <- unlist(lapply(1:(d - 1), function(c) phi^seq(1:(d - c)) * diag(sigma)[c]))
    
    # sigma[lower.tri(sigma)] <- do.call(c,lapply(1:(d-1),function(c)phi^seq(1:(d-c))*diag(sigma)[c]))
    sigma[upper.tri(sigma)] <- t(sigma)[upper.tri(t(sigma))]
    return(gamma^2*sigma)
  }
  power_equation <- function(x, power_par){t(sapply(1:nrow(power_par),
                                                    function(c) power_par[c,1]*x^power_par[c,2] ) )}
  
  
  mle <- function(par,data,prob){
    par1 <- par[1:2]
    par2 <- matrix(par[-c(1:2)],nrow = k)
    miu <- power_equation(times1,par2[1:k,])
    temp_S <- sapply(1:k,function(c)dmvnorm(data,miu[c,],get_SAD1_covmatrix(par1,d=30))*prob[c])
    LL <- sum(-log(rowSums(temp_S)))
    return(LL)
  }
  
  cat(paste0("Start biFunClu Calculation ","\n","Cluster_number=",k))
  while ( Delta > 1 && iter <= itermax ) {
    # initiation
    if(iter == 1){
      init_SAD_par <-input[[1]]
      init_curve_par <-input[[2]]
      pro <- input[[3]]
    }
    #E step, calculate the posterior probability
    old_par <- c(init_SAD_par,init_curve_par)
    LL_mem <- mle(par=old_par,data,prob=pro)
    miu <-power_equation(times1,init_curve_par[1:k,])
    
    mvn.c <- sapply(1:k, function(c) dmvnorm(data,miu[c,],get_SAD1_covmatrix(init_SAD_par,d=30))*pro[c])
    omega <- mvn.c/rowSums(mvn.c)
    #M step, calculate parameters
    pro <- colSums(omega)/sum(omega)
    
    new_par <- try(optim(old_par, mle, data = data, prob = pro, method = "Nelder-Mead", control = list(abstol = 1e-6, reltol = 1e-6)))
    
    if ('try-error' %in% class(new_par))
      break
    L_Value <- new_par$value
    init_SAD_par <- new_par$par[1:2]
    init_curve_par <- matrix(new_par$par[-c(1:2)],nrow = k)
    Delta <- abs(L_Value-LL_mem)
    if (Delta > 2000)
      break
    cat("iter=",iter,"LL=",L_Value,'\n')
    iter <- iter+1; LL_mem <- L_Value
  } 
  
  BIC <- 2*(L_Value)+log(nrow(data))*length(old_par)
  cat("Finish biFunClu Calculation")
  #plot-----------
  cluster <- apply(omega,1,which.max)
  clustered_ck <- data.frame(row.names(data),data[,1:30],cluster)
  clustered_data <- clustered_ck
  get_plot <- function(clustered_data){
    colnames(clustered_data) <- c("marker",seq(min(index_fs),max(index_fs),length=30),"cluster")
    long_df <- melt(clustered_data,c("marker","cluster"))
    colnames(long_df) <- c("marker","cluster","time","effect")
    p <-  ggplot()+geom_line(long_df,mapping=aes(as.numeric(as.character(time)),effect,group=marker,
                                                 colour= as.character(cluster)),alpha=1)+
      facet_wrap(long_df$cluster,scales = "fixed")+ 
      theme(legend.position="none") 
    #+ xlab("Time")+ylab("generic_effect")
    return(p)
  }
  
  p1 <- get_plot(clustered_ck)
  clustered_ck <- clustered_ck[,-1]
  return_object <- list(init_SAD_par,init_curve_par,pro,LL_mem,BIC,clustered_ck,p1) #cluster????
  cat("Finish biFunClu Calculation")
  names(return_object)<-c("SAD_par", "curve_par", "pro", "LL", 
                          "BIC", "clustered_ck","plot1")
  return(return_object)
  
}

set.seed(2024)
biFunClu_intial_pars <- get_init_par(data=data, k=n, times=times)
biFunClu_results<- get_cluster(data=data,k=n,input=biFunClu_intial_pars)




