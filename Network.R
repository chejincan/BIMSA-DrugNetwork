rm(list = ls()) 
library(mvtnorm)
library(pbapply)
library(parallel)
library(orthopolynom)
library(glmnet)
library(ggplot2)
library(reshape2)
setwd("E:/")

get_legendre_par <- function(times,order) {
  get_interaction <- function(data,col){
    n <- nrow(data)
    clean_data <- data
    gene_list <- list()
    m <- clean_data[,col]
    M <- clean_data[,-col]
    x_matrix <- M
    x_matrix <- as.matrix(x_matrix)
    name <- colnames(clean_data)
    ridge1_cv <- cv.glmnet(x = x_matrix, y = m,alpha =0)
    best_ridge_coef <- abs(as.numeric(coef(ridge1_cv, s = ridge1_cv$lambda.min))[-1])
    
    fit_res <- cv.glmnet(x = x_matrix, y = m,alpha = 1,
                         penalty.factor =1/best_ridge_coef,
                         keep = TRUE)
    best_alasso_coef1 <- coef(fit_res, s = fit_res$lambda.min)
    
    gene_list_one <- list()
    gene_list_one[[1]] <- name[col]
    gene_list_one[[2]] <- best_alasso_coef1@Dimnames[[1]][best_alasso_coef1@i[-1]+1]
    gene_list_one[[3]] <- best_alasso_coef1@x[-1]
    gene_list[[col]] <- gene_list_one
    gene_list_one
    return(gene_list_one)
  }
  #module_relationship <- pblapply(1:k,function(c)get_interaction(t(cluster_mean),c))
  #par=init_curve_par
  cluster_mean <- SS
  rownames(cluster_mean) <- 1:k
  module_relationship <- pblapply(1:k,function(c)get_interaction(t(cluster_mean),c))
  #----------------------
  get_effect <- function(pars,effect,times,order,y0){
    if ( length(pars) != order ) {warning("legendre_pars != legendre_order")}
    LOP <-  legendre.polynomials(order, normalized=F) #Legendre polynomials
    d_LOP_fit <-  sapply(1:length(pars),function(c)
      pars[c]*polynomial.derivatives(LOP)[[c+1]])
    h <- scaleX(times,u=-1,v=1)[2]-scaleX(times,u=-1,v=1)[1] #per step h
    if(order==1){
      f <- function(x,y){d_LOP_fit}
    }else{
      f <- function(x,y){dy=do.call(sum,polynomial.values(polynomials=d_LOP_fit,x=x));dy}
    }
    #rk4 for legendre with step=h
    LOP_rk4 <- function(x0,y0){
      #f <- function(x,y){dy=do.call(sum,polynomial.values(polynomials=d_LOP_fit,x=x));dy}
      k1 <- f(x0,y0) 
      k2 <- f(x0+h/2,y0+h/2*k1)
      k3 <- f(x0+h/2,y0+h/2*k2)
      k4 <- f(x0+h,y0+h*k3)
      y <- y0+h/6*(k1+2*(1-1/sqrt(2))*k2+2*(1+1/sqrt(2))*k3+k4)
      return(y)
    }
    #dy_LOP, the increasment of each step
    dy <- sapply(1:length(times),function(c)LOP_rk4(scaleX(times,u=-1,v=1)[1],y0))
    #dy_LOP*y= main or sub effect
    dy_fit <- effect*c(0,dy[1:(length(times)-1)])
    return(cumsum(dy_fit))
  }
  
  ode_optimize <- function(pars, ind, dep, times, data, order, effect) {
    if(order==1){
      ind_pars <- matrix(pars,ncol=order)[1,]
      dep_pars <- matrix(pars[-1],ncol=order)
    }else{
      ind_pars <- matrix(pars,ncol=order)[1,]
      dep_pars <- matrix(pars,ncol=order)[-1,]
    }
    initial_value <- data[, ind][1]
    ind_effect <- get_effect(ind_pars, data[, ind], times, order, initial_value) + initial_value 
    
    if (is.null(nrow(dep_pars))) {
      dep_effect <- get_effect(dep_pars, data[, dep], times, order, 0)
      y <- ind_effect + dep_effect
    } else {
      dep_effect <- sapply(1:length(dep), function(c)
        get_effect(dep_pars[c,], data[, dep[c]], times, order, 0))
      y <- ind_effect + rowSums(dep_effect)
    }
    
    ssr <- sum((data[,ind]-y)^2)
    #add penalty
    alpha=1e-5
    ridge <- sum((data[,ind]-y)^2+alpha*(sum(ind_pars^2)+0*sum(dep_pars^2)))
    
     if(min(ind_effect)>0){
     return(ridge)
     }else{
     return(ridge*500)
   }
   
    return(ridge)
  }
  
  
  get_value <- function(effect,data,times,order){
    #input
    ind <- data[[1]]
    dep <- data[[2]]
    ind_no <- as.numeric(which(colnames(effect)==ind))
    dep_no <- as.numeric(sapply(1:length(dep), function(c) which(colnames(effect)==dep[c])))
    init_pars <- rep(0.1,(length(ind_no)+length(dep_no))*order)
    result <- optim(init_pars,ode_optimize,ind=ind_no,dep=dep_no,
                    times=times,data=effect,order=order,
                    method = "BFGS")
    #control=list(maxit=200000,trace=T)
    par_after <- matrix(result$par,length(ind)+length(dep),order)
    return(par_after)
  }
  # Nelder-Mead
  core.number <- detectCores()
  cl <- makeCluster(getOption("cl.cores", core.number))
  clusterEvalQ(cl, {library(orthopolynom)})
  clusterExport(cl, c("get_value","ode_optimize","get_effect","times","get_interaction",
                      "cluster_mean","module_relationship","order"),envir=environment())
  lop_par <- pblapply(1:nrow(cluster_mean),function(c)get_value(t(cluster_mean),
                                                                module_relationship[[c]],times,order),cl=cl)
  stopCluster(cl)
  return(list(lop_par,module_relationship))
}



all_lop_par <- get_legendre_par(times=times,order=order)



get_output <- function(relationship,par,effect,times,order){
  get_effect <- function(pars,effect,times,order,y0){
    if ( length(pars) != order ) {warning("legendre_pars != legendre_order")}
    LOP <-  legendre.polynomials(order, normalized=F) #Legendre polynomials
    d_LOP_fit <-  sapply(1:length(pars),function(c)
      pars[c]*polynomial.derivatives(LOP)[[c+1]])
    h <- scaleX(times,u=-1,v=1)[2]-scaleX(times,u=-1,v=1)[1] #per step h
    if(order==1){
      f <- function(x,y){d_LOP_fit}
    }else{
      f <- function(x,y){dy=do.call(sum,polynomial.values(polynomials=d_LOP_fit,x=x));dy}
    }
    #rk4 for legendre with step=h
    LOP_rk4 <- function(x0,y0){
      #f <- function(x,y){dy=do.call(sum,polynomial.values(polynomials=d_LOP_fit,x=x));dy}
      k1 <- f(x0,y0) 
      k2 <- f(x0+h/2,y0+h/2*k1)
      k3 <- f(x0+h/2,y0+h/2*k2)
      k4 <- f(x0+h,y0+h*k3)
      y <- y0+h/6*(k1+2*(1-1/sqrt(2))*k2+2*(1+1/sqrt(2))*k3+k4)
      return(y)
    }
    #dy_LOP, the increasment of each step
    dy <- sapply(1:length(times),function(c)LOP_rk4(scaleX(times,u=-1,v=1)[1],y0))
    #dy_LOP*y= main or sub effect
    dy_fit <- effect*c(0,dy[1:(length(times)-1)])
    return(cumsum(dy_fit))
  }
  output <- list()
  output[[1]] <- relationship[[1]]  
  output[[2]] <- relationship[[2]]  
  if(order==1){
    output[[3]] <- par[1,]
    output[[4]] <- matrix(par[-1],ncol=order)
  }else{
    output[[3]] <- par[1,]
    output[[4]] <- par[2:nrow(par),]
  }
  ind_no <- as.numeric(which(colnames(effect)==output[[1]]))
  dep_no <- as.numeric(sapply(1:length(output[[2]]), 
                              function(c) which(colnames(effect)==output[[2]][c])))
  inital_value <- effect[,ind_no][1]
  ind_effect <- get_effect(as.numeric(output[[3]]),effect[,ind_no],times,order,inital_value)+inital_value
  if (length(dep_no)==1) {
    dep_effect <- get_effect(as.numeric(output[[4]]),effect[,dep_no],times,order,0)
  }else{
    dep_effect <- sapply(1:length(dep_no), function(c)
      get_effect(as.numeric(output[[4]][c,]),effect[,dep_no[c]],times,order,0))
    colnames(dep_effect) <- dep_no
  }
  #------------
  all_effect <- cbind(ind_effect,dep_effect)
  effect_mean <- apply(all_effect,2,mean)
  output[[5]] <- effect_mean
  output[[6]] <- all_effect
  return(output)
}
#t(sapply(1:k, function(c)legendre_fit(cluster_result$curve_par[c,c(1:5)])))
cluster_mean <-SS
rownames(cluster_mean) <- 1:k
module_relationship <- all_lop_par[[2]]
net <- pblapply(1:k,function(c)get_output(module_relationship[[c]],all_lop_par[[1]][[c]],
                                          t(cluster_mean),times=times,order=order))
all_net <- net
#all_netM<-all_net


get_after <- function(i){
  temp <- matrix(NA,nrow = length(i[[2]]),ncol=3) 
  temp[,1] <- i[[2]]
  temp[,2] <- i[[1]]
  temp[,3] <- i[[5]][2:(length(i[[2]])+1)]
  
  colnames(temp) <- c('from','to','dep_effect')
  temp <- data.frame(temp)
  temp[,3] <- as.numeric(as.character(temp[,3]))
  return(temp)
}

a <- NULL
for (i in 1:k) {
  a <- rbind(a, get_after(all_net[[i]]))
}
a

aF<-a
aF$NewColumn <- ifelse(aF$dep_effect < 0, "-", "+")

colnames(aF)[4] <- "effect_type"

aF$dep_effect <- abs(aF$dep_effect)


word_mapping <- rownames(SS)
word_mapping

aF[[1]] <- word_mapping[match(aF[[1]], seq_along(word_mapping))]

aF[[2]] <- word_mapping[match(aF[[2]], seq_along(word_mapping))]


print(aF)


