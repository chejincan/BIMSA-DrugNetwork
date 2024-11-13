rm(list = ls())
library(glmnet)
library(parallel)
library(deSolve)
library(orthopolynom)
library(pbapply)
library(ggplot2)
library(patchwork)
library(cowplot)
library(mvtnorm)  
options(scipen = 10)
setwd("E:/")
data<-read.csv("NonWGpost.csv", stringsAsFactors = FALSE)

df_fs <- df_fs[,order(colSums(df_fs))]
index_fs <- colSums(df_fs)

df_fs <- log10(df_fs+1)
index_fs<-log10(index_fs+1)

power_equation <- function(par,x){
  y <- par[1]*x^par[2]
  y
}
power_par <- function(y,times){
  
  set.seed(200)
  x <- as.numeric(times)
  y <- as.numeric(y)
  lmFit <- lm( log( y + runif(1, min = 0, max = 0.1))  ~ log(x))
  coefs <- coef(lmFit)
  a <- exp(coefs[1])
  b <- coefs[2]
  
  tmp <- c(a,b)
  par_est <- function(par,x){
    sum( (y - power_equation(par,x))^2 )
  }
  r <- optim(tmp,par_est,x = as.numeric(times),method = "Nelder-Mead")
  return(r$par)
}

power_fit <- function(x, par){t(sapply(1:nrow(par),function(c)par[c,1]*x^par[c,2]))}

par.mu1 <- t(sapply(1:nrow(df_fs),function(c)power_par(df_fs[c,],index_fs)))

df1 <- data.frame(cbind(index_fs,t(df_fs)))

times1 <- seq(min(index_fs),max(index_fs),length=30)

line_data1 <- data.frame(cbind(times1,t(power_fit(times1,par.mu1))))

colnames(line_data1)[-1] <- rownames(df_fs)


