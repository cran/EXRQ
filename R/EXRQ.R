#' Estimation for Quantile Power Transformation Model 
#' 
#' This function estimates the power transformation 
#' parameter at a single given quantile level
#' @details 
#' This function estimates the transformation parameter lam following 
#' the estimation method in Mu and He (2007) such that the conditional
#' quantile of the transformed response is linear in covariates.
#' The transformed response is defined as
#' \deqn{\Lambda(y) =
#' {(y+a)^\lambda-1}\lambda, if \lambda \neq 0; =log(y+a) if \lambda=0.
#' }
#' 
#' @param y a vector of length n representing the response
#' @param x a n x p matrix of n observations and p predictors
#' @param tau the quantile level of interest
#' @param lams a set of transformation parameters for grid search
#' @param a the location shift
#' @return A list is returned with the following components
#' @return lam: the estimated transformation parameter
#' @return coef: the estimated quantile coefficient from the power-transformed linear quantile regression
#' @references Mu, Y. and He, X. (2007). Power transformation toward a linear 
#' regression quantile. Journal of the American Statistical Association, 102, 269-279.
#' @export
#' @importFrom quantreg rq
#' 
#' 
PowT.1tau.func <- function(y, x, tau, lams= seq(-2, 2, 0.1), a){  
  n <- length(y)
  compare.x <- diag(n)
  for(i in 1:n){
    for(j in 1:n){
      compare.x[i,j] <- prod(x[i,]<x[j,])
    }
  }
  Vn = bhat <- NULL
  for(lam in lams){
    if(lam==0) {Lam.y <- log(y+a)} else {Lam.y <- ((y+a)^lam-1)/lam}
    idx.keep <- which(!is.na(Lam.y))
    Lam.y <- Lam.y[idx.keep]
    x2 <- x[idx.keep,]
    
    fit <- rq(Lam.y~x2, tau)
    res <- fit$res
    res <- round(res,10)
    bhat <- rbind(bhat, fit$coef)
    score <- tau - 1*(res<=0)
    Rn <- compare.x[idx.keep, idx.keep] * score
    Rn <- apply(Rn, 2, sum)/n
    Vn <- c(Vn, mean(Rn^2))
  }
  idx <- order(Vn)[1]
  lam <- lams[idx]
  coef <- bhat[idx,] 
  return(list(lam=lam, coef=bhat))
}

#' Quantile of the Pareto Distribution 
#' 
#' @param p the quantile level
#' @param gamma the shape parameter
#' @return the pth quantile
#' @export

qpareto = function(p, gamma)
{
  (1-p)^(-gamma)
}

#' Random Generation for the Pareto Distribution 
#' 
#' @param n number of observations
#' @param gamma the shape parameter
#' @return a vector of n i.i.d. random variables from the Pareto distribution
#' @export
rpareto = function(n, gamma)
{
  u = runif(n, 0, 1)
  (1-u)^(-gamma)
}

#' Hill Estimator of the Extreme Value Index
#' @references Chernozhukov, C., Fernandez-Val, I., and Galichon, A. (2010).
#' Quantile and probability curves without crossing. Econometrica, 78, 1093-1125.
#' @param x the estimated quantiles at intermediate quantile levels
#' @param taus the corresponding quantile levels
#' @param tol the tolerance level used for checking quantile crossing
#' @param min.prop the minimum proportion of quantiles that are estimated higher than the adjacent lower quantiles
#' @details 
#' The function estimates the extreme value index using Hill estimator based on the estimated 
#' intermediate quantiles. 
#' @export
#' @importFrom graphics plot
#' @importFrom quantreg rearrange
#' @importFrom stats stepfun
#' @return The estimated extreme value index is returned. If the proportion of cases with quantile crossing is too high, an NA is returned.

EVI.CFG.func = function(x, tol=1e-4, min.prop=0.3, taus)
{
  tt=mean(diff(x)>tol,na.rm=T);
  if(!is.na(tt)&tt>min.prop)
  {
    #rearrange based on Chernozhukov et al. (2006) to ensure monotonicity of quantile function
    idx = which(x>0 & !is.na(x))
    f=x[idx]
    taus2=taus[idx]
    tmp = plot(rearrange(stepfun(taus2[-length(taus2)],f),xmax=1))
    x=tmp$y
    out=mean(log(x/x[1]))
  }
  else out=NA
  return(out)
}




#' Two-Stage Extreme Conditional Quantile Estimator
#' 
#' This function provides the Two-Stage estimator in Wang, Li and He (2012) for conditional extreme quantiles based
#' on covariate-dependent extreme value index estimation. The intermediate conditional quantile is
#' estimated by quantile regression of the response on the original scale without any transformation.
#' The method is based on Hill estimator for the extreme value index and works for heavy-tailed distributions.
#' 
#' @references Wang, H., Li, D., and He, X. (2012). Estimation of high conditional quantiles for heavytailed 
#' distributions, Journal of the American Statistical Association, 107, 1453-1464.
#' @param y a vector of length n representing the response
#' @param x a n x p matrix of n observations and p predictors
#' @param tau.e the extreme quantile level of interest
#' @param xstar a m x p matrix of m observations and p predictors representing 
#' the covariate of interest
#' @param k the number of upper order statistics used in Hill estimator
#' @param tol the tolerance level used for checking quantile crossing
#' @importFrom quantreg rq
#' @importFrom grDevices dev.off pdf
#' @export
#' @return A list of the following commponents is returned
#' @return Q2Stage: the estimated (extrapolated) conditional extreme quantile of the response given x=xstar at the quantile level tau.e
#' @return gamma.x: the estimated covariate-dependent extreme value index (Hill estimator associated with x=xstar)
#' @seealso \code{\link{ThreeStage}}
#' @examples 
#' #A simulation example (sqrt transformation, heteroscedastic error)
#' library(EXRQ)
#' n=500
#' tau.e = c(0.99, 0.993, 0.995)
#' set.seed(12368819)
#' x1 = runif(n, -1, 1)   
#' x2 = runif(n, -1, 1)   
#' sqrty = 2 + x1 + x2 + (1+0.8*x1)*rpareto(n, 0.5)
#' x = as.matrix(cbind(x1, x2))
#' y = sqrty^2
#' xstar = rbind(c(-0.5,0),c(0,-0.5),c(0,0),c(0.5,0),c(0,0.5))
#' ## 2Stage method in Wang, Li and He (2012), no transformation
#' out.2stage <- TwoStage(y, x, xstar, tau.e, k=50)

TwoStage= function(y, x, xstar, tau.e, k, tol=1e-4)
{
  n = length(y)
  tau = 1-k/n
  p = ncol(x)
  nx = length(xstar)/p
  max.tau = (n-as.integer(n^(0.1)))/(n+1)
  taus = seq(tau, max.tau, length=k)
  rq1 = rq(y~x, taus)
  Q = cbind(1, xstar)%*%rq1$coef
  
  #estimate extreme value index at x=xstar
  pdf(file = NULL)
  gamma.x = apply(Q, 1, EVI.CFG.func, tol=tol, taus=taus)
  dev.off()
  
  #extrapolation
  Q2Stage = t(outer(((1-tau)/(1-tau.e)), (gamma.x), "^")) *Q[,1]
  
  out = list(Q2Stage=Q2Stage, gamma.x=gamma.x)
  return(out)
}


#' Three-Stage Extreme Conditional Quantile Estimator
#'
#' Provides the estimation of extreme conditional quantile using the three-stage estimation method in Wang and Li (2013). Specifically the function estimates
#' the tau.e-th conditional quantile of Y given x=xstar based on the power-transformed quantile regression model and extreme value theory. The method is based
#' on Hill estimator for the extreme value index and works for heavy-tailed distributions (on the original scale).
#' @param y a vector of n responses
#' @param x a n x p matrix of n observations and p predictors
#' @param xstar a m x p matrix of m observations and p predictors
#' @param tau.e the extreme quantile level of interest
#' @param grid.lam the set of lambda (transformation parameter) values for grid search
#' @param grid.k the grid for the number of upper order statistics involved in Hill estimator;
#' used for searching for the data-adaptive k. If the lenfth of grid.k is 1, then k is fixed at grid.k and no selection is performed.
#' @param tau.lam the quantile level used for estimating the transformation parameter
#' @param a location shift parameter in the power transformation (introduced to avoid negative y values)
#' @param tol the tolerance level for checking quantile crossing issue
#' @return A list is returned with the following components
#' @return lam: the estimated power-transformation parameter
#' @return k: the selected k, the number of upper order statistics involved in Hill estimator
#' @return gamma.x: the estimated x-dependent extreme value index (EVI)
#' @return cgmma: the pooled EVI estimation
#' @return Q3Stage: the three-stage estimator of the tau.e-th conditional quantile of Y given xstar based on the x-dependent EVI estimation
#' @return Q3StageP: the three-stage estimator of the tau.e-th conditional quantile of Y given xstar based on the pooled EVI estimation
#' @references Wang, H. and Li, D. (2013). Estimation of conditional high quantiles through power transformation. Journal of the American Statistical Association, 108, 1062-1074.
#' @importFrom quantreg rq
#' @export
#' @seealso \code{\link{TwoStage}}
#' @examples
#' #A simulation example (sqrt transformation, heteroscedastic error)
#' library(EXRQ)
#' n=500
#' tau.e = c(0.99, 0.993, 0.995)
#' set.seed(12368819)
#' x1 = runif(n, -1, 1)   
#' x2 = runif(n, -1, 1)   
#' sqrty = 2 + x1 + x2 + (1+0.8*x1)*rpareto(n, 0.5)
#' x = as.matrix(cbind(x1, x2))
#' y = sqrty^2
#' xstar = rbind(c(-0.5,0),c(0,-0.5),c(0,0),c(0.5,0),c(0,0.5))
#' ## 3Stage estimator
#' out.3stage <- ThreeStage(y, x, xstar, tau.e, grid.lam=seq(-0.5, 1.5, 0.1), grid.k=50, tau.lam=0.9)

ThreeStage = function(y, x, xstar, tau.e, grid.lam=seq(-2, 2, 0.1), grid.k, tau.lam, a=0, tol=1e-4)
{
  x = as.matrix(x)
  n = length(y)
  p = ncol(x)
  nx = length(xstar)/p
  
  max.tau = (n-as.integer(n^(0.1)))/(n+1)
  
  ###
  ######### Estimate the transformation parameter
  ###
  if(length(grid.lam)>1)
  {
    tmp = PowT.1tau.func(y, x, tau=tau.lam, lams= grid.lam, a)
    lam = tmp$lam
  }
  else if(length(grid.lam)==1) lam =grid.lam
  
  # the transformed Y
  if(lam==0) {Lam.y=log(y+a)} else {Lam.y = ((y+a)^lam-1)/lam}
  
  ###
  ########### If k=NULL, then select k among the grid points grid.k
  ###
  if(length(grid.k)==1) k=grid.k
  else if(length(grid.k)>1)
    k= select.k.func(y=y, x=x, Lam.y=Lam.y, lam=lam, a=a, max.tau=max.tau, grid.k=grid.k, n=n)
  
  ###
  ########## Obtain conditional quantiles of Y at the intermediate levels taus
  ###
  tau = 1-k/n
  taus = seq(tau, max.tau, length=k)
  rq1 = rq(Lam.y~x, taus)
  Lam.Q = cbind(1, xstar)%*%rq1$coef   #cond. quantile at the transformed scale
  
  ### calculate EVI on the original scale
  tt = est.gamma.func(taus=taus, Lam.Q, lam, a, tol)
  
  gamma.x = tt$gamma.x
  Q = tt$Q
  cgamma = mean(gamma.x, na.rm=T)
  
  ###
  ########## extrapolation to the extreme tails tau.e
  ###
  Q3Stage = t(outer(((1-tau)/(1-tau.e)), (gamma.x), "^")) *Q[1:nx,1]   # based on x-dependent EVI
  Q3StageP = outer(Q[1:nx,1], (((1-tau)/(1-tau.e))^cgamma), "*")        #based on pooled EVI
  
  out = list(lam=lam, k=k, Q3Stage=Q3Stage, Q3StageP=Q3StageP, gamma.x=gamma.x, cgamma=cgamma)
  return(out)
}


#' Selection of the Tuning Parameter k
#'
#' This function selects the tuning parameter k, the number of upper order statistics involved in Hill estimator of EVI among a grid of points
#' following the method described in Section 3.3 of Wang and Li (2013). The method selects k as the value that minimizes the discrepancy between
#' the estimated x-dependent EVI on the transformed scale and lam times the estimated x-dependent EVI on the original scale
#'
#' @param y a vector of n untransformed responses
#' @param x a n x p matrix of n observations and p predictors
#' @param Lam.y a vector of n power-transformed responses
#' @param lam the power-transformation parameter
#' @param a location shift parameter in the power transformation (introduced to avoid negative y values)
#' @param max.tau the upper bound of the intermediate quantile levels
#' @param grid.k the grid for the number of upper order statistics involved in Hill estimator
#' @param n the number of observations
#' @return the selected k is returned
#' @export
#' @references Wang, H. and Li, D. (2013). Estimation of conditional high quantiles through power transformation. Journal of the American Statistical Association, 108, 1062-1074.
#' @importFrom quantreg rq

select.k.func <- function(y, x, Lam.y, lam, a, max.tau, grid.k, n)
{
  obj = NULL
  grid.k = grid.k[grid.k>as.integer((1-max.tau)*n)+1]
  for(k in grid.k)
  {
    #estimate cond. quantiles of Lam.Y and Y given x=x
    tau = 1-k/n
    taus = seq(tau, max.tau, length=k)
    rq1 = rq(Lam.y~x, taus)
    Lam.Q = cbind(1, x)%*%rq1$coef   #cond. quantile at the transformed scale
    if(lam==0) {Q =exp(Lam.Q)-a} else {Q = (Lam.Q*lam+1)^(1/lam)-a}  # quantile at the original scale
    #Estimate the EVI of Y (Hill estimator)
    gamma.x = apply(Q, 1, function(x)  {x=sort(x[x>0]); mean(log(x/x[1]))})
    #Estimate the EVI of Lam.Y (Moment estimator)
    Mn1 = apply(Lam.Q, 1, function(x) {x=sort(x[x>0]); mean(log(x/x[1]))})  #Hill's estimator on transformed scale
    Mn2 = apply(Lam.Q, 1, function(x) {x=sort(x[x>0]); mean((log(x/x[1]))^2)})
    gamma.minus = 1-0.5/(1-Mn1^2/Mn2)
    gamma.star.x = Mn1 + gamma.minus #moment estimator (on the transformed scale)
    
    obj = c(obj, mean((lam*gamma.x-gamma.star.x)^2, na.rm=T))
  }
  idx = which.min(obj)
  k = grid.k[idx]
  return(k)
}

#' Estimation of the Extreme Value Index on the Original Scale
#'
#' This function estimates the extreme value index on the original scale based on the estimated intermediate conditional quantiles
#' on the transformed scale
#' @param taus a grid of intermediate high quantile levels
#' @param Lam.Q a vector of the same length as taus, representing the estimated intermediate conditional quantiles of Y (at taus) on the transformed scale
#' @param lam the power-transformation parameter
#' @param a location shift parameter in the power transformation (introduced to avoid negative y values)
#' @param tol the tolerance level for checking quantile crossing issue
#' @export
#' @importFrom stats quantile
#' @importFrom grDevices dev.off pdf
#' @return A list is returned with the following components.
#' @return gamma.x: the estimated EVI. If quantile crossing is too severe,
#' which suggests that the estimated intermediate conditional quantiles are unstable, then NA is returned.
#' @return Q: the estimated conditional quantile of Y on the original scale
#'
est.gamma.func = function(taus, Lam.Q, lam, a=0, tol)
{
  
  # quantile at the original scale
  if(lam==0) {Q = exp(Lam.Q)-a}
  else if(lam<0)
  {
    tt=Lam.Q*lam+1;
    tol2 = max(tol, quantile(as.vector(tt[tt>0]),0.002,na.rm=T))
    #if the tt for the baseline quantile is less than tol, discard this subject
    tt2 = 1*(abs(tt[,1])<tol)
    tt[which(tt2==1),] = NA
    #if more than 30% of other quantiles are less than tol, also discard this subject
    tt3 = rowMeans(tt<tol2)
    tt[which(tt3>0.3)] = NA
    #if other quantiles are less than tol, make those quantiles NA so they will not be used for calculating gamma
    tt[which(tt<tol2)] = NA
    Q = tt^(1/lam)-a
  }
  else if(lam>0)
  {
    tt=Lam.Q*lam+1;
    #if the tt for the baseline quantile is less than tol, discard this subject
    tt2 = 1*(abs(tt[,1])<tol)
    tt[which(tt2==1),] = NA
    #if more than 30% of other quantiles are less than tol, also discard this subject
    tt3 = rowMeans(tt<tol)
    tt[which(tt3>0.3)] = NA
    #if other quantiles are less than tol, make those quantiles NA so they will not be used for calculating gamma
    tt[which(tt<tol)] = NA
    Q = tt^(1/lam)-a
  }
  
  #estimate extreme value index at x=xstar (if quantile crossing is too severe then
  # return NA for the gamma.x)
  pdf(file = NULL)
  gamma.x = apply(Q, 1, EVI.CFG.func, tol=tol, taus=taus)
  dev.off()
  return(list(gamma.x=gamma.x, Q=Q))
}




#' Estimation of the C vector
#'
#' This function estimates the C vector involved in the function test.EVI for testing the constancy of EVI
#' @export
#' @param y a vector of n untransformed responses
#' @param x a n x p matrix of n observations and p predictors
#' @param tau an upper quantile level close to one
#' @param M a constant larger than one that is used for estimating the c vector and thus K(x) function. The default is two
#' @return A p-dimensional vector is returned.
#' @importFrom quantreg rq
#'
Estc.func = function(y,x,tau=0.99, M=2){
  #estimate the c vector; K(x)=(X^Tc)^{-EVI}
  #but this c estimation is very unstable
  beta2 = rq(y~x,1-M*(1-tau))$coefficients
  beta1 = rq(y~x, tau)$coefficients
  xbar = c(1, colMeans(x))
  c.est = (beta2-beta1)/as.numeric(xbar%*%(beta2-beta1))
  return(c.est)
}

#' Testing the Constancy of EVI Over Covariates
#'
#' This function tests whether the extreme value index of Y, gamma(x), is constant or varying across the covariate x by using the test procedure
#' described in Section 3.4 of Wang and Li (2013).
#' @export
#' @param y a vector of n untransformed responses
#' @param x a n x p matrix of n observations and p predictors
#' @param grid.lam a grid of points for power-transformation parameter
#' @param grid.k a grid of points for k, the number of upper order statistics involved in Hill estimator
#' @param tau.lam the quantile level used for estimating the transformation parameter
#' @param u.x the proportion to be trimmed in the x direction
#' @param a location shift parameter in the power transformation (introduced to avoid negative y values)
#' @param M a constant larger than one that is used for estimating the c vector and thus K(x) function. The default is two
#' @param tol  the tolerance level for checking quantile crossing issue
#' @return A list is returned with the following components
#' @return lam: the estimated power-transformation parameter
#' @return k: the selected tuning parameter k, the number of upper order statistics involved in Hill estimator
#' @return Tm: the proposed test statistic
#' @return scaledTm: the standardized test statistic
#' @return pval.iid: the p-value based on iid assumption, that is, assuming that K(x)=1
#' @return pval.nid: the p-value based on estimated K(x)=(X'C)^(1/EVI)
#' @return gamma.bar: the pooled EVI estimator
#' @return hat.gamma: a N-dimensional vector consisting of the estimated x-dependent EVI at x=xstar
#' @return xstar: a N x p matrix of N observations and p predictors
#' @importFrom quantreg rq
#' @importFrom mnormt rmnorm
#' @importFrom stats pchisq quantile runif var
#' @references  Wang, H. and Li, D. (2013). Estimation of conditional high quantiles through power transformation. Journal of the American Statistical Association, 108, 1062-1074.
#' @examples
#' library(EXRQ)
#' n=500
#' tau.e = c(0.99, 0.993, 0.995)
#' set.seed(12368819)
#' x1 = runif(n, -1, 1)   
#' x2 = runif(n, -1, 1)   
#' sqrty = 2 + x1 + x2 + (1+0.8*x1)*rpareto(n, 0.5)
#' x = as.matrix(cbind(x1, x2))
#' y = sqrty^2
#' out = testC.EVI(y, x, grid.lam=seq(-0.5, 1.5, 0.1), grid.k=50, tau.lam=0.9)
#' (Tval = out$scaledTm)
#' (pval.iid = out$pval.iid)
#' (pval.nid = out$pval.nid)


testC.EVI = function(y, x, grid.lam= seq(-2, 2, 0.1), grid.k, tau.lam=0.9, u.x=0, a=0, M=2, tol=1e-4)
{
  #same as test.func, but used k*mean() instead of k-(n-m)
  #test if the extreme value index of y is a constant or varying in x.
  #pretend that the errors are iid with K(z)=1 (so the limiting dist of the test statistic is simpler)
  
  
  ## outputs
  ##pval.iid: pvalue assuming K(X)=1
  ##pval.nid: by estimating K(X)
  
  n = length(y)
  p = length(x)/n
  max.tau = (n-as.integer(n^(0.1)))/(n+1)
  
  ###estimate the transformation function
  if(length(grid.lam)>1)
  {
    tmp = PowT.1tau.func(y, x, tau=tau.lam, lams= grid.lam, a)
    lam = tmp$lam
  }
  else if(length(grid.lam)==1)
    lam =grid.lam
  
  # the transformed Y
  if(lam==0) {Lam.y=log(y+a)} else {Lam.y = ((y+a)^lam-1)/lam}
  
  ### If k=NULL, then select k among the grid points grid.k
  if(length(grid.k)==1) k=grid.k
  else if(length(grid.k)>1)
    k= select.k.func(y, x, Lam.y, lam, a, max.tau, grid.k, n)
  xstar = x
  if(u.x>0)
  {
    idx = apply(xstar, 2, function(x) 1*(x<quantile(x,1-u.x) & x>quantile(x,u.x)))
    idx = apply(idx, 1, prod)
    idx = which(idx==1)
    xstar = xstar[idx,]
  }
  #estimate EVI of Y, on the original scale
  tau = 1-k/n
  taus = seq(tau, max.tau, length=k)
  rq1 = rq(Lam.y~x, taus)
  Lam.Q = cbind(1, xstar)%*%rq1$coef   #cond. quantile at the transformed scale
  
  ### calculate EVI on the original scale
  gamma.x = est.gamma.func(taus, Lam.Q, lam, a, tol)$gamma.x
  
  gamma.bar = mean(gamma.x, na.rm=T)
  
  mm = length(gamma.x)
  m = n - as.integer(n^(0.1))
  
  Tm = (k-(n-m))*mean((gamma.x-gamma.bar)^2, na.rm=T)
  Tm2 = Tm/(gamma.bar^2)
  
  (pval.iid = 1-pchisq(Tm2,p))
  ngamma = sum(!is.na(gamma.x))
  
  ### calculate p-value for noniid errors
  c.vec = Estc.func(Lam.y, x, tau=1-k/n, M)
  tau = 1-k/n
  xbar = colMeans(x)
  sigma = cbind(1,x)%*%c.vec
  min.sigma = quantile(sigma, 0.15)
  sigma = pmax(min.sigma, sigma)
  
  Z = cbind(1,x)
  tmp = Z*as.vector(sigma^(-1))
  Vw = var(tmp)
  H=t(tmp)%*%Z/n
  Hinv = solve(H)
  # obtain critival values by simulation
  V1 = t(Z)%*%Z/n
  w = rmnorm(n = 5000, mean = rep(0, p+1), V1)
  tmp2 = diag(w %*% Hinv %*% Vw %*% Hinv %*% t(w))
  (pval.nid = mean(tmp2>Tm2))
  
  out = list(lam=lam, k=k, Tm=Tm, scaledTm=Tm2, pval.iid=pval.iid, pval.nid=pval.nid, gamma.bar=gamma.bar, hat.gamma=gamma.x, xstar=xstar)
  
  return(out)
}


