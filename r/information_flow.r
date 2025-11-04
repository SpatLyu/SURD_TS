'
multi_causality_est <- function(X, max_lag=1, np=1, dt=1, series_temporal_order=NULL, significance_test=1) 
By  Yineng Rong (yinengrong@foxmail.com)


source("LK_Info_Flow.R")
IF_result=multi_causality_est(X, np=1, dt=1, max_lag=1, series_temporal_order=NULL, significance_test=1):

Estimate information flow, the information transfer from columns to raws among
M time series (stored as a MXM matrix).

% On input:
%    X: matrix storing the M time series (each as Nx1 column vectors)
%    max_lag: time order of lags (default 1)
%    >=1 lag time step
%    -1 Determine the lag order of the data based on the AIC criterion.
%   - 2 Determine the lag order of the data based on the BIC criterion.
%    np: integer >=1, time advance in performing Euler forward 
%	 differencing, e.g., 1, 2. Unless the series are generated
%	 with a highly chaotic deterministic system, np=1 should be
%	 used. 
%    (default 1)
%    dt: frequency of sampling (default 1)
%
%    series_teporal_order: Nx1 column vectors, records the timestamp of 
%    each sample, with a minimum sampling interval of dt (used for 
%    panel data, or missing measurement data).
%    (default [])
%    significance_test:  1  do the significance test (default)
%           =0  not (to save the computation time)
%
% On output:
%    a structure value IF_result with sub
%    IF:  information flow
%    nIF: normalized information flow
%    SEIF: standard error of information flow
%    errr.e90/e95/e99:standard error at 90/95/99% confidence level
%    dnoise: dH1_noise/dt with correlated noise
%    dnoise_old: dH1_noise/dt without correlated noise (##  been commented out.)
%    nIF: normalized information flow without correlated noise
%    p: p-value of information flow

% Citations: 
%    X.S. Liang, 2016: Information flow and causality as rigorous notions
%    ab initio. Phys. Rev. E, 94, 052201.
%    X.S. Liang, 2014: Unraveling the cause-effect relation between time series. Phys. Rev. E 90, 052150.
%    X.S. Liang, 2015: Normalizing the causality between time series. Phys. Rev. E 92, 022126.
%    X.S. Liang, 2021: Normalized Multivariate Time Series Causality Analysis and Causal Graph Reconstruction. Entropy. 23. 679.
% 
% Note: This is an alternative and simplified version of
%   multi_causality_est.m, which was originally coded by X.S. Liang
%   Here all the causal relations are inferred once for all.
% 
% Author: Yineng Rong (yinengrong@foxmail.com)
'  

library(stats)
library(MASS)

multi_causality_est <- function(X, max_lag=1, np=1, dt=1, series_temporal_order=NULL, significance_test=1) {
    res <- list()
    dm <- dim(X)
#    m <- dm[1]
#    n <- dm[2]
  m = nrow(X)
  n = ncol(X)


  if (is.null(series_temporal_order)) {
    series_temporal_order <- seq(from = 1, to = (dt * n), by = dt)
  }

  if (max_lag < 0) {
    max_lag = lag_order(X, series_temporal_order/dt, max_lag)
  }

  q1 = max_lag + 1
  XX = array(0, dim = c(m, q1, n+q1-1))

  for (k in 1:q1) {
    XX[, k, k:(k+n-1)] <- X
  }

  XX <- XX[, , 1:n]


  panel_data = pre_panel_data(XX, series_temporal_order/dt, max_lag)
  XL = panel_data$X0
  XR = panel_data$XL

  m2 = nrow(XL)
  n2 = ncol(XL)


  X1 = rbind(XR, matrix(1, nrow = 1, ncol = n2))
  if (max_lag>1){
    X2 = rbind(XL, XR[1:(m2 * (max_lag-1)), ], matrix(1, nrow = 1, ncol = n2))
    }
    else {
      X2 = rbind(XL, matrix(1, nrow = 1, ncol = n2))
    }




  A0 = X2 %*% ginv(X1)
  A = (A0 - diag(m2*max_lag+1)) / dt
  At = A[1:(nrow(A)-1), 1:(ncol(A)-1)]
  C = cov(t(X1[1:(nrow(X1)-1), ]))#cov(X1[1:(nrow(X1)-1), ])


  IF = ((C / diag(C))) * At

  E = X2 - A0 %*% X1
  SIG = E %*% t(E) / (n - np - m2 - 1)

  B = sqrt(abs(SIG[1:(nrow(SIG)-1), 1:(ncol(SIG)-1)]) / dt)
  dH1_noise = colSums(B^2) / (2 * diag(C))
  Z = rowSums(abs(IF)) + abs(dH1_noise)
  Z = matrix(rep(Z, m2*max_lag), nrow = m2*max_lag, ncol = m2*max_lag)
  nIF = IF / Z




  if (significance_test == 1) {
    diag_SIG <- diag(SIG[-nrow(SIG), -ncol(SIG)])
    reshaped_sig <- matrix(diag_SIG, ncol = 1, nrow = m2*max_lag)
    diag_x <- diag(solve(X1[-nrow(X1), ] %*% t(X1[-nrow(X1), ])))
    reshaped_x <- matrix(diag_x,ncol = m2*max_lag, nrow = 1)
    se_a <- sqrt(reshaped_sig %*% reshaped_x)


    SE_IF <- se_a * abs(C / diag(C)) / dt  # standarized error of IF
    p <- (1 - pnorm(abs(IF / SE_IF))) * 2  # p-value

    z90 = 1.65
    z95 = 1.96
    z99 = 2.56

    e90_IF = SE_IF * z90
    e95_IF = SE_IF * z95
    e99_IF = SE_IF * z99
    res=list('IF' = temporal_lag(t(IF), m2,max_lag), 
         'nIF' = temporal_lag(t(nIF), m2,max_lag), 
         'dnoise' = dH1_noise[1:(length(dH1_noise)-m2)], 
         'max_lag' = max_lag, 
         'SEIF' = temporal_lag(t(SE_IF), m2,max_lag), 
         'err_e90' = temporal_lag(t(e90_IF), m2,max_lag), 
         'err_e95' = temporal_lag(t(e95_IF), m2,max_lag), 
         'err_e99' = temporal_lag(t(e99_IF), m2,max_lag),  
         'p' = temporal_lag(t(p), m2,max_lag)
    )
  } else {
    res=list('IF' = temporal_lag(IF, m2,max_lag), 
         'nIF' = temporal_lag(nIF, m2,max_lag), 
         'dnoise' = temporal_lag(dH1_noise, m2,max_lag), 
         'max_lag' = max_lag 
    )
  }

  return(res)
}

pre_panel_data <- function(XX, t, q){
  dims = dim(XX)
  n = dims[1]
  m = dims[3]

  if(length(t) == 0){
    t = 1:m
  }

  tt = matrix(-pi^2, nrow = q+2, ncol = m+q)

  for(k in 1:(q+1)){
    tt[k, k:(k+m-1)] = t
  }

  tt[q+2, ] = tt[q+1, ] - 1
  dt = abs(diff(tt) + 1)
  max_dt = apply(dt, 2, max)
  dims <- c(1, length(max_dt))
  max_dt <- array(max_dt,dim=dims)
  ind = which(max_dt <0.0000001)
    dims <- c(1, length(ind))
  ind <- array(ind,dim=dims)
  
  
  M = length(ind)

  nq = n * q

  X0 = matrix(nrow = n, ncol = M)

  for(i in 1:n){
    X0[i, ] = XX[i, 1, ind]
  }
  XL = array(dim = c(nq, M))

  sub_XX <- XX[, 2:(q+1), ind]
  #reordered_XX <- aperm(sub_XX, c(2, 1, 3)) 
  XL <- array(data = as.vector(sub_XX), dim = c(nq, M))
  
  return(list(X0=X0, XL=XL))
}



infocrit <- function(L, k, m){
    if(m - k - 1 <= 0){
        aic <- NaN
    } else {
        aic <- -2 * L + 2 * k * (m/(m - k - 1))
    }
    bic <- -2 * L + k * log(m)
    return(list(aic=aic, bic=bic))
}



lag_order <- function(X, t, option) {

  n <- nrow(X)
  m <- ncol(X)

  X <- sweep(X, 1, rowMeans(X))

  morder <- 1:min(max(floor(m/(n^2 + n)), 1), 20)
  nummo <- length(morder)

  aic <- rep(NaN, nummo)
  bic <- rep(NaN, nummo)

  q <- max(morder)
  q1 <- q + 1

  XX <- array(0, c(n, q1, m + q))
  for (k in 0:q) {
    XX[, k + 1, (k + 1):(k + m)] <- X
  }

  for (i in 1:nummo) {
    q <- morder[i]

    if (q >= m) {
      print(paste('WARNING: model order too large (must be <', m, ')'))
      next
    }

    XX <- XX[, , 1:m]

    pre_panel_data_output <- pre_panel_data(XX, t, q)
    X0 <- pre_panel_data_output$X0
    XL <- pre_panel_data_output$XL

    M <- ncol(X0)

    A <- ginv(t(XL)) %*% t(X0)
    E <- X0 - (XL %*% A)

    DSIG <- det(t(E) %*% E / (M - 1))

    if (DSIG <= 0) {
      print('WARNING: residuals covariance not positive definite')
      next
    }

    infocrit_output <- infocrit(-(M/2) * log(DSIG), q*n*n, M)
    aic[i] <- infocrit_output$aic
    bic[i] <- infocrit_output$bic
  }

  if (option == -1) {
    max_lag <- min(which(aic == min(aic, na.rm = TRUE)))
  } else {
    max_lag <- min(which(bic == min(bic, na.rm = TRUE)))
  }

  return(max_lag)
}

temporal_lag <- function(X, m1, max_lag) {
  # Removing causality from future to past when max_lag > 1
  #return(array(X[1:m1,], dim = c(m1, m1, max_lag)))
  XX <- aperm((array(X[1:m1,], dim = c(m1, max_lag, m1))),c(1,3,2))
  return (XX)
}