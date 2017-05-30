alts = function (x,  y, alpha1 = 0.5, alpha2 = 1.5, k1 = 15, k2 = 6, nn = TRUE, intercept = FALSE){
  ##################################################################################################
  ########################## Adapted Least Trimmed Square ##########################################
  ##################################################################################################
  ### INPUT
  ### x: Independent variables
  ### y: Dependent variables
  ### alpha1 : between 0 and 1
  ### alpha2 : larger than 1
  ### k1 : larger than 1
  ### k2 : larger than 1 and k2 <= k1
  ### nn : non-negative coefficients
  ### intercept : include intercept in the model
  ##################################################################################################
  ### OUTPUT 
  ### beta : estimation of coefficients
  ### number_outlier : number of outliers
  ### outlier_detect : index of detected outliers
  ### X.new : good observed points for independent variables
  ### Y.new : good observed points for dependent variables
  ### k1 : modified k1 (if the input value is not appropriate)
  ### k2 : modified k2 (if the input value is not appropriate)
  if(nn){
    library (nnls)
    if(intercept){
      m1 = nnls (cbind(1, x), y)
    }else{
      m1 = nnls (x, y)
    }
  }else{
    if(intercept){
      m1 = lm (y ~ x)
    }else{
      m1 = lm (y ~ x - 1) 
    }
  }
  res = abs (resid(m1))
  n   = length (y)
  order_id = order (res, decreasing = F)
  res_int  = res[order_id]
  Y_int = y[order_id]
  X_int = x[order_id, ]
  index1 = min (which (res_int > k1 * median (res_int)))
  k_up   = n - index1
  ######### modify k1 if too large ##########
  if(k_up <= 0) {
    while (k_up <= 0) {
      k1 = k1 - 1
      index1 = min (which (res_int > k1 * median (res_int)))
      k_up   = n - index1
    }
    if(k1 < 1){
      stop("there is no significant outlier")
    }
  }
  y_out_up = Y_int[c ((n - k_up + 1) : n)]
  k.low.ex = alpha1 * k_up
  k_low    = ceiling (k.low.ex)
  out_id   = (1:n > (n - k_low))
  kep_id   = (1:n < (n - k_low + 1))
  Y_new    = Y_int[kep_id]
  X_new    = X_int[kep_id, ]
  repeat{
    if(nn){
      if(intercept){
        m_new = nnls (cbind(1, X_new), Y_new)
        beta_new = m_new$x
        Y_hat    = cbind(1, x) %*% beta_new
      }else{
        m_new = nnls (X_new, Y_new)
        beta_new = m_new$x
        Y_hat    = x %*% beta_new
      }
    }else{
      if(intercept){
        m_new = lm (Y_new ~ X_new)
        beta_new = m_new$coefficients
        Y_hat    =  cbind (1, x) %*% beta_new
      }else{
        m_new = lm (Y_new ~ X_new - 1) 
        beta_new = m_new$coefficients
        Y_hat    = x %*% beta_new
      }
    }
    res_new  = abs (y - Y_hat)
    order_id = order (res_new, decreasing = F)
    res_new_ord = res_new[order_id]
    Y_ord    = y[order_id]
    X_ord    = x[order_id, ]
    index1   = min (which (res_new_ord > k2 * median (res_new_ord)))
    k_up     = n - index1
    if(k_up <= 0) {
      while (k_up <= 0) {
        k2 = k2 - 1
        index1 = min (which (res_int > k2 * median (res_int)))
        k_up   = n - index1
      }
      if(k2 < 1){
        stop("there is no significant outlier")
      }
    }
    k.low.ex = alpha2 * k.low.ex
    k_low    = min (ceiling(k.low.ex), k_up)
    out_id   = (1:n > (n - k_low))
    kep_id   = (1:n < (n - k_low + 1))
    Y_new    = Y_ord[kep_id]
    X_new    = X_ord[kep_id, ]
    if (prod(out_id == ((1:n) > n-k_up)) == 1){
      break
    }
  }
  if(nn){
    if(intercept){
      model = nnls (cbind(1, X_new), Y_new)
      }else{
      model = nnls (X_new, Y_new)
      }
    coefficients = model$x
  }else{
    if(intercept){
      model = lm (Y_new ~ X_new)
      }else{
      model = lm (Y_new ~ X_new - 1)
      }
    coefficients = model$coefficients
  }
  number_outlier = sum(out_id)
  outlier_id     = order_id[out_id]
  result         = list(beta = coefficients, number_outlier = number_outlier, outlier_detect = outlier_id, 
                        X.new = X_new, Y.new = Y_new, k1 = k1, k2 = k2)
  return(result)
}

tuningBIC = function (x, y, n, p, alpha1 = 0.1, alpha2 = 1.1, k1 = 10, up = 10, low = 0.8, nn = TRUE, 
        intercept = FALSE){
    ##################################################################################################
    ########################## Tuning Parameter using adjusted BIC ###################################
    ##################################################################################################
    ### INPUT
    ### x : independent variables
    ### y : dependent variables
    ### n : sample size
    ### p : number of independent variables
    ### up : upper bound for k2
    ### low : lower bound for k2
    ### alpha1 : between 0 and 1
    ### alpha2 : larger than 1
    ### k1 : larger than 1
    ### nn : non-negative coefficients
    ### intercept : include intercept in the model
    ###################################
    ### OUTPUT
    ### k2 : parameter k2 for ALTS
    para  = NULL
    BIC.alts  = NULL
    for (j in seq (low, up, 0.1)) {
        reg1  = alts (x = x, y = y, alpha1 = alpha1, alpha2 = alpha2, k1 = k1, k2 = j, nn = nn, intercept = intercept)
        if (intercept){
            res = reg1$Y.new - cbind (1, reg1$X.new) %*% reg1$beta
        }else{
            res = reg1$Y.new - reg1$X.new %*% reg1$beta
        }
        sse   = t (res) %*% res
        k     = reg1$number_outlier + p + 1
        BIC2  = (n - p) * log (sse / (n - p)) + k * (log (n-p) + 1)
        BIC.alts = rbind (BIC.alts, BIC2)
    }
    ind2    = which.min (BIC.alts)
    seq     = seq (low, up, 0.1)     
    k2      = seq[ind2]
    return (k2)
}


sample.sim = function (n, p, sig, a1, a2, nneg = TRUE, intercept = FALSE){
    ### n : sample size
    ### p : number of predictors
    ### sig : variance of y
    ### a1  : proportion of outliers
    ### a2  : proportion of leverage points in outliers
    ### nneg : non-negtive or not
    ### intercept : including intercept or not
    u     = matrix (runif (n * p, 0, 30), n, p)
    sigma = matrix (0.5, p, p)
    diag (sigma) = 1
    x     = u %*% (sigma ^ (1/2))
    e     = rnorm (n, 0, sig)
    tau   = matrix(0, n, 1)
    loc.a1  = sample (n, floor(a1 * n))
    tau[loc.a1] = rnorm (length(loc.a1), 20, 5)
    loc.a2  = sample (loc.a1, a2 * length(loc.a1))
    lev = 2 * max(x)
    x[loc.a2, ] = rnorm(length(loc.a2) * p, lev, 1)
    if (nneg){
        if (intercept){
            beta  = runif(p + 1, 0, 1)
            y     = cbind(1, x) %*% beta + tau + e
        }else{
            beta  = runif(p, 0, 1)
            y     = x %*% beta + tau + e
        }
    }else{
        if(intercept){
            beta  = rnorm(p + 1, 0, 1)
            y     = cbind(1, x) %*% beta + tau + e
        }else{
            beta  = rnorm(p, 0, 1)
            y     = x %*% beta + tau + e
        }
    }
    loc   = loc.a1
    result  = list (y = y, x = x, loc = loc, beta = beta)
    return(result)
}

##### Sample with 20 % outliers including 20% leverage point (non-negative)
##### generate sample with outlier
sam1 = sample.sim (n = 200, p = 5, sig = 1, a1 = 0.2, a2 = 0.2)
##### tuning parameter k2
k = tuningBIC (sam1$x, sam1$y, 200, 5, alpha1 = 0.1, alpha2 = 1.1, k1 = 10, up = 10, low = 0.8)
##### ALTS
m1 = alts(sam1$x,  sam1$y, alpha1 = 0.2, alpha2 = 1.1, k1 = 7, k2 = k)


##### Sample with 10 % outliers including 5% leverage point
##### generate sample with outlier
sam2 = sample.sim (n = 300, p = 10, sig = 1, a1 = 0.1, a2 = 0.05, nneg = FALSE)
##### tuning parameter k2
k = tuningBIC (sam2$x, sam2$y, 300, 10, alpha1 = 0.1, alpha2 = 1.1, k1 = 10, up = 10, low = 0.8, nn = FALSE, 
        intercept = FALSE)
##### ALTS
m2 = alts(sam2$x,  sam2$y, alpha1 = 0.2, alpha2 = 1.1, k1 = 7, k2 = k, nn = FALSE, intercept = FALSE)



##### Sample with 15 % outliers including 10% leverage point (with intercept)
##### generate sample with outlier
sam3 = sample.sim (n = 500, p = 15, sig = 1, a1 = 0.15, a2 = 0.1, nneg = FALSE, intercept = TRUE)
##### tuning parameter k2
k = tuningBIC (sam3$x, sam3$y, 500, 15, alpha1 = 0.1, alpha2 = 1.1, k1 = 10, up = 10, low = 0.8, nn = FALSE, 
        intercept = TRUE)
##### ALTS
m3 = alts(sam3$x,  sam3$y, alpha1 = 0.2, alpha2 = 1.1, k1 = 10, k2 = k, nn = FALSE, intercept = TRUE)


