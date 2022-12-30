

olsfc_BPI_F <- function(x, h, maxlag = 0, nolag = NULL , sim = 2000) {
  # Set up data frame for modelling
  n <- length(x)
  modeldata <- data.frame(
    x = as.numeric(x),
    trend = seq_along(x),
    f1s_day = sin(2*pi*seq_along(x)/(24)),
    f1c_day = cos(2*pi*seq_along(x)/(24)),
    f2s_day = sin(4*pi*seq_along(x)/(24)),
    f2c_day = cos(4*pi*seq_along(x)/(24)),
    f3s_day = sin(6*pi*seq_along(x)/(24)),
    f3c_day = cos(6*pi*seq_along(x)/(24)),
    f1s_week= sin(2*pi*seq_along(x)/(7*24)),
    f1c_week = cos(2*pi*seq_along(x)/(7*24)),
    f2s_week = sin(4*pi*seq_along(x)/(7*24)),
    f2c_week = cos(4*pi*seq_along(x)/(7*24))
  )
  if(maxlag > 0)
  {
    lagnames <- paste0("lag", seq(maxlag))
    for (i in seq_along(lagnames))
      modeldata[[lagnames[i]]] <- c(rep(NA, i), x[seq(n - i)])
  } else {
    lagnames <- NULL
  }
  # Set up formula for linear model
  if (length(nolag) == 0)
    nolag <- seq(maxlag) 
  form <- "x ~ trend + f1s_day + f1c_day + f2s_day + f2c_day + f3s_day + f3c_day + f1s_week + f1c_week + f2s_week + f2c_week"
  for (i in nolag)
    form <- paste0(form, " + ", "lag", i)
  form <- as.formula(form)
  
  # Fit model
  fit <- lm(form, data = modeldata)
  # Adjust the residuals
  s.resid1 <- residuals(fit)
  leverage <- influence(fit)$hat
  s.resid <- residuals(fit)/sqrt(1-leverage)
  s.resid <- s.resid - mean(s.resid)
  # Set up forecast vector
  fc <- msts(numeric(h),
             seasonal.periods=c(24,24*7),
             start = c(2, 97)
  )
  lwr <- msts(numeric(h),
              seasonal.periods=c(24,24*7),
              start = c(2, 97)
  )
  upr <- msts(numeric(h),
              seasonal.periods=c(24,24*7),
              start = c(2, 97)
  )
  # Set up new data for forecasting
  trend <- length(x) + seq(h)
  f1s_day = sin(2*pi*(length(x) + seq(h))/(24))
  f1c_day = cos(2*pi*(length(x) + seq(h))/(24))
  f2s_day = sin(4*pi*(length(x) + seq(h))/(24))
  f2c_day = cos(4*pi*(length(x) + seq(h))/(24))
  f3s_day = sin(6*pi*(length(x) + seq(h))/(24))
  f3c_day = cos(6*pi*(length(x) + seq(h))/(24))
  f1s_week = sin(2*pi*(length(x) + seq(h))/(7*24))
  f1c_week = cos(2*pi*(length(x) + seq(h))/(7*24))
  f2s_week = sin(4*pi*(length(x) + seq(h))/(7*24))
  f2c_week = cos(4*pi*(length(x) + seq(h))/(7*24))
  newdata <- cbind.data.frame(trend = trend[1], f1s_day = f1s_day[1], f1c_day = f1c_day[1], f2s_day = f2s_day[1], 
                              f2c_day = f2c_day[1], f3s_day = f3s_day[1], f3c_day = f3c_day[1], f1s_week = f1s_week[1], 
                              f1c_week = f1c_week[1], f2s_week = f2s_week[1], f2c_week = f2c_week[1])
  for (i in seq_along(lagnames))
    newdata[[lagnames[i]]] <- tail(x,i)[1]
  newdata.frame <- matrix(NA, ncol = ncol(newdata), nrow = h)
  ep.drawslist <- list()
  # Compute forecasts
  for (i in seq_along(fc))
  {
    fore <- predict(fit , newdata = newdata)
    ## PI using bootstrap
    #ep.draws <- replicate( n = sim, unname(fore + sample(s.resid,size=1)))
    ep.draws <- fore + sample(s.resid,size=sim, replace = TRUE)
    PI <- fore + quantile(ep.draws, probs = c(0.025, 0.975))
    fc[i] <- unname(fore)
    lwr[i] <- unname(PI[1])
    upr[i] <- unname(PI[2])
    ep.drawslist[[length(ep.drawslist)+1]] <- ep.draws
    # Update newdata
    if(maxlag > 0)
    {
      newdata[lagnames[seq(maxlag)]] <- c(fc[i], newdata[lagnames[1:(maxlag-1)]]) 
    }
    newdata[['trend']] <- trend[i+1]
    if(i == max(seq_along(fc)))
    {
      newdata[['f1s_day']] <- 'NA'
      newdata[['f1c_day']] <- 'NA'
      newdata[['f2s_day']] <- 'NA'
      newdata[['f2c_day']] <- 'NA'
      newdata[['f3s_day']] <- 'NA'
      newdata[['f3c_day']] <- 'NA'
      newdata[['f1s_week']] <- 'NA'
      newdata[['f1c_week']] <- 'NA'
      newdata[['f2s_week']] <- 'NA'
      newdata[['f2c_week']] <- 'NA'
    }
    else{
      newdata[['f1s_day']] <- f1s_day[i+1]
      newdata[['f1c_day']] <- f1c_day[i+1]
      newdata[['f2s_day']] <- f2s_day[i+1]
      newdata[['f2c_day']] <- f2c_day[i+1]
      newdata[['f3s_day']] <- f3s_day[i+1]
      newdata[['f3c_day']] <- f3c_day[i+1]
      newdata[['f1s_week']] <- f1s_week[i+1]
      newdata[['f1c_week']] <- f1c_week[i+1]
      newdata[['f2s_week']] <- f2s_week[i+1]
      newdata[['f2c_week']] <- f2c_week[i+1]
    }
  }
  return(list(cbind.data.frame('fc' = fc, 'lwr' = lwr, 'upr' = upr), s.resid1, ep.drawslist))
}


