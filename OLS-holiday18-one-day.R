library(tidyverse)
library(forecast)
source('smatrix.R')

## reading data
hour.all <- read.csv('hour.all.csv', header = TRUE)[,-1]

## Drop highway No.5
hour.all <- subset(hour.all, Freeway!= 'No5')

## Coding station, freeway and area
hour.all$Station.code <- as.numeric(as.factor(hour.all$Station))
hour.all$Freeway.code <- as.numeric(as.factor(hour.all$Freeway))
hour.all$Area.code <- as.numeric(as.factor(hour.all$Area))


## fixing the station codes
hour.all$Station.code <- sprintf("%03d", as.numeric(hour.all$Station.code))
hour.all$vehicle.type <- sprintf("%02d", as.numeric(hour.all$vehicle.type))

## Creating category column
hour.all$cat <- paste(hour.all$sensor, hour.all$direction, hour.all$vehicle.type,  hour.all$Station, hour.all$Freeway, hour.all$Area, hour.all$Region)

### removing just zero series
hour.all <- hour.all %>%
  group_by(cat) %>%
  filter(any(sum(traffic != 0)))
#
## removing 80% zero series
hour.all <- hour.all %>%
  group_by(cat) %>%
  filter(mean(traffic == 0) <= 0.8)

## sorting all the series
hour.all <- hour.all[order(hour.all$cat),]


## Each series in one column
hourmatrix <- hour.all$traffic %>%
  matrix(nrow = 20424, ncol = 1590) %>%
  as.data.frame() 

## Create columns' name
names <- unique(paste(substr(hour.all$Region, 1, 1), hour.all$Station.code,  hour.all$Freeway.code, hour.all$direction, hour.all$vehicle.type, sep = "")) 
colnames(hourmatrix) <- names

## Hierarchy structure
name_length <- str_length(names)
grouping_hts <- rbind(
  #Region
  str_sub(names, start = name_length - 7, end = name_length - 7),
  #Station
  str_sub(names, start = name_length - 7, end = name_length - 4),
  #Freeway
  str_sub(names, start = name_length - 3, end = name_length - 3),
  #Direction
  str_sub(names, start = name_length - 2, end = name_length - 2),
  #Vehicle
  str_sub(names, start = name_length - 1, end = name_length),
  #Region x Direction
  paste(str_sub(names, start = name_length - 7, end = name_length - 7), str_sub(names, start = name_length - 2, end = name_length - 2), sep = ""),
  #Region x Freeway
  paste(str_sub(names, start = name_length - 7, end = name_length - 7), str_sub(names, start = name_length - 3, end = name_length - 3), sep = ""),
  #Region x Vehicle
  paste(str_sub(names, start = name_length - 7, end = name_length - 7), str_sub(names, start = name_length - 1, end = name_length), sep = ""),
  #Freeway x Direction
  paste(str_sub(names, start = name_length - 3, end = name_length - 3), str_sub(names, start = name_length - 2, end = name_length - 2), sep = ""),
  #Freeway x Vehicle
  paste(str_sub(names, start = name_length - 3, end = name_length - 3), str_sub(names, start = name_length - 1, end = name_length), sep = ""),
  #Vehicle x Direction
  paste(str_sub(names, start = name_length - 1, end = name_length), str_sub(names, start = name_length - 2, end = name_length - 2), sep = "")
)

library(hts)
hourgts <- gts(hourmatrix, groups = grouping_hts) 

allyhourall <- aggts(hourgts) 

## Only the last set of long holidays - One-day
holiday.18 <- allyhourall[20065:20424,]

allyholiday.18 <- forecast::msts(holiday.18, seasonal.periods=c(24,24*7))



## base forecasts & 2000 sample paths 
n <- 360
t <- 336
h <- 24
result.all <- list()
sample.path <- list()
error.train.all <- NULL
start.time <- Sys.time()
for(j in 1:ncol(allyholiday.18)){
  result <- NULL
  sample.test <- NULL
  error.train <- NULL
  for(i in 0:((n-t-h)/h)){
    train <- window(allyholiday.18[,j], start = c(1, ((i*h)+1)), end = c(1, (t + (i*h))))
    valid <- window(allyholiday.18[,j], start = c(1,( t + (i*h) + 1)), end = c(1, (t + (i*h) +h)))
    set.seed(123)
    f <- olsfc_BPI_F(train, h, maxlag = 24, nolag = c(1,24))
    sample.test1 <- f[[3]]
    sample.test2 <- do.call(rbind, sample.test1)
    sample.test <- rbind(sample.test, sample.test2)
    m <- data.frame(f[[1]],apply(f[[1]],2,as.numeric))[,4:6]
    result <- rbind(result, m)
    error.train <- bind_rows(error.train, as.data.frame(f[[2]]))
    print(j)
  }
  result.all[[length(result.all)+1]] <- result
  sample.path[[length(sample.path)+1]] <- sample.test
  error.train.all <- bind_cols(error.train.all, error.train)
}
end.time <- Sys.time()
time.taken.base <- end.time - start.time


#### saving results
fc.OLS <- NULL
for(i in 1:length(result.all)){
  fc.OLS <- cbind(fc.OLS, result.all[[i]][,1])
}

## ###
## reconciling base forecasts and sample paths
## ###

### reading data in a list (each dataframe is for one series)

library(rio)

mydataf <- lapply(1:ncol(sample.path[[1]]), function(y) as.data.frame(sapply(sample.path, function(x) x[, y])))


# name the data frames
names(mydataf) <- 1:length(mydataf)

fc.OLS.res <- list()
for(i in 0:0){
  fc.OLS.res[[i+1]] <- error.train.all[(1+(312*i)):(312+(312*i)),]
}

library(Matrix)
start.time <- Sys.time()
list.G <- list()
for(i in 1:length(fc.OLS.res)){
  res <- as.matrix(fc.OLS.res[[i]])
  n <- nrow(res)
  covm <- crossprod(stats::na.omit(res)) / n
  tar <- diag(apply(res, 2, compose(crossprod, stats::na.omit))/n)
  corm <- cov2cor(covm)
  xs <- scale(res, center = FALSE, scale = sqrt(diag(covm)))
  xs <- xs[stats::complete.cases(xs),]
  v <- (1/(n * (n - 1))) * (crossprod(xs^2) - 1/n * (crossprod(xs))^2)
  diag(v) <- 0
  corapn <- cov2cor(tar)
  d <- (corm - corapn)^2
  lambda <- sum(v)/sum(d)
  lambda <- max(min(lambda, 1), 0)
  W <- lambda * tar + (1 - lambda) * covm
  gmat <- GmatrixG(hourgts$groups)
  smatrix <- as.matrix(SmatrixM(gmat))
  R <- t(smatrix)%*%solve(W)
  P <- Matrix::solve(R%*%smatrix)%*%R
  SP <- smatrix%*%P
  list.G[[length(list.G)+1]] <- SP
  print(i)
}
## Redonciling sample paths
list.sample.rec <- list()
for(i in 1:length(mydataf)){
  fit.test <- split(mydataf[[i]], rep(1:1,each=(24)))
  result.path <- NULL
  for(j in 1:length(list.G)){
    fc.rec <- matrix(NA, nrow = 24, ncol = ncol(fc.OLS.res[[1]]))
    for(k in 1:nrow(fit.test[[1]])){
      f.1 <- matrix(as.numeric(fit.test[[j]][k,]), ncol = 1, nrow = ncol(fit.test[[1]]))
      fc.rec [k,] <- list.G[[j]] %*% f.1
    }
    result.path <- bind_rows(result.path, as.data.frame(fc.rec))
  }
  list.sample.rec[[length(list.sample.rec)+1]] <- result.path
  print(i)
}

end.time <- Sys.time()
time.taken.sample.path.rec <- end.time - start.time

# name the data frames
names(list.sample.rec) <- 1:length(list.sample.rec)

## Computing prediction intervals - 95% and 99%
mydataf2 <- lapply(1:ncol(list.sample.rec[[1]]), function(y) as.data.frame(sapply(list.sample.rec, function(x) x[, y])))

quan975 <- matrix(NA, ncol = length(mydataf2), nrow = nrow(mydataf2[[1]]))

for(i in 1:length(mydataf2)){
  quan975[,i] <-  apply(mydataf2[[i]], 1, quantile, probs= c(0.975))
}

quan25 <- matrix(NA, ncol = length(mydataf2), nrow = nrow(mydataf2[[1]]))

for(i in 1:length(mydataf2)){
  quan25[,i] <-  apply(mydataf2[[i]], 1, quantile, probs= c(0.025))
}

quan95 <- matrix(NA, ncol = length(mydataf2), nrow = nrow(mydataf2[[1]]))

for(i in 1:length(mydataf2)){
  quan95[,i] <-  apply(mydataf2[[i]], 1, quantile, probs= c(0.95))
}

quan05 <- matrix(NA, ncol = length(mydataf2), nrow = nrow(mydataf2[[1]]))

for(i in 1:length(mydataf2)){
  quan05[,i] <-  apply(mydataf2[[i]], 1, quantile, probs= c(0.05))
}






## Reconciling base forecasts
colnames(fc.OLS) <- colnames(allyhourall)

fit.list <- list()
for(i in 0:0){
  fit.list[[i+1]] <- fc.OLS[(1+(24*i)):(24+(24*i)),]
}
fc.OLS.res <- list()
for(i in 0:0){
  fc.OLS.res[[i+1]] <- error.train.all[(1+(312*i)):(312+(312*i)),]
}


list <- list()
list.G <- list()
start.time <- Sys.time()
for(i in 1:length(fc.OLS.res)){
  res <- as.matrix(fc.OLS.res[[i]])
  n <- nrow(res)
  covm <- crossprod(stats::na.omit(res)) / n
  tar <- diag(apply(res, 2, compose(crossprod, stats::na.omit))/n)
  corm <- cov2cor(covm)
  xs <- scale(res, center = FALSE, scale = sqrt(diag(covm)))
  xs <- xs[stats::complete.cases(xs),]
  v <- (1/(n * (n - 1))) * (crossprod(xs^2) - 1/n * (crossprod(xs))^2)
  diag(v) <- 0
  corapn <- cov2cor(tar)
  d <- (corm - corapn)^2
  lambda <- sum(v)/sum(d)
  lambda <- max(min(lambda, 1), 0)
  W <- lambda * tar + (1 - lambda) * covm
  gmat <- GmatrixG(hourgts$groups)
  smatrix <- as.matrix(SmatrixM(gmat))
  R <- t(smatrix)%*%solve(W)
  P <- Matrix::solve(R%*%smatrix)%*%R
  SP <- smatrix%*%P
  
  fit.test <- fit.list[[i]]
  fc.rec <- matrix(NA, nrow = 24, ncol = ncol(allyhourall))
  for(j in 1:nrow(fit.test)){
    f.1 <- matrix(as.numeric(fit.test[j,]), ncol = 1, nrow = ncol(fit.test))
    fc.rec [j,] <- SP %*% f.1
  }
  colnames(fc.rec) <- colnames(allyhourall) 
  list[[length(list)+1]] <- fc.rec
  list.G[[length(list.G)+1]] <- P
  print(i)
}

end.time <- Sys.time()
time.taken.forecast.rec <- end.time - start.time

result.fit.all.rec <- do.call(rbind.data.frame, list)


### print computation times
time.taken.base
time.taken.sample.path.rec
time.taken.forecast.rec

