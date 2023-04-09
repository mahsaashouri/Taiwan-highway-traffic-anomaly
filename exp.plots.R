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


week.data <- data.frame(
  date = seq(as.POSIXct("2021-04-23 23:00:00", tz="CET"), as.POSIXct("2021-04-30 22:00:00", tz="CET"), by="hour"),
  'data' = tail(allyhourall[,'N2951S31'], 168))



ggplot(data = week.data) + geom_line(aes(x = date, y = data)) +
  ylab('Number of cars')+xlab('Date')  +
  scale_x_datetime( breaks = scales::date_breaks("24 hour"),
                    labels = scales::date_format("%a-%d\n%H:%M",  tz="CET")) +
  theme_light() +
  theme(
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    legend.title=element_text(size=20), 
    legend.text=element_text(size=20))


year.data <- data.frame(
  date = seq(as.POSIXct("2021-01-10 23:00:00", tz="CET"), as.POSIXct("2021-04-30 22:00:00", tz="CET"), by="hour"),
  'data' = tail(allyhourall[,'S2693N42'], 2639))

ggplot(data = year.data) + geom_line(aes(x = date, y = data)) +
  ylab('Number of big trucks')+xlab('Date')  +
  theme_light() +
  theme(
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    legend.title=element_text(size=20), 
    legend.text=element_text(size=20))

## ACF and PACF plots

ggAcf(year.data$data, main='') + 
  theme(
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    legend.title=element_text(size=20), 
    legend.text=element_text(size=20))

ggPacf(year.data$data , main='') + 
  theme(
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    legend.title=element_text(size=20), 
    legend.text=element_text(size=20))

