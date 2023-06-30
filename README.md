#  Government policy evaluation on consecutive holidays: traffic flow of Taiwanese highways

This study examines the Taiwanese highway traffic behavior on consecutive holidays. 
This investigation is essential for evaluating the Taiwanese freeway bureau's traffic control and management strategies. 
For this reason, a prediction-based detection method is proposed for finding highway traffic anomalies using reconciled ordinary least squares (OLS) 
forecasts and bootstrapped prediction intervals. Two fundamental features of traffic flow time series, seasonality, and spatial autocorrelation, are captured by adding 
Fourier terms in OLS models, spatial aggregation, and reconciliation step. Our OLS approach is a simple method that models complex traffic datasets with reasonable 
accuracy while avoiding the computational problem. Our dataset includes hourly time series of Taiwan's consecutive holidays in 2019, 2020, 
and 2021 (73 days - 1752 hours). Our time series is the number of vehicles collected by the Electronic Tolling Collection (ETC) system on three main highways. 

The following figure shows an examples of Taiwanese highway hourly time series (from 2021-01-10 23:00:00 to 2021-04-30 22:00:00) 
in northern region, Zhongli runway - Hukou, south direction, and car as vehicle types:

![alt text](<https://github.com/mahsaashouri/Taiwan-highway-traffic-anomaly/blob/main/N2971S31.png>)
 
Analyses of Taiwan's consecutive holidays showed strong variations in anomalies across different directions and highways. Specifically, we detected some areas and highways comprising a high number of traffic anomalies (north direction-central and southern regions-highways No. 1 and 3, south direction-southern region-highway No.3), and others with generally normal traffic (east and west direction). These results could provide important decision-support information to traffic authorities.


### Note:
For a complete dataset, please contact: mah.ashouri62@gmail.com
