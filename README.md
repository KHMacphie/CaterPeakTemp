# CaterPeakTemp
Data and code for "Warmer springs lead to earlier and higher peaks of arboreal caterpillars"

SlidingWindow.R (Code)
Code to model the phenological distribution of caterpillars for each site by year combination (site-year) and use parameter estimates in a multivariate meta-analysis sliding window to identify the period when temperature best predicts each phenological parameter.

SpatioTemporalModel.R (Code)
Stan code to run main spatio-temporal model and code to produce results and figures

SpaceVsTimeModel.R (Code)
Data manipulation to calculate site mean temperatures and annual deviations, stan code to run space versus time model and code to produce results

data_cater.csv (Data)
Branch beating data with site-year specific temperatures for each phenological parameter

sy_daily_temp.csv (Data)
Daily mean temperatures for each site-year for use in the sliding window

sy_temp_by_para.csv (Data)
Site-year mean temperatures during the windows identified for each phenological parameter in the sliding window 
