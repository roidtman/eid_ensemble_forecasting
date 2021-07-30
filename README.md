# Trade-offs between individual and ensemble forecasts of an emerging infectious disease (Oidtman *et al.*)


## Getting started

The R scripts are written assuming you have the following folder structure:

```
eid_ensemble_forecasting
│   README.md
└─── code
└─── output
└─── data
```
In the code/data_assimilation_script_example folder, there are example R Markdown files to run the data assimilation algorithm, which produces 20,000 forecasts for 15 time points. For the entirety of the analysis, there were 12 spatially coupled models for which 20,000 forecasts were generated for each of the 15 time points and there were 4 non-spatial models for which 20,000 forecasts were generated for each of the 15 time points for each of the 31 departments considered in the analysis. At each time point, a forecast was generated for the rest of the time series, so the length of each forecast became shorter as more data was assimilated into the model. In addition to the forecasts, model fits were also generated at each assimilation period. Model fits were only for the four week periods of each assimilation period.


There are 16 models which we produce forecasts for, which we internally number 1-16. 


Models 1-2, 5-6, 9-10, 13-14 all use a static R.
Models 3-4, 7-8, 11-12, 15-16 all use a dynamics R.
Models 5-8 are non-spatial, while 1-4 and 9-16 are all spatially coupled.
Models 1-4 use a CDR-informed mobility matrix, models 9-12 use a gravity model, models 13-16 use a radiation model.
The odd numbered models assume 1 introduction and the even numbered models assume 2 introductions.


In the code/data_assimilation_script_example folder, we have example scripts for models that were run on the Databricks computing platform. For models 9-16, we included each script. For the non-spatial models, we included one example for La Guajira (i.e., dept 1). We did not include models that make use of the CDR-informed mobility matrix (i.e., models 1-4). The mobile phone data set used in this study is proprietary and subject to strict privacy regulations. The access to this data set was granted after reaching a non-disclosure agreement with the proprietor, who anonymized and aggregated the original data before giving access to the authors. The mobile phone is available on request after negotiation of a non-disclosure agreement with the company. The contact person is Enrique Frías-Martínez (enrique.friasmartinez@telefonica.com). 


All output csv files from those forecasting algorithms are included in the output folder, with numbers on folders corresponding to model number. Within each no_tm_ini_# folder, there is the I_A folder (model fits), I_F folder (model forecasts), particle folder (particles that generate the I_A and I_F) and RData files for the processed forecasts and model fits. Those that are indexed with rho denote they correspond to cases, while those without rho denote they correspond to infections. 

### Software and packages

We used R version 3.6.3, "Holding the Windsock" to process  
We used R on Databricks to run forecasting algorithms. 

R packages (with version number) necessary for these analyses:
* viridis (0.5.1)
* lubridate (1.7.4)
* sn (1.6.1)
* dplyr (0.8.5)
* stats (3.6.3)
* yarrr (0.1.5)
* vioplot (0.3.4)
* sf (0.8.1)
* maptools (0.9.9)
* fields (10.3)
* stringr (1.4.0)
* rgeos (0.5.2)
* rgdal (1.4.8)
* mvtnorm (1.1.0)

## Processing the forecast output
The file code/processing/process_forecasts.R processes forecast output for each spatially coupled model and produces the I_A_processed.RData, I_F_I_A_processed_rho.RData, and I_F_processed.RData for each respective model. The file code/processing/process_non-spatial_forecasts.R processes forecast output for each department for the non-spatial models and produces one RData file for each model encompassing all departments. It is necessary to first process forecast outputs into RData files for subsequent analyses. 

**For spatially coupled models**: To run the processing_forecasts.R, you must identify the model output (and file path) you are hoping to process. You can do this directly by changing the file_path argument on line 17. Path names correspond to the folder names in the output folder. Spatial models are 1-4, 9-16. 

**For non-spatial models**: To run the processing_non-spatial_forecasts.R, you must identify the model output (and file path) you are hoping to process. You can do this directly by changing the file_path argument on line 17. Path names correspond to the folder names in the output folder. Non-spatial models are 5-8. 

Warning: the processed RData output files for these models are very large. 

## Producing ensemble weights and other analyses
* The file code/1_analysis_forecasting_ensemble_through_time_EM.R produces the EM ensemble weights. 
* The file code/1_analysis_forecasting_ensemble_EW.R produces the equally-weighted ensemble. 
* The file code/1_analysis_forecasting_output_per_model.R provides sanity checks and exploratory data analysis plots. 
* The file code/1_analysis_particle_retention_rate.R provides analyses on particle retention and explores *why* forecasts are improving.

## Producing figures
All other code/produce_figure files provide analyses and processing necessary for figures. 



## Questions about code and software can be directed to Rachel Oidtman (rjoidtman@gmail.com)

