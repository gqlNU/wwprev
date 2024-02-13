# Integrating wastewater and randomised prevalence survey data for national COVID surveillance

A statistical framework to nowcast prevalence of COVID-19 at a fine spatio-temporal resolution by integrating fine space-time scale wastewater measurements, an affordable non-traditional disease metric for COVID-19, with prevalence data available at a coarse-spatial resolution collected from prevalence surveys run at a reduced scale.

<img src="https://github.com/gqlNU/wwprev_dev/assets/6213918/5924a733-3a67-4d0e-b106-94c2dede3e3a" width="600">


##  Installing the package

You first need to obtain a personal access token (PAT) by following the instruction via this link:  https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/creating-a-personal-access-token.

Enter the token as a string to the object `PAT` then run the following code to install the package:

```R
PAT <- " " # enter your PAT token here
devtools::install_github("gqlNU/wwprev", auth_token = PAT)
```

You also need to install the `R` package `publicWW` [] using the following code:

```R
devtools::install_github("gqlNU/publicWW", auth_token = PAT) # the same PAT used above
```


##  Steps in obtaining the nowcast prevalence at fine spatial scale

### 1. Model fitting.
This script [inst/scripts/fitting.R](inst/scripts/fitting.R) fits a model from the data integration framework to the weekly prevalence estimates and the weekly wastewater viral concentration estimates. The default setting simulates the scenario where the two sets of estimates are available at both the Lower Tier Local Authority (LTLA) level and the national level for the first 20 weeks and the weekly prevalence estimates are only available at the national level for the subsequent 20 weeks. The aim of the data integration framework is to disaggregate the national level prevalence to the LTLA level for the latter 20 weeks, i.e. the nowcast period.

### 2. Disaggregation
Taking the output from Step 1, the script [inst/scripts/disaggregating.R](inst/scripts/disaggregating.R) produces the nowcast prevalence at the LTLA level.

### 3. Summary of nowcast
The script [inst/scripts/summarising.R](inst/scripts/summarising.R) assesses the quality of the nowcast prevalence against the observed but held-out LTLA prevalence.


##  Scripts used to produce some of the results in the manuscript

- [inst/scripts/tau_map_onGitHub.R](inst/scripts/tau_map_onGitHub.R) produces the map of correlation between prevalence and wastewater viral concentration across the LTLAs in England. The LTLA and the regions boundaries are available within this package and were obtained via [https://geoportal.statistics.gov.uk/datasets/77c428630a664bf08301a0d07a2adbc9/] and [https://geoportal.statistics.gov.uk/datasets/35ec07df9b774e938dd5b5595f10dfcc/] respectively (last accessed 13/02/2024).

- [inst/scripts/ww_only_data_integration.R](inst/scripts/ww_only_data_integration.R) compares the LTLA-weekly nowcast prevalence from two models, the data integration model that uses both the wastewater data and national-level prevalence and one that only uses the wastewater data.  


##  Data sources used in the study
- Weekly-LTLA level debiased prevalence estimates from Nicholson et al. 2022 [https://www.nature.com/articles/s41564-021-01029-0] (included in the `wwprev` package)
- Weekly-LTLA and weekly-national estimates of wastewater viral concentration from Li et al. (2023) [https://www.sciencedirect.com/science/article/pii/S0160412023000387] (included in the `wwprev` package)
- English Indices of Multiple Deprivation 2019: [https://www.gov.uk/government/statistics/english-indices-of-deprivation-2019]
- Ethnicity based on the 2021 Census: [https://www.ons.gov.uk/visualisations/dvc2203/map/datadownload.xlsx]

<!--
(the data) used for Figure 3 shown in [https://www.ons.gov.uk/peoplepopulationandcommunity/culturalidentity/ethnicity/bulletins/ethnicgroupenglandandwales/census2021](this report)
-->
