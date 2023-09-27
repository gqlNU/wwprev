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




<!--
##  Data sources
- Weekly-LTLA level prevalence estimates from Nicholson et al. 2022
- English Indices of Multiple Deprivation 2019: [https://www.gov.uk/government/statistics/english-indices-of-deprivation-2019]
- Ethnicity based on the 2021 Census: (the data)[https://www.ons.gov.uk/visualisations/dvc2203/map/datadownload.xlsx] used for Figure 3 shown in (this report)[https://www.ons.gov.uk/peoplepopulationandcommunity/culturalidentity/ethnicity/bulletins/ethnicgroupenglandandwales/census2021]

- English Indices of Multiple Deprivation 2019: [https://www.gov.uk/government/statistics/english-indices-of-deprivation-2019]
- Ethnicity based on the 2021 Census: [https://www.ons.gov.uk/peoplepopulationandcommunity/culturalidentity/ethnicity/bulletins/ethnicgroupenglandandwales/census2021#how-ethnic-composition-varied-across-england-and-wales]
-->
