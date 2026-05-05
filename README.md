# Bhutan climate downscaling (250 m)

This repository contains scripts used to generate a high-resolution (250 m) climate dataset for Bhutan based on CMIP6 general circulation models (GCMs).

## Overview

The dataset includes monthly climate variables for both historical (1980–2014) and future (2015–2100) periods under four CMIP6 scenarios (SSP126, SSP245, SSP370, SSP585). The downscaling was performed using a delta (change factor) approach, combining GCM-derived climate change signals with observation-based climatology.

## Variables

### Primary variables
- Precipitation (pr)
- Minimum temperature (tasmin)
- Maximum temperature (tasmax)
- Vapour pressure (vp)

### Derived variables
- Bioclimatic variables (BIO1–BIO19)

Vapour pressure was derived from CMIP6 specific humidity (huss) and surface pressure (ps), aggregated to monthly means, and subsequently downscaled to 250 m resolution. It is included to support derivation of humidity-related variables (e.g., relative humidity) required for process-based models such as CLIMEX.

## Data Structure

### 1. Monthly climate variables
The primary variables (pr, tasmin, tasmax, vp) are organised as:

```
variable/
  experiment/
    model/
```

Example:

```
pr/
  ssp126/
    ACCESS-CM2/
```

Each file contains monthly data for a single year:

```
<variable>_<model>_<experiment>_<period>_<year>.tif
```

Example:

```
pr_ACCESS-CM2_ssp126_2015-2100_2015.tif
```

Each GeoTIFF file contains 12 layers corresponding to monthly values.

### 2. Bioclimatic variables
Bioclimatic variables (BIO1–BIO19) are provided as derived products.

**Historical:**

```
bioclim/
  historical/
    <model>/
      bio1.tif ... bio19.tif
```

**Future:**

```
bioclim/
  future/
    <scenario>/
      <period>/
        <model>/
          bio1.tif ... bio19.tif
```

These are generated for each model, scenario, and time period.

## Methods (Summary)
- Daily CMIP6 data were aggregated to monthly values  
- Historical climatologies were computed for each model  
- Change factors were calculated (multiplicative for precipitation and appropriate scaling for temperature and related variables)  
- Change factors were applied to observation-based climatology  
- Outputs were downscaled to 250 m resolution  

## Vapour Pressure (Important Note)
Vapour pressure is available for one GCM (NorESM2-LM), as it was the only model with complete specific humidity (huss) and surface pressure (ps) data across all scenarios. This variable was generated to support derivation of relative humidity and related applications.

## Requirements
The scripts use the following R packages:

- terra
- stringr
- base R functions

## Usage
The scripts are provided as reference implementations of the workflow used to generate the dataset.  
Users will need to adapt file paths and input data locations to their own systems.

## Data Availability
The dataset is available at:  
[INSERT YOUR DATA REPOSITORY LINK HERE]

---

## Citation
If you use this dataset or code, please cite:

Dorji, S. (Year). *A high-resolution climate dataset for Bhutan*. [Journal / DOI]

---

## Contact
Sangay Dorji
Email: [dorjismo@gmail.com]
