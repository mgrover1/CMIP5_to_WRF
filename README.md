# Dynamical Downscaling Using Pseudo Global Warming Method: A Case Study"

### Calculates Climate Change Deltas Following the methodology from [Trapp and Hoogewind (2016)](https://journals.ametsoc.org/doi/full/10.1175/JCLI-D-15-0623.1)

![Figure1_Trapp_Hoogewind](/Figures/Equations/Climo_Delta_Equation.PNG)
Figure 1 from Trapp and Hoogewind (2016)

![Figure2_Trapp_Hoogewind](/Figures/Equations/Climo_Delta_Equation2.PNG)
Figure 2 from Trapp and Hoogewind (2016)

![Summary Climate Delta](/Figures/Equations/Summary_Climate_Delta.png)
Summary of Climate Delta Calculation - this process is followed for all variables listed below

### CMIP5 Models Used - Downloaded from [CMIP5 Data Archive](https://esgf-node.llnl.gov/projects/cmip5/)

- MIROC5
- GFDL
- NCAR CCSM4

All hourly data are downloaded and stored in the following directory on Keeling ```/data/meso/a/mgrover4/```
with directories using naming convention ```model_6hour```

Monthly data are stored in ```/data/jtrapp/a/jtrapp/CMIP5``` in directories MIROC, GFDL, or NCAR

### Overall Workflow

- (Optional) Download CMIP data from esgf
    - See shell scripts within the data directories
    - Data for May and March are already included

- Calculate Climate Change Deltas
    - Notebooks Used
        - [Atmos Deltas](/Notebooks/Atmos_Deltas.ipynb)
        - [Soil Deltas](/Notebooks/Soil_Deltas.ipynb)
    - Notebooks can be found on Keeling in ```/data/keeling/a/mgrover4/b/CMIP5_to_WRF/Notebooks/```

- Add Climate Deltas to NAM fields
    - Notebook Used
        - [Write to Grib](/Notebooks/Write_to_Grib.ipynb)

- Export to GRIB file to be used by WRF
    - Notebook Used
        - [Write to Grib](/Notebooks/Write_to_Grib.ipynb)

### Climate Change Delta Variables
- 3D Atmos (6 hourly, linearly interpolated to 3 hourly)
    - UA (U-Component of Wind at Isobaric Levels)
    - VA (V-Component of Wind at Isobaric Levels)
    - TA (Temperature of Atmosphere at Isobaric Levels)
    - HUS (Specific Humidity at Isobaric Levels)

- 2-D (Most data 3 hourly, some 6 hourly data were interpolated to 3 hourly)
    - UAS (U-Component of Wind at 10m Above Surface)
    - VAS (V-Component of Wind at 10m Above Surface)
    - TAS (Temperature of Wind at 2m Above Surface)
    - HUSS (Specific Humidity at 2m Above Surface)
    - PS (Pressure at Surface)
    - PSL (Mean Sea-Level Pressure at Surface)
    - TS (Temperature at Surface (Skin Temperature))

- Soil/Surface (Monthly data)
    - TSL (Temperature at Soil Levels) - Land Mask Added
    - MRLSL (Specific Humidity at Soil Levels) - Land Mask Added
