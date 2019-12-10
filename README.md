# CMIP5_to_WRF

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

### Overall Workflow
- Calculate Climate Change Deltas
    - Notebooks Used
        - [3D Atmos](/Notebooks/3D-Vars.ipynb)   
        - [2D Atmos](/Notebooks/2D-Vars.ipynb)
        - [Surface](/Notebooks/Surface_Variables.ipynb)

- Add to NAM fields
    - Notebook Used
        - [Write to Grib](/Notebooks/Write_to_Grib.ipynb)
 
- Export to GRIB file to be used by WRF
    - Notebook Used
        - [Write to Grib](/Notebooks/Write_to_Grib.ipynb)

### Climate Change Delta Variables
- 3D Atmos
    - UA (U-Component of Wind at Isobaric Levels)
    - VA (V-Component of Wind at Isobaric Levels)
    - TA (Temperature of Atmosphere at Isobaric Levels)
    - ZG (Geopotential Heights at Isobaric Levels)
    - HUR (Relative Humidity at Isobaric Levels)
    
- 2-D
    - UAS (U-Component of Wind at 10m Above Surface)
    - VAS (V-Component of Wind at 10m Above Surface)
    - TAS (Temperature of Wind at 2m Above Surface)
    - HURS (Relative Humidity at 2m Above Surface)
    - HUSS (Specific Humidity at 2m Above Surface)
    - PS (Pressure at Surface)
    - PSL (Mean Sea-Level Pressure at Surface)
    - TS (Temperature at Surface (Skin Temperature))

- Soil/Surface
    - TSL (Temperature at Soil Levels) - Land Mask Added
    - MRLSL (Specific Humidity at Soil Levels) - Land Mask Added
