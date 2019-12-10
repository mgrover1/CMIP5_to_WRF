# CMIP5_to_WRF

### Calculates Climate Change Deltas Following the methodology from Trapp and Hoogewind (2016)

![Figure1_Trapp_Hoogewind](/Figures/Equations/Climo_Delta_Equation.PNG)

![Figure2_Trapp_Hoogewind](/Figures/Equations/Climo_Delta_Equation2.PNG)

### CMIP5 Models Used - Downloaded from [CMIP5 Data Archive](https://esgf-node.llnl.gov/projects/cmip5/)

- MIROC5
- GFDL
- NCAR CCSM4

### Overall Workflow
- Calculate Climate Change Deltas
    - Notebooks Used
        - 3D Atmos
        - 2D Atmos
        - Surface

- Add to NAM fields
    - Notebook Used
        - Write to Grib
 
- Export to GRIB file to be used by WRF
    - Notebook Used
        - Write to Grib

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
