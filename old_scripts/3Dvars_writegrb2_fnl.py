#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')
import numpy as np
from mpl_toolkits.basemap import Basemap
from ncepgrib2 import Grib2Encode, Grib2Decode
import pygrib
import os
import Nio

# Note:  this works on conte, but not on coates (the ncepgrib2 doesn't exist?)
# usage:  ./RH_writegrb2.py

hourlp = ['00','06','12','18','00','06','12','18','00','06','12','18','00','06','12','18','00','06','12','18','00']
day0lp = ['14','14','14','14','15','15','15','15','16','16','16','16','17','17','17','17','18','18','18','18','19']
daylp = ['14','14','14','14','15','15','15','15','16','16','16','16','17','17','17','17','18','18','18','18','19']

for nt in range(len(daylp)):


#loop through...
 day = daylp[nt]
 hour = hourlp[nt]
 day0 = day0lp[nt]

# set these
 year = '2004'
 month = '09'

 ref_year = '%s' %year    #valid year
 ref_month = '%s' %month       #valid month
 ref_day = '%s' %day       #valid day
 ref_day0 = '%s' %day0       #valid calendar day
 ref_hour = '%s' %hour       #valid hour
 ref_min = 00
 ref_sec = 00

 print ref_day
 print ref_hour

 ppath = '/data/jtrapp/a/jtrapp/ivan/MIROCSLdataprep091404/'

#open netcdf file and read in vars and dims
 RH_file = ppath+'%s%s%s%s00.combine.%s%s%s.RH3D.nc' %(year,month,day,hour,year,month,day0)
 TMP_file = ppath+'%s%s%s%s00.combine.%s%s%s.TMP3D.nc' %(year,month,day,hour,year,month,day0)
 UGRD_file = ppath+'%s%s%s%s00.combine.%s%s%s.UGRD3D.nc' %(year,month,day,hour,year,month,day0)
 VGRD_file = ppath+'%s%s%s%s00.combine.%s%s%s.VGRD3D.nc' %(year,month,day,hour,year,month,day0)
 HGT_file = ppath+'%s%s%s%s00.combine.%s%s%s.HGT3D.nc' %(year,month,day,hour,year,month,day0)


 nc1 = Nio.open_file(RH_file,'r')
 nc2 = Nio.open_file(TMP_file,'r')
 nc3 = Nio.open_file(UGRD_file,'r')
 nc4 = Nio.open_file(VGRD_file,'r')
 nc5 = Nio.open_file(HGT_file,'r')

 RH = nc1.variables['RH'][:]
 lat = nc1.variables['lat_3'][:]
 lon = nc1.variables['lon_3'][:]
 plevs = nc1.variables['lv_ISBL7'][:]
 plevs2 = nc2.variables['lv_ISBL3'][:]
# Note:  the plevs here are in hPa, but the grib2 is considered to be in Pa
# so in the subsequent usages of p, multiply by 100
 print nc1.variables['lv_ISBL7']
 TMP = nc2.variables['TMP'][:]
 UGRD = nc3.variables['UGRD'][:]
 VGRD = nc4.variables['VGRD'][:]
 HGT = nc5.variables['HGT'][:]

 numlat = nc1.dimensions['lat_3']
 numlon = nc1.dimensions['lon_3']
 numPlevs = nc1.dimensions['lv_ISBL7']
# print nc1.variables['RH']

#varAttVal = RH.attributes['initial_time']
#varAttsAndVals = RH.attributes
#print ['varAttVal']
 nc1.close()
 nc2.close()
 nc3.close()
 nc4.close()
 nc5.close()

 RH = np.ma.array(RH,mask=False)
 HGT = np.ma.array(HGT,mask=False)
 TMP = np.ma.array(TMP,mask=False)
 UGRD = np.ma.array(UGRD,mask=False)
 VGRD = np.ma.array(VGRD,mask=False)


#open file to write to

 cpath = '/data/jtrapp/a/jtrapp/ivan/MIROCSLdata091404/'

 grbfile=cpath+'%s%s%s%s00.cdas1.%s%s%s.pgrbh.grb2' %(year,month,day,hour,year,month,day0)
#if file does not already exist, create it
 if os.path.isfile(grbfile) == False:
   f=open(grbfile,'wb')
#if file already exists, we are just going to append records to it
 else:
   f=open(grbfile,'a+')
#### add something that removes the file if it already exists...

##########################
# Identification Section #
##########################
# from decoded nam file
# grid id section [   7    0    2    1    1 2013    5   20    0    0    0    0    1]


 orig_id = 7           # 7=NWS NCEP
 sub_id = 4            # 4=EMC
 grb_master_table = 11 #use newest version.
 grb_local_table = 0   #local tables not used
 sig_ref_time = 2      #verifying time of forecast

 ref_min = 00
 ref_sec = 00
 prod_status = 2       #research products
 type_data = 1         # 1=forecast
 idsect = [orig_id,sub_id,grb_master_table,grb_local_table,sig_ref_time,ref_year,ref_month,ref_day,ref_hour,ref_min,ref_sec,prod_status,type_data]


###################
# Grid Definition #
###################

 grid_def = 0                 #Source of grid definition, 0: specified in Code Table (grid_def_template_number)
 num_gridpts = numlat*numlon  #Number of grid points
 num_octets = 0               #0: regular grid
 interp_optpts = 0            #there is no appended list
 grid_def_template_number = 0 #lat/lon grid
 gdsinfo = [grid_def, num_gridpts, num_octets, interp_optpts, grid_def_template_number]

 earthshape = 6               #Earth assumed spherical with radius=6371229 m
 sf_Earthradius = 0           #scale factor for radius of Earth

 gdtmpl = [earthshape, sf_Earthradius, 0, 0, 0, 0, 0, numlon, numlat, 0, 0 , lat[0]*1e6,lon[0]*1e6,48,lat[numlat-1]*1e6,lon[numlon-1]*1e6,\
         (lon[1]-lon[0])*1e6,(lat[numlat-1]-lat[numlat-2])*1e6,0]

 drtnum = 40                  #Data Representation Template Number- 40: Grid Point Data -JPEG2000 Compression

 for varID in ['R3D','T3D','U3D','V3D','H3D']:

    pdtnum = 0                   #Analysis or forecast at a horizontal level or in a horizontal layer at a point in time.
    type_gen = 2

    if varID == 'R3D':
      for i,p in enumerate(plevs):
        fixed_sfc = 100            #type of surface, 100: pressure levels (Pa)
        drtmpl = [0, 0, 0, 7, 0, 0, 255] #data packing/compression
        discipline =  0            #0:meteorological, 1:hydrological, 2:land surface, 3:space, 10:ocean
        parameter_cat = 1          # 0: temp, 1: moisture, 2: momentum, 3: mass
        parameter_num = 1          # 1: relative humidity (%)
#        pdtmpl = [parameter_cat, parameter_num, type_gen, 255, 255, 0, 0, 1, 0, fixed_sfc, 0, p, 255, 0, 0]
# need to convert p to Pa
        pdtmpl = [parameter_cat, parameter_num, type_gen, 255, 255, 0, 0, 1, 0, fixed_sfc, 0, p*100, 255, 0, 0]
        grib2 = Grib2Encode(discipline,idsect)
        grib2.addgrid(gdsinfo,gdtmpl)
        grib2.addfield(pdtnum,pdtmpl,drtnum,drtmpl,RH[i,:,:])
        grib2.end()
        f.write(grib2.msg)

    elif varID == 'T3D':
      for i,p in enumerate(plevs2):
        fixed_sfc = 100            #type of surface, 100: pressure levels (Pa)
        drtmpl = [1158811648, 0, 1, 10, 0, 0, 255]
        discipline =  0            #0:meteorological, 1:hydrological, 2:land surface, 3:space, 10:ocean
        parameter_cat = 0          # 0: temp, 1: moisture, 2: momentum, 3: mass
        parameter_num = 0          # 0: temp in K
        pdtmpl = [parameter_cat, parameter_num, type_gen, 255, 255, 0, 0, 1, 0, fixed_sfc, 0, p*100, 255, 0, 0]

        grib2 = Grib2Encode(discipline,idsect)
        grib2.addgrid(gdsinfo,gdtmpl)
        grib2.addfield(pdtnum,pdtmpl,drtnum,drtmpl,TMP[i,:,:])
        grib2.end()
        f.write(grib2.msg)

    elif varID == 'U3D':
      for i,p in enumerate(plevs2):
        fixed_sfc = 100            #type of surface, 100: pressure levels (Pa)
        drtmpl = [-986906624, 0, 2, 13, 0, 0, 255]
        discipline =  0            #0:meteorological, 1:hydrological, 2:land surface, 3:space, 10:ocean
        parameter_cat = 2          # 0: temp, 1: moisture, 2: momentum, 3: mass
        parameter_num = 2          # 2: u-component of wind speed m/s
        pdtmpl = [parameter_cat, parameter_num, type_gen, 255, 255, 0, 0, 1, 0, fixed_sfc, 0, p*100, 255, 0, 0]

        grib2 = Grib2Encode(discipline,idsect)
        grib2.addgrid(gdsinfo,gdtmpl)
        grib2.addfield(pdtnum,pdtmpl,drtnum,drtmpl,UGRD[i,:,:])
        grib2.end()
        f.write(grib2.msg)

    elif varID == 'V3D':
      for i,p in enumerate(plevs2):
        fixed_sfc = 100            #type of surface, 100: pressure levels (Pa)
        drtmpl = [-986906624, 0, 2, 13, 0, 0, 255]
        discipline =  0            #0:meteorological, 1:hydrological, 2:land surface, 3:space, 10:ocean
        parameter_cat = 2          # 0: temp, 1: moisture, 2: momentum, 3: mass
        parameter_num = 3          # 3: v-component of wind speed m/s
        pdtmpl = [parameter_cat, parameter_num, type_gen, 255, 255, 0, 0, 1, 0, fixed_sfc, 0, p*100, 255, 0, 0]
        grib2 = Grib2Encode(discipline,idsect)
        grib2.addgrid(gdsinfo,gdtmpl)
        grib2.addfield(pdtnum,pdtmpl,drtnum,drtmpl,VGRD[i,:,:])
        grib2.end()
        f.write(grib2.msg)

    elif varID == 'H3D':
      for i,p in enumerate(plevs2):
        fixed_sfc = 100            #type of surface, 100: pressure levels (Pa)
        drtmpl = [1235854176, 4, 3, 16, 0, 0, 255]
        discipline =  0            #0:meteorological, 1:hydrological, 2:land surface, 3:space, 10:ocean
        parameter_cat = 3          # 0: temp, 1: moisture, 2: momentum, 3: mass
        parameter_num = 5          # 3: v-component of wind speed m/s
        pdtmpl = [parameter_cat, parameter_num, type_gen, 255, 255, 0, 0, 1, 0, fixed_sfc, 0, p*100, 255, 0, 0]
        grib2 = Grib2Encode(discipline,idsect)
        grib2.addgrid(gdsinfo,gdtmpl)
        grib2.addfield(pdtnum,pdtmpl,drtnum,drtmpl,HGT[i,:,:])
        grib2.end()
        f.write(grib2.msg)

    else:
      print varID
      print "something went wrong"

 f.close()

 print "done writing 3D fields to grib file"
