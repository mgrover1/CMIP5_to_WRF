#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')
import numpy as np
from mpl_toolkits.basemap import Basemap
from ncepgrib2 import Grib2Encode, Grib2Decode
import pygrib
import os
import Nio

# Note:  this version is designed specifically for FNL (and other lat/lon grids)


hourlp = ['00','06','12','18','00','06','12','18','00','06','12','18','00','06','12','18','00','06','12','18','00']
day0lp = ['14','14','14','14','15','15','15','15','16','16','16','16','17','17','17','17','18','18','18','18','19']
daylp = ['14','14','14','14','15','15','15','15','16','16','16','16','17','17','17','17','18','18','18','18','19']


for nt in range(len(daylp)):

# loop through...
 day = daylp[nt]
 hour = hourlp[nt]
 day0 = day0lp[nt]

# set these
 year = '2004'
 month = '09'

 ref_year = '%s' %year    #valid year
 ref_month = '%s' %month       #valid month
 ref_day = '%s' %day       #valid day
 ref_hour = '%s' %hour       #valid hour
 ref_day0 = '%s' %day       #valid day
 ref_min = 00
 ref_sec = 00

 print ref_day
 print ref_hour

 ppath = '/data/jtrapp/a/jtrapp/ivan/MIROCSLdataprep091404/'


#open netcdf file and read in vars and dims
 GRND_file = ppath+'%s%s%s%s00.combine.%s%s%s.GRNDsfc.nc' %(year,month,day,hour,year,month,day0)
 ICE_file = ppath+'%s%s%s%s00.combine.%s%s%s.ICEsfc.nc' %(year,month,day,hour,year,month,day0)
 LAND_file =ppath+'%s%s%s%s00.combine.%s%s%s.LANDsfc.nc' %(year,month,day,hour,year,month,day0)
 PRES_file = ppath+'%s%s%s%s00.combine.%s%s%s.PRESsfc.nc' %(year,month,day,hour,year,month,day0)
 SOILW_file = ppath+'%s%s%s%s00.combine.%s%s%s.SOILWsfc.nc' %(year,month,day,hour,year,month,day0)
 TSOIL_file = ppath+'%s%s%s%s00.combine.%s%s%s.TSOILsfc.nc' %(year,month,day,hour,year,month,day0)
 SKIN_file = ppath+'%s%s%s%s00.combine.%s%s%s.SKINsfc.nc' %(year,month,day,hour,year,month,day0)
 PSL_file = ppath+'%s%s%s%s00.combine.%s%s%s.PSLsfc.nc' %(year,month,day,hour,year,month,day0)
 SPFH_file = ppath+'%s%s%s%s00.combine.%s%s%s.SPFHsfc.nc' %(year,month,day,hour,year,month,day0)
 TMP_file = ppath+'%s%s%s%s00.combine.%s%s%s.TMPsfc.nc' %(year,month,day,hour,year,month,day0)
 UGRD_file = ppath+'%s%s%s%s00.combine.%s%s%s.UGRDsfc.nc' %(year,month,day,hour,year,month,day0)
 VGRD_file = ppath+'%s%s%s%s00.combine.%s%s%s.VGRDsfc.nc' %(year,month,day,hour,year,month,day0)
 RH_file = ppath+'%s%s%s%s00.combine.%s%s%s.RHsfc.nc' %(year,month,day,hour,year,month,day0)

 nc1 = Nio.open_file(SPFH_file,'r')
 nc2 = Nio.open_file(TMP_file,'r')
 nc3 = Nio.open_file(UGRD_file,'r')
 nc4 = Nio.open_file(VGRD_file,'r')
 nc5 = Nio.open_file(PRES_file,'r')
 nc6 = Nio.open_file(PSL_file,'r')
 nc7 = Nio.open_file(SKIN_file,'r')
 nc8 = Nio.open_file(TSOIL_file,'r')
 nc9 = Nio.open_file(SOILW_file,'r')
 nc10 = Nio.open_file(LAND_file,'r')
 nc11 = Nio.open_file(GRND_file,'r')
 nc12 = Nio.open_file(RH_file,'r')
 SPFH = nc1.variables['SPFH'][:]
 lat = nc3.variables['lat_3'][:]
 lon = nc3.variables['lon_3'][:]
 TMP = nc2.variables['TMP'][:]
 UGRD = nc3.variables['UGRD'][:]
 VGRD = nc4.variables['VGRD'][:]
 PRES = nc5.variables['PRES'][:]
 PSL = nc6.variables['PSL'][:]
 SKIN = nc7.variables['SKIN'][:]
 TSOIL = nc8.variables['TMP'][:]
 SOILW = nc9.variables['SOILW'][:]
#depth1 = nc8.variables['lv_DBLL0'][:]
#depth0 = nc8.variables['lv_DBLL0_l0'][:]
 LAND = nc10.variables['LAND'][:]
 GRND = nc11.variables['GRND'][:]
 RH = nc12.variables['RH'][:]
 numlat = nc3.dimensions['lat_3']
 numlon = nc3.dimensions['lon_3']
# just hard code these...
# numlat = 181
# numlon = 360

 nc1.close()
 nc2.close()
 nc3.close()
 nc4.close()
 nc5.close()
 nc6.close()
 nc7.close()
 nc8.close()
 nc9.close()
 nc10.close()
 nc11.close()
 nc12.close()

# this is needed (only for the surface data), otherwise incorrect bitmapping
 SPFH = np.ma.array(SPFH,mask=False)
 TMP = np.ma.array(TMP,mask=False)
 UGRD = np.ma.array(UGRD,mask=False)
 VGRD = np.ma.array(VGRD,mask=False)
 PSL = np.ma.array(PSL,mask=False)
 PRES = np.ma.array(PRES,mask=False)
 SKIN = np.ma.array(SKIN,mask=False)
 LAND = np.ma.array(LAND,mask=False)
 GRND = np.ma.array(GRND,mask=False)
 SOILW = np.ma.array(SOILW,mask=False)
 TSOIL = np.ma.array(TSOIL,mask=False)
 RH = np.ma.array(RH,mask=False)

# levels for the soil moisture and temperature
#levs1 = [.0,.10,.40,1.0]
#levs2 = [.10,.40,1.0,2.0]
 levs1 = [0,10]
 levs2 = [10,200]


#open file to write to
 cpath = '/data/jtrapp/a/jtrapp/ivan/MIROCSLdata091404/'

 grbfile=cpath+'%s%s%s%s00.cdas1.%s%s%s.sfluxgrbf.grb2' %(year,month,day,hour,year,month,day0)

 if os.path.isfile(grbfile) == False:
   f=open(grbfile,'wb')
 else:
   f=open(grbfile,'a+')

##########################
# Identification Section #
##########################
 orig_id = 7           # 7=NWS NCEP
 sub_id = 4            # 4=EMC
 grb_master_table = 11 #use newest version.
 grb_local_table = 0   #local tables not used
 sig_ref_time = 2      #verifying time of forecast

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
 grid_def_template_number = 0 #lat/lon grid (Equidistant Cylindrical
 gdsinfo = [grid_def, num_gridpts, num_octets, interp_optpts, grid_def_template_number]

#gdtmpl contains data values for the specified Grid Def. Template
#http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_temp3-0.shtml
 earthshape = 6               #Earth assumed spherical with radius=6371229 m
 sf_Earthradius = 0           #scale factor for radius of Earth
# from KH's RH script
#gdtmpl = [earthshape, sf_Earthradius, 0, 0, 0, 0, 0, numlon, numlat, 0, 0 , lat[numlat-1]*1e6, lon[0]*1e6,48,lat[0]*1e6, \
#         lon[numlon-1]*1e6,(lon[1]-lon[0])*1e6,(lat[numlat-1]-lat[numlat-2])*1e6,scanmodeflag]
# from my CDAS script...
 gdtmpl = [earthshape, sf_Earthradius, 0, 0, 0, 0, 0, numlon, numlat, 0, 0 , lat[0]*1e6,lon[0]*1e6,48,lat[numlat-1]*1e6,lon[numlon-1]*1e6,\
         (lon[1]-lon[0])*1e6,(lat[numlat-1]-lat[numlat-2])*1e6,0]
# gdtmpl = [earthshape, sf_Earthradius, 0, 0, 0, 0, 0, numlon, numlat, 0, 0 , lat[numlat-1]*1e6, lon[0]*1e6,48,lat[0]*1e6, \
#         lon[numlon-1]*1e6,(lon[1]-lon[0])*1e6,(lat[numlat-1]-lat[numlat-2])*1e6,64]
 drtnum = 40

 for varID in ['T2','Q2','U10','V10','PSL','PRES','SKIN','LAND','GRND','ST','SM','RH']:
           pdtnum = 0 #Analysis or forecast at a horizontal level or in a horizontal layer at a point in time.
           type_gen = 2 #
          ######################
          # Product Definition #
          ######################

           if varID == 'T2':
             fixed_sfc = 103 # 103: specified height level above ground m
             scalef = 0 #scale factor
             height = 2 #2m
             drtmpl = [1185291264, 0, 2, 14, 0, 0, 255]
             parameter_cat = 0 # 0: temp, 1: moisture, 2: momentum, 3: mass
             parameter_num = 0 # temp in K
             pdtmpl = [parameter_cat, parameter_num, type_gen, 255, 255, 0, 0, 1, 0, fixed_sfc,scalef,height,255,0,0]
             grib2 = Grib2Encode(0,idsect)
             grib2.addgrid(gdsinfo,gdtmpl)
             grib2.addfield(pdtnum,pdtmpl,drtnum,drtmpl,TMP)
             grib2.end()
             f.write(grib2.msg)

           elif varID == 'Q2': #newQonPlevs:
             fixed_sfc = 103 # 103: specified height level above ground m
             scalef = 0 #scale factor
             height = 2 #2m
             drtmpl = [1065353216, 0, 5, 12, 0, 0, 255]
             parameter_cat = 1 # 0: temp, 1: moisture, 2: momentum, 3: mass
             parameter_num = 0 # specific humidity in kg/kg
             pdtmpl = [parameter_cat, parameter_num, type_gen, 255, 255, 0, 0, 1, 0, fixed_sfc,scalef,height,255,0,0]
             grib2 = Grib2Encode(0,idsect)
             grib2.addgrid(gdsinfo,gdtmpl)
             grib2.addfield(pdtnum,pdtmpl,drtnum,drtmpl,SPFH)
             grib2.end()
             f.write(grib2.msg)

           elif varID == 'U10':
             fixed_sfc = 103 # 103: specified height level above ground m
             scalef = 0 #scale factor
             height = 10 #2m
             drtmpl = [-986906624, 0, 2, 13, 0, 0, 255]
             parameter_cat = 2 # 0: temp, 1: moisture, 2: momentum, 3: mass
             parameter_num = 2 # 2: m/s
             pdtmpl = [parameter_cat, parameter_num, type_gen, 255, 255, 0, 0, 1, 0, fixed_sfc,scalef,height,255,0,0]
             grib2 = Grib2Encode(0,idsect)
             grib2.addgrid(gdsinfo,gdtmpl)
             grib2.addfield(pdtnum,pdtmpl,drtnum,drtmpl,UGRD)
             grib2.end()
             f.write(grib2.msg)

           elif varID == 'V10':
             fixed_sfc = 103 # 103: specified height level above ground m
             scalef = 0 #scale factor
             height = 10 #2m
             drtmpl = [-986148864, 0, 2, 13, 0, 0, 255]
             parameter_cat = 2 # 0: temp, 1: moisture, 2: momentum, 3: mass
             parameter_num = 3 # 2: m/s
             pdtmpl = [parameter_cat, parameter_num, type_gen, 255, 255, 0, 0, 1, 0, fixed_sfc,scalef,height,255,0,0]
             grib2 = Grib2Encode(0,idsect)
             grib2.addgrid(gdsinfo,gdtmpl)
             grib2.addfield(pdtnum,pdtmpl,drtnum,drtmpl,VGRD)
             grib2.end()
             f.write(grib2.msg)

           elif varID == 'PSL':
             fixed_sfc = 101 # ground surface
             scalef = 0 #scale factor
             height = 0 #ground surface
             drtmpl = [1231667856, 1, 1, 16, 0, 0, 255]
             parameter_cat = 3 # 0: temp, 1: moisture, 2: momentum, 3: mass
             parameter_num = 1 #192 pressure reduced to MSL
             pdtmpl = [parameter_cat, parameter_num, type_gen, 255, 255, 0, 0, 1, 0, fixed_sfc,scalef,height,255,0,0]
             grib2 = Grib2Encode(0,idsect)
             grib2.addgrid(gdsinfo,gdtmpl)
             grib2.addfield(pdtnum,pdtmpl,drtnum,drtmpl,PSL)
             grib2.end()
             f.write(grib2.msg)

           elif varID == 'PRES':
             fixed_sfc = 1
             scalef = 0 #scale factor
             height = 0 #ground surface
             drtmpl = [1224376704, 4, 1, 16, 0, 0, 255]
             parameter_cat = 3 # 0: temp, 1: moisture, 2: momentum, 3: mass
             parameter_num = 0 #
             pdtmpl = [parameter_cat, parameter_num, type_gen, 0, 81, 0, 0, 1, 0, fixed_sfc,scalef,height,255,0,0]
 #[  3   0   2   0  81   0   0   1   0   1   0   0 255   0   0]
             grib2 = Grib2Encode(0,idsect)
             grib2.addgrid(gdsinfo,gdtmpl)
             grib2.addfield(pdtnum,pdtmpl,drtnum,drtmpl,PRES)
             grib2.end()
             f.write(grib2.msg)

           elif varID == 'SKIN':
             fixed_sfc = 1 # ground or water surface
             scalef = 0 #scale factor
             height = 0 #ground surface
             drtmpl = [1157423104, 0, 1, 11, 0, 0, 255]
             parameter_cat = 0 # 0: temp, 1: moisture, 2: momentum, 3: mass
             parameter_num = 0 # temp in K
             pdtmpl = [parameter_cat, parameter_num, type_gen, 255, 255, 0, 0, 1, 0, fixed_sfc,scalef,height,255,0,0]
             grib2 = Grib2Encode(0,idsect)
             grib2.addgrid(gdsinfo,gdtmpl)
             grib2.addfield(pdtnum,pdtmpl,drtnum,drtmpl,SKIN)
             grib2.end()
             f.write(grib2.msg)

           elif varID == 'LAND':
             fixed_sfc = 1 # ground or water surface
             scalef = 0 #scale factor
             height = 0 #ground surface
             drtmpl = [0, 0, 0, 1, 0, 0, 255]
             parameter_cat = 0 # 0: temp, 1: moisture, 2: momentum, 3: mass
             parameter_num = 0 # temp in K
             pdtmpl = [parameter_cat, parameter_num, type_gen, 255, 255, 0, 0, 1, 0, fixed_sfc,scalef,height,255,0,0]
             discipline = 2 #0:meteorological, 1:hydrological, 2:land surface, 3:space, 10:ocean
             grib2 = Grib2Encode(discipline,idsect)
             grib2.addgrid(gdsinfo,gdtmpl)
             grib2.addfield(pdtnum,pdtmpl,drtnum,drtmpl,LAND)
             grib2.end()
             f.write(grib2.msg)

           elif varID == 'GRND':
             fixed_sfc = 1 # ground or water surface
             scalef = 0 #scale factor
             height = 0 #ground surface
             drtmpl = [-969024513, 4, 2, 16, 0, 0, 255]
             parameter_cat = 3 # 0: temp, 1: moisture, 2: momentum, 3: mass
             parameter_num = 5 # height (m)
             pdtmpl = [parameter_cat, parameter_num, type_gen, 255, 255, 0, 0, 1, 0, fixed_sfc,scalef,height,255,0,0]
             discipline = 0 #0:meteorological, 1:hydrological, 2:land surface, 3:space, 10:ocean
             grib2 = Grib2Encode(discipline,idsect)
             grib2.addgrid(gdsinfo,gdtmpl)
             grib2.addfield(pdtnum,pdtmpl,drtnum,drtmpl,GRND)
             grib2.end()
             f.write(grib2.msg)

           elif varID == 'ST':
             for i in range(len(levs1)):
#            for i in range(0, 4):
#             for i,d in enumerate(levs1):
#            for l in range(len(lev)):
#              if l == 0:
#                height1 = .0
#                height2 = .10
#
#              elif l == 1:
#                height1 = .10
#                height2 = .40
#
#              elif l == 2:
#                height1 = .40
#                height2 = 1.00
#
#              else:
#                height1 = 1.00
#                height2 = 2.00
               print i, levs1[i], levs2[i]

               fixed_sfc = 106 # 106: depth below land surface
#              scale_factor = 0 #scale factor
               scale_factor = 2 #scale factor
               drtmpl = [1158152192, 0, 1, 10, 0, 0, 255]
               discipline = 0 #2 #0:meteorological, 1:hydrological, 2:land surface, 3:space, 10:ocean
               parameter_cat = 0 # 3 #0 #
               parameter_num = 0 #2 # temp in K
               pdtmpl = (parameter_cat, parameter_num, type_gen, 255, 255, 0, 0, 1, 0, fixed_sfc,scale_factor,levs1[i],fixed_sfc,scale_factor,levs2[i])
#              pdtmpl = [parameter_cat, parameter_num, type_gen, 255, 255, 0, 0, 1, 0, fixed_sfc,scale_factor,height1,fixed_sfc,scale_factor,height2]
               grib2 = Grib2Encode(discipline,idsect)
               grib2.addgrid(gdsinfo,gdtmpl)
#              grib2.addfield(pdtnum,pdtmpl,drtnum,drtmpl,TSOIL[l,:,:])
               grib2.addfield(pdtnum,pdtmpl,drtnum,drtmpl,TSOIL[i,:,:])
               grib2.end()
               f.write(grib2.msg)

           elif varID == 'SM':
#            for l in range(len(lev)):
#            for i in range(0, 1):
             for i in range(len(levs1)):
#              if l == 0:
#                height1 = .0
#                height2 = .10
#
#              elif l == 1:
#                height1 = .10
#                height2 = .40
#
#              elif l == 2:
#                height1 = .40
#                height2 = 1.00
#
#              else:
#                height1 = 1.00
#                height2 = 2.00

               fixed_sfc = 106 # 106: depth below land surface
#              scale_factor = 0 #scale factor
               scale_factor = 2 #scale factor
               drtmpl = [1107820544, 0, 3, 10, 0, 0, 255]
               discipline = 2 #2 #0:meteorological, 1:hydrological, 2:land surface, 3:space, 10:ocean
               parameter_cat = 0 # soil products
               parameter_num = 192 #volumetric soil moisture content fraction
#              pdtmpl = [parameter_cat, parameter_num, type_gen, 255, 255, 0, 0, 1, 0, fixed_sfc,scale_factor,height1,fixed_sfc,scale_factor,height2]
               pdtmpl = (parameter_cat, parameter_num, type_gen, 255, 255, 0, 0, 1, 0, fixed_sfc,scale_factor,levs1[i],fixed_sfc,scale_factor,levs2[i])
               grib2 = Grib2Encode(discipline,idsect)
               grib2.addgrid(gdsinfo,gdtmpl)
#               grib2.addfield(pdtnum,pdtmpl,drtnum,drtmpl,SOILW[l,:,:])
               grib2.addfield(pdtnum,pdtmpl,drtnum,drtmpl,SOILW[i,:,:])
               grib2.end()
               f.write(grib2.msg)

           elif varID == 'RH': #newQonPlevs:
             fixed_sfc = 103 # 103: specified height level above ground m
             scalef = 0 #scale factor
             height = 2 #2m
             drtmpl = [0, 0, 0, 7, 0, 0, 255]
##            drtmpl = [1065353216, 0, 5, 12, 0, 0, 255]
             parameter_cat = 1 # 0: temp, 1: moisture, 2: momentum, 3: mass
             parameter_num = 1 # relative humidity in %
             pdtmpl = [parameter_cat, parameter_num, type_gen, 255, 255, 0, 0, 1, 0, fixed_sfc,scalef,height,255,0,0]
             grib2 = Grib2Encode(0,idsect)
             grib2.addgrid(gdsinfo,gdtmpl)
             grib2.addfield(pdtnum,pdtmpl,drtnum,drtmpl,RH)
             grib2.end()
             f.write(grib2.msg)



           else:
             print varID
             print "something went wrong"

 f.close()

 print "done writing sfc fields to grib file"
