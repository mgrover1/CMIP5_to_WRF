; script to compute decadal average of 3D atms fields from gfdl
; and then put on the ncepfnl grid
; use this particular version to loop through the different variables 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

; set up all of the plotting parameters, for later use
wks = gsn_open_wks("ncgm","3dplot-gfdl")
gsn_define_colormap(wks,"wh-bl-gr-ye-re")  ; choose a colormap

;vars = (/"ua","va","ta","zg"/)
vars = (/"hur"/)

; note: use a separate run of this script to average/grid the RH field 
; (because it has a
; different pressure-level distribution than the other variables)


; loop through variables...
;do nv = 0,3
do nv = 0,0

vname=vars(nv)

;Set Path
main1 = "/data/jtrapp/a/jtrapp/CMIP5/GFDL/"
main2 = "/data/jtrapp/a/crrllsm2/ivandata/ncepfnl/onedeg_ncepfnl/"

; eventually will set full file name from shell script 

ffile1 = vname+"_Amon_GFDL-CM3_rcp85_r1i1p1_209101-209512.nc"
ffile2 = vname+"_Amon_GFDL-CM3_rcp85_r1i1p1_209601-210012.nc"

hfile1 = vname+"_Amon_GFDL-CM3_historical_r1i1p1_199001-199412.nc"
hfile2 = vname+"_Amon_GFDL-CM3_historical_r1i1p1_199501-199912.nc"


; sample name file
nfile = "fnl_20040914_00_00.grib1"

;Pointer
ff1 = addfile(main1+ffile1, "r")
ff2 = addfile(main1+ffile2, "r")
hf1 = addfile(main1+hfile1, "r")
hf2 = addfile(main1+hfile2, "r")

nmf = addfile(main2+nfile, "r")

lat2d  = nmf->lat_3  ;(ygrid, xgrid)
lon2d  = nmf->lon_3  ;(ygrid, xgrid)

dimlatc   = dimsizes(lat2d)
dimlonc   = dimsizes(lon2d)
  nlat    = dimlatc(0)
  mlon    = dimlonc(0)

; read one nanl variable to get metadata
if(vname.eq."hur") then
   nanlvar = nmf->R_H_3_ISBL
else
   nanlvar = nmf->HGT_3_ISBL
end if

; set up all of the plotting parameters, for later use

res                        = True               ; plot mods desired for origin
res@cnFillOn               = True               ; color fill
res@cnLinesOn              = False              ; no contour lines
res@gsnSpreadColors        = True               ; use total colormap
res@mpGridAndLimbOn        = True
res@mpGridLineDashPattern  = 2                  ; lat/lon lines as dashed
res@pmTickMarkDisplayMode  = "Always"           ; turn on tickmarks

; state borders...
res@mpPerimOn           = True
res@mpOutlineBoundarySets = "GeophysicalAndUSStates"

;  res@gsnAddCyclic          = False    ; data already has cyclic point
				       ; this must also be set for any zoom

; note that the gsn_csm_*map_ce templates automatically set 
; res@mpLimitMode="LatLon" for you. If you are plotting a different projection,
; you may have to set this resource.

;  res@mpMinLatF            = -60      ; range to zoom in on
;  res@mpMaxLatF            =  30.
;  res@mpMinLonF            =  30.
;  res@mpMaxLonF            = 120.


; create new arrays, using hardcoded MIROC5 metadata
; Note:  3D files are (time,lev,lat,lon), 2D are (time,lat,lon)
;nlat = 128
;nlon = 256
;nlev = 17  
;      plev = 100000, 92500, 85000, 70000, 60000, 50000, 40000, 30000, 25000, 
;      20000, 15000, 10000, 7000, 5000, 3000, 2000, 1000 
; the following is for vertical interpoloation (pressure level
; to pressure level)

;linlog = 2   ; ln(p) interpolation
; will need to check this...DCS says that linlong=2 didn't work in real.exe
linlog = -2   ; ln(p) interpolation; in this case, set negative to allow for "upward" extrapolation 

pdim = 0    ; location in array of p-coordinate
; these are the pressure levels in the GFDL data
pin = (/100000, 92500, 85000, 77500, 70000, 60000, 50000, 40000, 30000, \
    25000, 20000, 15000, 10000, 7000, 5000, 3000, 2000, 1000, 700, 500, 300, \
    200, 100 /)

; NCEPFNL
if(vname.eq."hur") then
printVarSummary(vname)
pout = (/100000,97500,95000,92500,90000,85000,\
       80000,75000,70000,65000,\
       60000,55000,50000,45000,40000,\
       35000,30000,25000,20000,\
       15000,10000/)
else
printVarSummary(vname)
pout = (/ 100000,97500,95000,92500,90000,85000,\
       80000,75000,70000,65000,\
       60000,55000,50000,45000,40000,\
       35000,30000,25000,20000,\
       15000,10000,7000,5000,3000,2000,1000/)
end if

pout!0 = "p"
pout@units = "Pa"

dimlp   = dimsizes(pout)
nlevs    = dimlp(0)

print("pout levs="+nlevs)
dimlp   = dimsizes(pin)
ninlevs    = dimlp(0)
print("pin levs="+ninlevs)


nyr = 5   ; 3D files are in 1/2 decadal segments
nmax = nyr*12
print("nmax= # "+nmax)

np = 0

; get the data for the first September
nf = 1.
ny = 8
fvaravg  = ff1->$vname$(ny,:,:,:)  ; time, plev, lat, lon
hvaravg  = hf1->$vname$(ny,:,:,:)  ; time, plev, lat, lon


; this should loop through all subsequent Septembers in the dataset
do ny = 20,nmax-1,12

   print("mon # "+ny)

; check for any missing data...
   if(any(.not.ismissing(ff1))) then

   fvar  = ff1->$vname$(ny,:,:,:)  ; time, plev, lat, lon
   hvar  = hf1->$vname$(ny,:,:,:)  ; time, plev, lat, lon

   lat   = ff1->lat
   lon   = ff1->lon
   plev   = ff1->plev

   fvaravg = fvaravg + fvar
   hvaravg = hvaravg + hvar

   nf = nf +1.

  else
    print("Some of the data are missing.  Skipping this file.")
  end if

end do    ; all months

; loop through all Septembers in the second file
do ny = 8,nmax-1,12

   print("mon # "+ny)

; check for any missing data...
   if(any(.not.ismissing(ff2))) then

   fvar  = ff2->$vname$(ny,:,:,:)  ; time, plev, lat, lon
   hvar  = hf2->$vname$(ny,:,:,:)  ; time, plev, lat, lon

   lat   = ff2->lat
   lon   = ff2->lon
   plev   = ff2->plev

   fvaravg = fvaravg + fvar
   hvaravg = hvaravg + hvar

   nf = nf +1.

  else
    print("Some of the data are missing.  Skipping this file.")
  end if

end do    ; all months

print("nfiles in average= "+nf)

; form averages
fvaravg=fvaravg/nf
hvaravg=hvaravg/nf

; get difference (future-historical)
dvaravg=fvaravg-hvaravg

; regrid the different field to the FNL grid, using the following
; hard-coded parameters
;
; Bilinear remapping is used

fvar_regrid = new((/ninlevs,nlat,mlon/),"float",1e+20)
hvar_regrid = new((/ninlevs,nlat,mlon/),"float",1e+20)
printVarSummary(fvar_regrid)
printVarSummary(fvaravg)

do nnn = 0, ninlevs-1

print("interp lev="+nnn)

Opt                = True
Opt@DstGridLat     = lat2d
Opt@DstGridLon     = lon2d
Opt@RemoveWgtFile  = True
Opt@DstForceOverwrite = True
Opt@SrcForceOverwrite = True
; this is to deal with the missing data
Opt@SrcMask2D      = where(.not.ismissing(hvaravg(nnn,:,:)),1,0)

; regrid the future & past, then take difference
fvar_regrid(nnn,:,:) = ESMF_regrid(fvaravg(nnn,:,:),Opt)
hvar_regrid(nnn,:,:) = ESMF_regrid(hvaravg(nnn,:,:),Opt)

end do

;dvar_vint = where(ismissing(dvar_vint),0,dvar_vint)

; get difference (future-historical)
dvar_regrid=fvar_regrid-hvar_regrid


copy_VarMeta(fvar_regrid,dvar_regrid)
printVarSummary(dvar_regrid)

  res@gsnAddCyclic          = False    ; data already has cyclic point
                                       ; this must also be set for any zoom

;  res@mpLimitMode       = "Corners"         ; choose range of map
;  res@mpLeftCornerLatF  = lat2d(0,0)
;  res@mpLeftCornerLonF  = lon2d(0,0)
;  res@mpRightCornerLatF = lat2d(nlat-1,mlon-1)
;  res@mpRightCornerLonF = lon2d(nlat-1,mlon-1)
;  res@tfDoNDCOverlay = True
;  res@mpProjection           = "LambertConformal"
;  res@mpLambertParallel1F    = lat2d@mpLambertParallel1F
;  res@mpLambertParallel2F    = lat2d@mpLambertParallel2F
;  res@mpLambertMeridianF     = lat2d@mpLambertMeridianF
;do nnn = 0, ninlevs-1
;plot = gsn_csm_contour_map(wks,fvar_regrid(2,:,:),res)
;plot = gsn_csm_contour_map(wks,hvar_regrid(2,:,:),res)
;plot = gsn_csm_contour_map(wks,dvar_regrid(nnn,:,:),res)
;plot = gsn_csm_contour_map(wks,dvar_regrid_m(nnn,:,:),res)
;end do

; now vertically interpolate
; to be consistent with the cmip5 data, interpolate assuming
; that pressure dereases upward, and then reverse, to be 
; consistent with cmip5 data

temp = int2p_n_Wrap(pin,dvar_regrid,pout,linlog,pdim)

; reverse direction of vertical coordinate so that it's top to bottom
dvar_vint = temp(::-1,:,:)
;copy_VarMeta(nanlvar,dvar_vint)
printVarSummary(dvar_vint)

; one final step:  set missing values to zero (i.e., a missing is simply a zero difference)
; this will ensure that the differences can be added everywhere to the original CDAS
; removed - afraid this will cause odd gradients...

;dvar_vint = where(ismissing(dvar_vint),0,dvar_vint)

; go through and plot all levels
do n = 0,nlevs-1

plot = gsn_csm_contour_map(wks,dvar_vint(n,:,:),res)

end do

;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
; finally, do a basic dump of the difference field to a
; netcdf file
;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

  ncfile = vname + "_GFDLncepfnldiff.nc"
  system ("rm " + ncfile)               ; remove any pre-existing file
  ncdf   = addfile(ncfile ,"c")   ; open output netCDF file

  ncdf->$vname$  = dvar_vint              
  delete(ncfile)
  delete(ncdf)

  ncfile = vname + "_GFDLncepfnl.nc"
  system ("rm " + ncfile)               ; remove any pre-existing file
  ncdf   = addfile(ncfile ,"c")   ; open output netCDF file

  ncdf->fvaravg  = fvar_regrid
  ncdf->hvaravg  = hvar_regrid

  delete(ncfile)
  delete(ncdf)

  delete(res)


end do


end
