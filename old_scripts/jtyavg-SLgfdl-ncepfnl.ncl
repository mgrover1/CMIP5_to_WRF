; script to compute decadal average of 2D (sfc) atms fields from CMIP5
; the assumption here is that the script is being applied to
; GFDL output
; updated:  10/27/15
;
; this version is for soil temp and moisture

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

;Set Path
main1 = "/data/jtrapp/a/jtrapp/CMIP5/GFDL/"
main2 = "/data/jtrapp/a/crrllsm2/ivandata/ncepfnl/onedeg_ncepfnl/"

; FNL soil levels (these are layers)
NANLlevs = (/0.1, 2./)

vars = (/"tsl","mrlsl"/)

do nv = 0,1
;do nv = 1,1

; set variable name
vname=vars(nv)

; eventually will set full file name from shell script

ffile1 = vname+"_Lmon_GFDL-CM3_rcp85_r1i1p1_209101-209512.nc"
ffile2 = vname+"_Lmon_GFDL-CM3_rcp85_r1i1p1_209601-210012.nc"

hfile1 = vname+"_Lmon_GFDL-CM3_historical_r1i1p1_199001-199412.nc"
hfile2 = vname+"_Lmon_GFDL-CM3_historical_r1i1p1_199501-199912.nc"


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
landmask = nmf->LAND_3_SFC ;(land = 1, sea = 0)
dimlatc   = dimsizes(lat2d)
dimlonc   = dimsizes(lon2d)
nlat    = dimlatc(0)
mlon    = dimlonc(0)

; set up all of the plotting parameters, for later use
wks = gsn_open_wks("ncgm","soilplot-gfdl")
gsn_define_colormap(wks,"wh-bl-gr-ye-re")  ; choose a colormap

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

nyr = 5   ; assume a 1/2 decadal analysis
nmax = nyr*12

; get the data for the first May
nf = 1.
ny = 8

fvaravg  = ff1->$vname$(ny,:,:,:)  ; time, lat, lon
hvaravg  = hf1->$vname$(ny,:,:,:)  ; time, lat, lon

slevs   = (/dble2flt(hf1->depth)/)   ; soil levels, increasing downward
printVarSummary(slevs)
print("gfdl soil levs="+slevs)
nlevs    = dimsizes(slevs)
print("gfdl nlevs="+nlevs)

; this should loop through all subsequent Mays in the dataset
do ny = 20,nmax-1,12

   print("mon # "+ny)

; check for any missing data...
   if(any(.not.ismissing(ff1))) then

   fvar  = ff1->$vname$(ny,:,:,:)  ; time, lat, lon
   hvar  = hf1->$vname$(ny,:,:,:)  ; time, lat, lon

   lat   = ff1->lat
   lon   = ff1->lon

   timef = ff1->time(ny)
   timeh = hf1->time(ny)

   print("timef and timeh= "+timef+" "+timeh)

   fvaravg = fvaravg + fvar
   hvaravg = hvaravg + hvar

   nf = nf +1.

  else
    print("Some of the data are missing.  Skipping this file.")
  end if

end do    ; all months

; this should loop through all subsequent Septembers in the dataset
do ny = 8,nmax-1,12

   print("mon # "+ny)

; check for any missing data...
   if(any(.not.ismissing(ff2))) then

   fvar  = ff2->$vname$(ny,:,:,:)  ; time, lat, lon
   hvar  = hf2->$vname$(ny,:,:,:)  ; time, lat, lon

   lat   = ff2->lat
   lon   = ff2->lon

   timef = ff2->time(ny)
   timeh = hf2->time(ny)

   print("timef and timeh= "+timef+" "+timeh)

   fvaravg = fvaravg + fvar
   hvaravg = hvaravg + hvar

   nf = nf +1.
  else
    print("Some of the data are missing.  Skipping this file.")
  end if

end do


print("nfiles in average= "+nf)

; form averages
fvaravg=fvaravg/nf
hvaravg=hvaravg/nf

; check plot...
;plot = gsn_csm_contour_map_ce(wks,fvaravg(3,:,:),res)
;plot = gsn_csm_contour_map_ce(wks,hvaravg(3,:,:),res)

; get % increase instead of straight difference 
; note that the check for zero here is for soilmoisture, which
; is set to zero over water; the assumption is that if hvar is zero
; so also will be fvar, and therefore the quantity will be 1.
dvaravg = 1.+ (fvaravg-hvaravg)/where(hvaravg.ne.0., hvaravg , 1.)
copy_VarMeta(fvaravg,dvaravg)
;plot = gsn_csm_contour_map_ce(wks,dvaravg(3,:,:),res)
; lev 9 = 0.9 m

; regrid the different field to the nam grid, using the following
; hard-coded parameters
; Bilinear remapping is used

fvar_regrid = new((/nlevs,nlat,mlon/),"float",1e+20)
hvar_regrid = new((/nlevs,nlat,mlon/),"float",1e+20)

do nnn = 0, nlevs-1

;print("interp lev="+nnn)

Opt                = True
Opt@DstGridLat     = lat2d
Opt@DstGridLon     = lon2d
Opt@RemoveWgtFile  = True
Opt@DstForceOverwrite = True
Opt@SrcForceOverwrite = True
; this is to deal with the missing data
Opt@SrcMask2D      = where(.not.ismissing(hvaravg(nnn,:,:)),1,0)

; regrid the future & past, then take difference
;fvar_regrid = ESMF_regrid(fvaravg,Opt)
;hvar_regrid = ESMF_regrid(hvaravg,Opt)
fvar_regrid(nnn,:,:) = ESMF_regrid(fvaravg(nnn,:,:),Opt)
hvar_regrid(nnn,:,:) = ESMF_regrid(hvaravg(nnn,:,:),Opt)

end do

; get % difference (future-historical)
; same comment regarding water points
dvar_regrid = 1.+ (fvar_regrid-hvar_regrid)/where(hvar_regrid.ne.0., hvar_regrid , 1.)
;dvar_regrid=1.+(fvar_regrid-hvar_regrid)/hvar_regrid

; redundant setting, but put 'dvar' change to 1 over ocean points
; i.e., setting it to no change

do nnn = 0, nlevs-1

dvar_regrid(nnn,:,:)   = where(landmask.eq.1., dvar_regrid(nnn,:,:), 1.)

end do

copy_VarMeta(fvar_regrid,dvar_regrid)
printVarSummary(dvar_regrid)

;  res@gsnAddCyclic          = False    ; data already has cyclic point
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

;plot = gsn_csm_contour_map(wks,hvar_regrid(3,:,:),res)
;plot = gsn_csm_contour_map(wks,fvar_regrid(3,:,:),res)
;plot = gsn_csm_contour_map(wks,dvar_regrid(3,:,:),res)

; now, vertically average the gfdl levels to the FNL layers
; Note:   depth = 0.01, 0.04, 0.08, 0.125, 0.175, 0.25, 0.35, 0.5, 0.7,
;                 0.9, 1.2, 1.6,
;                 2, 2.4, 2.8, 3.5, 4.5, 5.5, 6.75, 8.75 ;

; just 2 in FNL
dvar_vgrid = new((/2,nlat,mlon/),"float",1e+20)

; for these, just average about the FNL 'level'
; .1 m
dvar_vgrid(0,:,:)=(dvar_regrid(2,:,:) + dvar_regrid(3,:,:))/2. 

; 2 m
dvar_vgrid(1,:,:)=(dvar_regrid(11,:,:) + dvar_regrid(12,:,:) + \
                   dvar_regrid(13,:,:))/3.


do k = 0,1
   copy_VarMeta(dvar_regrid(k,:,:),dvar_vgrid(k,:,:))
end do


dvar_vgrid!0 = "depth"
dvar_vgrid&depth = NANLlevs
printVarSummary(dvar_vgrid)

;plot = gsn_csm_contour_map_ce(wks,dvar_regrid(3,:,:),res)
;plot = gsn_csm_contour_map_ce(wks,dvar_vgrid(0,:,:),res)


;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
; finally, do a basic dump of the difference field to a
; netcdf file
;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

  ncfile = vname + "_GFDLncepfnldiff.nc"
  system ("rm " + ncfile)               ; remove any pre-existing file
  ncdf   = addfile(ncfile ,"c")   ; open output netCDF file

  ncdf->$vname$  = dvar_vgrid             

  delete(ncfile)
  delete(ncdf)

  ncfile = vname + "_GFDLncepfnl.nc"
  system ("rm " + ncfile)               ; remove any pre-existing file
  ncdf   = addfile(ncfile ,"c")   ; open output netCDF file

  ncdf->fvar_regrid  = fvar_regrid
  ncdf->hvar_regrid  = hvar_regrid

  delete(ncfile)
  delete(ncdf)

end do


end
