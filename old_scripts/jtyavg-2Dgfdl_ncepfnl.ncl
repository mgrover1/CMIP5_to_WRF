; script to compute decadal average of 2D (sfc) atms fields from CMIP5
; the assumption here is that the script is being applied to
; gfdl output
; latest update JT 8/21/17
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

;Set Path
main1 = "/data/jtrapp/a/jtrapp/CMIP5/GFDL/"
main2 = "/data/jtrapp/a/crrllsm2/ivandata/ncepfnl/onedeg_ncepfnl/"


vars = (/"uas","vas","tas","hurs","huss","ts","ps","psl"/)

do nv = 0,7
;do nv = 0,0

; set variable name
vname=vars(nv)

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


; FNL grid is 1 deg lat by 1 deg lon
flat  = nmf->lat_3
flon  = nmf->lon_3

printVarSummary(flat)
printVarSummary(flon)
dimlatc   = dimsizes(flat)
dimlonc   = dimsizes(flon)
nlat    = dimlatc(0)
mlon    = dimlonc(0)


; set up all of the plotting parameters, for later use
wks = gsn_open_wks("pdf","2dplot-gfdl")
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

np = 0

; get the data for the first September
nf = 1.
ny = 8
fvaravg  = ff1->$vname$(ny,:,:)  ; time, lat, lon
hvaravg  = hf1->$vname$(ny,:,:)  ; time, lat, lon


; this should loop through all subsequent Septembers in the dataset
do ny = 20,nmax-1,12

   print("mon # "+ny)

; check for any missing data...
   if(any(.not.ismissing(ff1))) then

   fvar  = ff1->$vname$(ny,:,:)  ; time, lat, lon
   hvar  = hf1->$vname$(ny,:,:)  ; time, lat, lon

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

; this should loop through all subsequent Mays in the dataset
do ny = 8,nmax-1,12

   print("mon # "+ny)


; check for any missing data...
   if(any(.not.ismissing(ff2))) then

   fvar  = ff2->$vname$(ny,:,:)  ; time, lat, lon
   hvar  = hf2->$vname$(ny,:,:)  ; time, lat, lon

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

end do    ; all months

print("nfiles in average= "+nf)

; form averages
fvaravg=fvaravg/nf
hvaravg=hvaravg/nf

; get difference (future-historical)
dvaravg=fvaravg-hvaravg
copy_VarMeta(fvaravg,dvaravg)

; check plot...
plot = gsn_csm_contour_map(wks,fvaravg,res)
plot = gsn_csm_contour_map(wks,hvaravg,res)
plot = gsn_csm_contour_map(wks,dvaravg,res)

; regrid the different field to the nam grid, using the following
; hard-coded parameters
; Bilinear remapping is used

Opt                = True
Opt@DstGridLat     = flat
Opt@DstGridLon     = flon
Opt@RemoveWgtFile  = True
Opt@DstForceOverwrite = True
Opt@SrcForceOverwrite = True
; this is to deal with the missing data
Opt@SrcMask2D      = where(.not.ismissing(hvaravg),1,0)

; regrid the future & past, then take difference
fvar_regrid = ESMF_regrid(fvaravg,Opt)
hvar_regrid = ESMF_regrid(hvaravg,Opt)

; get difference (future-historical)
dvar_regrid=fvar_regrid-hvar_regrid

copy_VarMeta(fvar_regrid,dvar_regrid)
printVarSummary(dvar_regrid)

  

plot = gsn_csm_contour_map(wks,hvar_regrid,res)
plot = gsn_csm_contour_map(wks,fvar_regrid,res)
plot = gsn_csm_contour_map(wks,dvar_regrid,res)

;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
; finally, do a basic dump of the difference field to a
; netcdf file
;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

  ncfile = vname + "_GFDLncepfnldiff.nc"
  system ("rm " + ncfile)               ; remove any pre-existing file
  ncdf   = addfile(ncfile ,"c")   ; open output netCDF file

  ncdf->$vname$  = dvar_regrid             

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
