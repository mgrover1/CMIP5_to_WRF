; script to combine the decadal average of all fields from FNL
; wih the CMIP5
; updated:  5/10/16
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

hours = (/"00","06","12","18","00","06","12","18","00","06","12","18","00","06","12","18","00","06","12","18","00"/)
days = (/"14","14","14","14","15","15","15","15","16","16","16","16","17","17","17","17","18","18","18","18","19"/)
;case = (/"minmSL","gfnmSL","ncnmSL"/)
;type = (/"MIROC5","GFDL","NCAR"/)
;type = (/"GFDL"/)
;type2= (/"GFDL"/)
;case = (/"GFDLSL"/)
type = (/"NCAR"/)
type2= (/"NCAR"/)
case = (/"NCARSL"/)

ref_year="2004"
ref_yy="04"
ref_month="09"

nameCMIP = (/"zg","ta","hur","va","ua","hurs","ts","ps","tas","psl",\
             "tsl","mrlsl","pt","pt","pt","pt","vas","uas"/)


nameCDAS = (/"HGT_3_ISBL",\
            "TMP_3_ISBL","R_H_3_ISBL",\
            "V_GRD_3_ISBL","U_GRD_3_ISBL",\
            "R_H_3_HTGL","TMP_3_SFC","PRES_3_SFC",\
            "TMP_3_HTGL",\
            "PRMSL_3_MSL",\
            "TMP_3_DBLY","SOILW_3_DBLY","ICE_C_3_SFC",\
            "WEASD_3_SFC","HGT_3_SFC","LAND_3_SFC",\
            "V_GRD_3_HTGL","U_GRD_3_HTGL"/)
nameC = (/"HGT3D","TMP3D","RH3D","VGRD3D","UGRD3D",\
          "RHsfc","SKINsfc","PRESsfc",\
          "TMPsfc","PSLsfc",\
          "TSOILsfc","SOILWsfc","ICECsfc","WEASDsfc","GRNDsfc","LANDsfc",\
          "VGRDsfc","UGRDsfc"/)

; this is for the variable name in the nc file
nameF = (/"HGT","TMP","RH","VGRD","UGRD",\
          "RH","SKIN","PRES",\
          "TMP","PSL",\
          "TMP","SOILW","ICEC","WEASD","GRND","LAND",\
          "VGRD","UGRD"/)

do nc=0,0

dcase = case(nc)
ctype = type(nc)
ctype2= type2(nc)

do nt=0,20

ref_day=days(nt)
ref_hour=hours(nt)

print("day ="+ref_day+" hour="+ref_hour)

date=ref_month+days(0)+ref_yy
;date = ref_year+ref_month+days(0)+"_"+ref_hour+"_00"

; loop through variables
do nv = 0,17
CMIPname=nameCMIP(nv)
CDASname=nameCDAS(nv)
Cname=nameC(nv)
Fname=nameF(nv)

print("Cname= "+Cname)

;Set Path
main1 = "/data/jtrapp/a/crrllsm2/ivandata/ncepfnl/onedeg_ncepfnl/"
main2 = "./"
main3 = "/data/jtrapp/a/jtrapp/ivan/"+dcase+"dataprep"+date+"/"

; get time/date variables from shell script
file1 = "fnl_"+ref_year+ref_month+ref_day+"_"+ref_hour+"_00.grib1"
f1 = addfile(main1+file1, "r")
;CDASvar  = f1->$CDASname$  ; time, plev, lat, lon

;if(CDASname.eq."V_GRD_3_HTGL".or.\
;                                 CDASname.eq."U_GRD_3_HTGL") then
;CDASvar  = f1->$CDASname$(0,:,:)  ; lat, lon
;else
CDASvar  = f1->$CDASname$  ; lat, lon
;end if


if (CMIPname.eq."pt") then
; pass thru
  Cvar = CDASvar

else if(CMIPname.eq."tsl".or.CMIPname.eq."mrlsl") then
     print("in mrlsl loop: "+CMIPname)
; for the soil variables, multiply instead of add
; not sure why this is needed...
     file2 = CMIPname+"_"+ctype2+"ncepfnldiff.nc"
     f2 = addfile(main2+file2, "r")
     CMIPvar  = f2->$CMIPname$  ;  plev, lat, lon
     Cvar = CDASvar * CMIPvar
;     Cvar = CDASvar
     delete(CMIPvar)
     delete(f2)

else
; combine
printVarSummary(CDASvar)

  if (ctype.ne."NCAR") then
     file2 = CMIPname+"_"+ctype+"ncepfnldiff.nc"
     f2 = addfile(main2+file2, "r")
     CMIPvar  = f2->$CMIPname$  ;  plev, lat, lon
     printVarSummary(CMIPvar)
     Cvar = CDASvar + CMIPvar
     delete(CMIPvar)
     delete(f2)
  end if
; account for lack of surface winds in NCAR
  if (ctype.eq."NCAR") then
     if(CMIPname.eq."uas") then
       Cvar = CDASvar
     else if(CMIPname.eq."vas") then
       Cvar = CDASvar
     else
       file2 = CMIPname+"_"+ctype+"ncepfnldiff.nc"
       f2 = addfile(main2+file2, "r")
       CMIPvar  = f2->$CMIPname$  ;  plev, lat, lon
       Cvar = CDASvar + CMIPvar
       delete(CMIPvar)
       delete(f2)
     end if
     end if
  end if
end if
end if

; try copying from the original cdas to the modified cdas
copy_VarMeta(CDASvar,Cvar)
;copy_VarMeta(CDASvar_regrid,Cvar)

; check for RH > 100
if(Cname.eq."RHsfc".or.Cname.eq."RH3D") then
  Cvar= where(Cvar.gt.100., 100., Cvar)
end if
; check for RH < 0
if(Cname.eq."RHsfc".or.Cname.eq."RH3D") then
  Cvar= where(Cvar.lt.0., 0., Cvar)
end if

;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
; now, dump out the combined field to a
; netcdf file
;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

ncfile= main3+ref_year+ref_month+ref_day+ref_hour+"00.combine."\
          +ref_year+ref_month+ref_day+"."+Cname+".nc"
system ("mkdir " + main3)            ; create dir
system ("rm " + ncfile)               ; remove any pre-existing file
ncdf   = addfile(ncfile ,"c")   ; open output netCDF file

ncdf->$Fname$  = Cvar
delete(CDASvar)
delete(Cvar)
delete(f1)

end do ; end variables
; before leaving this time, compute sfc specific humidity
rhfile= main3+ref_year+ref_month+ref_day+ref_hour+"00.combine."\
         +ref_year+ref_month+ref_day+".RHsfc.nc"
tfile= main3+ref_year+ref_month+ref_day+ref_hour+"00.combine."\
         +ref_year+ref_month+ref_day+".TMPsfc.nc"
pfile= main3+ref_year+ref_month+ref_day+ref_hour+"00.combine."\
         +ref_year+ref_month+ref_day+".PRESsfc.nc"

rhf = addfile(rhfile, "r")
tf = addfile(tfile, "r")
pf = addfile(pfile, "r")
rh  = rhf->RH  ; lat, lon
t  = tf->TMP  ; lat, lon
p  = pf->PRES  ; lat, lon

q    = mixhum_ptrh (p*.01, t,rh , 2 ) ; spec hum (kg/kg)
ncfile= main3+ref_year+ref_month+ref_day+ref_hour+"00.combine."\
         +ref_year+ref_month+ref_day+".SPFHsfc.nc"
system ("rm " + ncfile)               ; remove any pre-existing file
ncdf   = addfile(ncfile ,"c")   ; open output netCDF file
ncdf->SPFH  = q


end do ; end times
end do ; end cases

end
