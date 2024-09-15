---
jupyter:
  jupytext:
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.14.0
  kernelspec:
    display_name: pjrpy3
    language: python
    name: pjrpy3
---

```python
import sys
print(sys.version)
%matplotlib inline
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as mcolors
import numpy as np
import xarray as xr
xr.set_options(keep_attrs=True)
import cartopy.crs as ccrs

# load some useful functions written or acquired by phil rasch
%run -i ./pjrlib_v2
```

```python
# process the UKESM data
```

```python
def fix_UKMO_ds(filename, dsin: xr.Dataset) -> xr.Dataset:
    """
    Rename some variables and rescale 
    :param ds: some xarray dataset
    :return: a normalized xarray dataset, or the original one
    """

    #print('filename is ', filename)
    name_dict = dict()
    ds = dsin.copy()
    ds = _normalize_lat_lon(ds) # harmonize the lat, lon fields
    ds.coords['lon'] = (ds.coords['lon'] + 180) % 360 - 180
    ds = ds.sortby(ds.lon)
    #print('ds',ds)
    for Vname in ds:
        #print('Vname',Vname)
        if (('PBL' in filename) & ('height_after' in Vname)):
            ds[Vname].attrs['UKESM name'] = Vname
            ds[Vname].attrs['units'] = 'kg/m2'
            ds[Vname].attrs['long_name'] = 'PBL height'
            ds = ds.rename({Vname:'PBLH'})
            #print('fixed PBLH')
            
        if (('LWP' in filename) & (Vname == 'm01s02i391')):
            ds[Vname].attrs['UKESM name'] = Vname
            ds[Vname].attrs['units'] = 'kg/m2'
            ds[Vname].attrs['long_name'] = 'Liq Water Path'
            ds = ds.rename({Vname:'TGCLDLWP'})
            #print('fixed LWP')
            
        if (('AOD' in filename) & (Vname == 'unknown')):
            ds[Vname].attrs['UKESM name'] = Vname
            ds[Vname].attrs['units'] = '1'
            ds[Vname].attrs['long_name'] = 'AOD 550nm'
            ds = ds.rename({Vname:'AOD'})
    
        if (('_r_e' in filename) & (Vname == 'unknown')):
            ds[Vname].attrs['UKESM name'] = Vname
            ds[Vname].attrs['long_name'] = 'Effective Radius at Cloud-top'
            ds[Vname].attrs['units'] = '$\mu m$'
            ds = ds.rename({Vname:'REPJR'})
            #print('fixed RE')

        if (('precip_rate' in filename) & (Vname == 'precipitation_flux')):
            ds[Vname].attrs['UKESM name'] = Vname
            ds[Vname].attrs['units'] = 'm/s'
            ds[Vname] = ds[Vname]/1000.
            ds[Vname].attrs['long_name'] = 'Total Precipitation'
            ds = ds.rename({Vname:'PRECT'})
            #print('fixed PRECT')
           
        if (('p_surface_Pa' in filename) & (Vname == 'surface_air_pressure')):
            ds[Vname].attrs['UKESM name'] = Vname
            ds[Vname] = ds[Vname]/100.
            ds[Vname].attrs['units'] = 'hPa'
            ds[Vname].attrs['long_name'] = ds[Vname].attrs['standard_name']
            ds = ds.rename({Vname:'PS'})
            #print('fixed PS')
        
        if (('cloud_fraction' in filename) & (Vname == 'cloud_area_fraction_in_atmosphere_layer')):
            ds[Vname].attrs['UKESM name'] = Vname
            #ds[Vname] = ds[Vname]/100.
            #ds[Vname].attrs['units'] = 'hPa'
            #ds[Vname].attrs['long_name'] = ds[Vname].attrs['standard_name']
            ds = ds.rename({Vname:'CLOUD'})
        
        if (('T_surface_K' in filename) & (Vname == 'surface_temperature')):
            ds[Vname].attrs['UKESM name'] = Vname
            ds[Vname].attrs['long_name'] = ds[Vname].attrs['standard_name']
            ds = ds.rename({Vname:'TS'})
            #print('fixed TS')
            
        if (('net_ToA_SW_W' in filename) & (Vname == 'unknown')):
            ds[Vname].attrs['UKESM name'] = Vname
            ds[Vname].attrs['units'] = 'W/m2'
            ds[Vname].attrs['long_name'] = 'Net TOA Shortwave'
            ds = ds.rename({Vname:'FSNT'})
            #print('fixed FSNT')
            
        if (('net_ToA_SW_C' in filename) & (Vname == 'unknown')):
            ds[Vname].attrs['UKESM name'] = Vname
            ds[Vname].attrs['units'] = 'W/m2'
            ds[Vname].attrs['long_name'] = 'Net TOA Shortwave Clear'
            ds = ds.rename({Vname:'FSNTC'})
            #print('fixed FSNT'
            
        if (('_SW_Clear' in filename) & (Vname == 'toa_outgoing_shortwave_flux_assuming_clear_sky')):
            ds[Vname] = -ds[Vname]
            ds[Vname].attrs['UKESM name'] = Vname
            ds[Vname].attrs['units'] = 'W/m2'
            ds[Vname].attrs['long_name'] = 'outgoing SW assuming clearsky (upward -ive)'
            ds = ds.rename({Vname:'MFSUTC'})
              

    return ds
```

```python
Vdict = {'net_ToA_LW_W_m2':'toa_outgoing_longwave_flux'
        ,'net_ToA_SW_W_m2':'FSNT'
        ,'net_ToA_SW_Clear_W_m2':'FSNTC'
        ,'AOD_550nm':'AOD'
        ,'LWP_kg_m2':'TGCLDLWP'
        ,'p_surface_Pa':'PS'
        ,'T_surface_K':'TS'
        ,'precip_rate_kg_m2_sec':'PRECT'
        ,'PBL_depth_metres':'PBLH'
        ,'cloudtop_r_e_microns':'REPJR'
        ,'cloud_fraction':'CLOUD'
        ,'Outgoing_SW_Clear_W_m2':'MFSUTC'
        }


```

```python
pref_fn = 'UKESM'

def make_ind1 (REG_ID, Varname, filetype=None):
    case_start1 = '~/NetCDF_Files/UKESM1_data/'+REG_ID+'_20450101_20490101_mean_'
    case_start1 = '~/NetCDF_Files/UKESM1_data_v2/Coupled_50Tg/'+REG_ID+'_coupled_50Tgy_20410101_20500101_mean_'
    case_end1 = ".nc"
    fstring1 ='%s%s%s' 
    
    if filetype == 'Fixed_SST':
        # fixed SST simulations
        case_start1 = '~/NetCDF_Files/UKESM1_data/'+REG_ID+'_AtmosOnly_19840101_19880101_mean_'
        pref1=REG_ID+'_50Tgpy_FixedSST'
        case_start1 = '~/NetCDF_Files/UKESM1_data_v2/AtmosOnly_25Tg_1979-1989/'+REG_ID+'_AtmosOnly_25Tgy_19790101_19890101_mean_'
        pref1=REG_ID+'_25Tgpy_FixedSST'
        #case_start1 = '~/NetCDF_Files/UKESM1_data_v2/AtmosOnly_50Tg_1979-1989/'+REG_ID+'_AtmosOnly_50Tgy_19790101_19890101_mean_'
        #pref1=REG_ID+'_50Tgpy_FixedSST'
        case_end1 = ".nc"
        fstring1 ='%s%s%s' 

    ind1 = fstring1 % (case_start1,Varname,case_end1)
    return ind1

def make_ind2(REG_ID, Varname, filetype='Fixed_SST'):    
    case_start2 = '~/NetCDF_Files/UKESM1_data/CTL_20450101_20490101_mean_'
    case_start2 = '~/NetCDF_Files/UKESM1_data_v2/Coupled_Control/CTL_coupled_20410101_20500101_mean_'
    case_end2 = ".nc"
    fstring2 ='%s%s%s' 

    if filetype == 'Fixed_SST':
        case_start2 = '~/NetCDF_Files/UKESM1_data/CTL_AtmosOnly_19840101_19880101_mean_'
        case_start2 = '~/NetCDF_Files/UKESM1_data_v2/AtmosOnly_Control_1979-1989/'+'CTL_AtmosOnly_19790101_19890101_mean_'
        case_end2 = ".nc"
        fstring2 ='%s%s%s' 
        
    ind2 = fstring2 % (case_start2,Varname,case_end2)

    return ind2


```

```python
# accumulate differences over 3 areas

difftitle=""

#Varlist = np.array(['p_surface_Pa'])
#Varlist = np.array(['Outgoing_SW_Clear_W_m2','p_surface_Pa','T_surface_K','precip_rate_kg_m2_sec','PBL_depth_metres','cloudtop_r_e_microns','AOD_550nm','LWP_kg_m2','net_ToA_LW_W_m2','net_ToA_SW_W_m2'])
Varlist = np.array(['Outgoing_SW_Clear_W_m2','precip_rate_kg_m2_sec','PBL_depth_metres','cloudtop_r_e_microns','AOD_550nm','LWP_kg_m2','net_ToA_LW_W_m2','net_ToA_SW_W_m2',
                    'net_ToA_SW_Clear_W_m2','cloud_fraction'])
Varlist = np.array(['net_ToA_SW_Clear_W_m2','net_ToA_SW_W_m2'])
#Varlist = np.array(['net_ToA_SW_Clear_W_m2','net_ToA_SW_W_m2'])
#Varlist = np.array(['AOD_550nm','LWP_kg_m2','net_ToA_LW_W_m2','net_ToA_SW_W_m2'])
#Varlist = np.array(['T_surface_K'])
Varlist = np.array(['LWP_kg_m2','net_ToA_SW_Clear_W_m2','net_ToA_SW_W_m2','cloud_fraction'])
#Varlist = np.array(['precip_rate_kg_m2_sec'])

FSNT1 = None
FSNT2 = None
FSNTC1 = None
FSNTC2 = None

# specify regions (assume lon always specified as west, then east limit)
xreg = np.array([[-150.,-110.],[-110,-70],[-25.,15.],[170.,-120.],[-170.,-90.]])%360.
yreg = np.array([[0.,30.],     [-30.,0.], [-30.,0.], [30.,50.],   [-50.,-30.] ])
namereg = ['NEP','SEP','SEA','NP','SP']
#xreg = [[0.,360.]]
#yreg = [[-90.,91.]]


reglist = np.array(['R1_NEP','R2_SEP','R3_SEA'])
#reglist = np.array(['R2_SEP'])


filetype = None
filetype = 'Coupled'
filetype = 'Fixed_SST'

for Varname in Varlist:
    print()
    print('-------------------------------'+Varname)
    nreg = 0 # is it the start of the region summation
    for REG_ID in reglist:
        ind1 = make_ind1(REG_ID,Varname,filetype)
        print('ind1 opening',ind1)
        DS1 = xr.open_mfdataset(ind1)
        # update the dataset
        DS1 = fix_UKMO_ds(ind1, DS1)
        VN = Vdict[Varname]
        V1 = xr_getvar_sl(VN,DS1,method='maxb850')
        ind2 = make_ind2(REG_ID,Varname,filetype)
        print('ind2 opening',ind2)
        DS2 = xr.open_mfdataset(ind2)
        DS2 = fix_UKMO_ds(ind2, DS2)
        V2 = xr_getvar_sl(VN,DS2,method='maxb850')

        DV = V1-V2
        weights = None
        if 'area' in DS1:
            area = DS1['area']
        elif 'area' in DS2:
            area = DS2['area']
        else:
            lat = V1['lat'].values
            lon = V1['lon'].values
            area = make_fvarea(lon,lat)
        weights = V1.copy()
        weights.data =area
        weights.attrs['units']='steradians'

        V1A = V1.weighted(weights).mean()
        sV1A = ' (%5.2f)' % V1A
        V2A = V2.weighted(weights).mean()
        sV2A = ' (%5.2f)' % V2A
        DVA = V1A-V2A
        sDVA = ' (%5.2f)' % DVA
        if nreg == 0:
            V1S = V1
            V2S = V2
            DVS = DV
            nreg = 1
        else:
            V1S = V1S + V1
            V2S = V2S +V2
            DVS = DVS + DV
            nreg = nreg+1
    
    print(' all regions summed')
    #V1S = V1S/nreg
    #V2S = V2S/nreg

    DVA = DVS.weighted(weights).mean()
    sDVA = ' ({:5.2f}{:s})'.format(DVA.values, DVA.units)

    fname = pref_fn+'_'+Varname+'_'+DV.name+'-D.pdf'

    if VN == 'FSNT':
        FSNT1=V1S
        FSNT2=V2S
        UDFSNT = DVS
    elif VN == 'FSNTC':
        FSNTC1=V1S
        FSNTC2=V2S
        UDFSNTC = DVS
    elif VN == 'TGCLDLWP':
        UDLWP = DVS
    elif VN == 'CLOUD':
        UDCLDLOW = DVS
    
    print('field processing complete',DVS.long_name)

uwts = weights

if ((FSNT1 is None) or (FSNTC1 is None)):
    print ('fields for SWCRE not requested')
else:
    SWCRE1 = FSNT1-FSNTC1
    SWCRE1.attrs['long_name'] = 'TOA shortwave CRE'
    SWCRE1 = SWCRE1.rename('SWCRE')
    SWCRE2 = FSNT2-FSNTC2
    SWCRE2.attrs['long_name'] = 'TOA shortwave CRE'
    SWCRE2 = SWCRE2.rename('SWCRE')
    DSWCRE = SWCRE1-SWCRE2
    DSWCREA = DSWCRE.weighted(weights).mean()
    sDSWCREA = ' ({:5.2f}{:s})'.format(DSWCREA.values, DSWCREA.units)
    fname = pref_fn+'_'+Varname+'_'+'SWCRE'+'-D.pdf'
    #pltfld(DSWCRE, difftitle+sDSWCREA,fname)

```

```python
# now do the CESM and E3SM
```

```python
def makeind(modname,Varname):
    if modname == 'E3SM':
        case_start1 = "~/data_for_protocol_paper/E3SMv2/20230426.v2.LR.F2010.MCB-SSLT-EM.R1-3.test01_ANN_000101_001112_climo_fv192x288.nc"
        case_end1 = ""
        pref1='E3SM_50Tgpy'
        fstring1 ='%s%.0s%.0s' 
        #fstring1 ='%s%s%s' 

        case_start2 = "~/data_for_protocol_paper/E3SMv2/20220930.v2.LR.F2010.E1_CNTL.eam.h0.1-14."
        case_end2 = ".nc"
        fstring2 ='%s%.0s%.0s' 
        fstring2 ='%s%s%s' 
        pref2='E3SMcontrol'
    else:
        case_start1 = "~/data_for_protocol_paper/CESM2/F2010climo.ss_NEP_SEP_SEA.1.5Tg.cam.h0."
        case_end1 = ".1-25.nc"
        pref1='CESM_7.5Tgpyr'
        fstring1 ='%s%.0s%.0s' 
        fstring1 ='%s%s%s' 

        case_start2 = "~/data_for_protocol_paper/CESM2/Fixed_SST.cam.h0.1-20."
        case_end2 = ".nc"
        fstring2 ='%s%.0s%.0s' 
        fstring2 ='%s%s%s' 
        pref2='CESMcontrol'

    ind1 = fstring1 % (case_start1,Varname,case_end1)
    ind2 = fstring2 % (case_start2,Varname,case_end2)
    return ind1, ind2

regtag = ""
weights = None

Varlist = np.array(['CLDTOT','TGCLDLWP','FSNT','FSNTC'])


# specify regions (assume lon always specified as west, then east limit)
xreg = np.array([[-150.,-110.],[-110,-70],[-25.,15.],[170.,-120.],[-170.,-90.]])%360.
yreg = np.array([[0.,30.],     [-30.,0.], [-30.,0.], [30.,50.],   [-50.,-30.] ])
namereg = ['NEP','SEP','SEA','NP','SP']

for pref1 in ['CESM','E3SM']:
    print ('processing',pref1)
    
    for Varname in Varlist:
        print()
        print('-------------------------------'+Varname) 

        ind1, ind2 = makeind (pref1,Varname)
        print('opening ind1 ',ind1)
        DS1 = xr.open_mfdataset(ind1)
        #print('xxx',DS1.time)
        DS1 = center_time(DS1)
        Var1 = xr_getvar(Varname,DS1)
        V1 = Var1.mean(dim='time',keep_attrs=True)

        print('opening ind2',ind2)
        #DS2 = xr.open_mfdataset(ind2)
        DS2 = xr.open_mfdataset(ind2)

        #DS2 = center_time(DS2)
        Var2 = xr_getvar(Varname,DS2)

        V2 = Var2.mean(dim='time',keep_attrs=True)

        DV = V1-V2
        
        if Varname == 'CLDTOT':
            if pref1 == 'E3SM':
                EDCLDTOT = DV
            else:
                CDCLDTOT = DV
        elif Varname == 'CLDLOW':
            if pref1 == 'E3SM':
                EDCLDLOW = DV
            else:
                CDCLDLOW = DV

        elif Varname == 'TGCLDLWP':
            if pref1 == 'E3SM':
                EDLWP = DV
            else:
                CDLWP = DV
        elif Varname == 'FSNT':
            if pref1 == 'E3SM':
                EDFSNT = DV
            else:
                CDFSNT = DV
        elif Varname == 'FSNTC':
            if pref1 == 'E3SM':
                EDFSNTC = DV
            else:
                CDFSNTC = DV
                
        
        print('DV range', DV.min().values, DV.max().values)
        if 'area' in DS1:
            area = DS1['area']
        elif 'area' in DS2:
            area = DS2['area']
        else:
            print('calculating weights')
            lat = Var1['lat'].values
            lon = Var1['lon'].values
            area = make_fvarea(lon,lat)
        weights = V1.copy()
        weights.data =area
        weights.attrs['units']='steradians'

        print(Varname, V1.attrs['long_name'],'Range V1 and V2 ',V1.min().values, V1.max().values, V2.min().values, V2.max().values)
        V1A = V1.weighted(weights).mean()
        sV1A = ' (%5.2f)' % V1A
        V2A = V2.weighted(weights).mean()
        sV2A = ' (%5.2f)' % V2A
        DVA = V1A-V2A
        #sDVA = ' (%5.2f)' % DVA
        sDVA = ' ({:5.2f}{:s})'.format(DVA.values, DVA.units)

        pref2 = 'pref2'
        print('area avgs '+pref1+' %5.2f' % (V1A.values),' '+pref2+' %5.2f' % (V2A.values),' Delta %5.2f' % (DVA.values))



        print('field processing complete for', pref1)


```

```python
def s3load_zarr(s3path):
    """ return an xarray DataSet given an s3path """
    return xr.open_dataset(fsspec.get_mapper(s3path),engine = 'zarr')

```

```python
def make_fvarea(lon,lat,alt=None):
    """make_fvarea(lon,lat)

    make finite volume areas, given 1-D lon and lat arrays
    """
    pio180 = np.pi/180.
    yw = np.zeros(len(lat)+1)
    # walls in degrees
    yw[1:-1] = (lat[1:] + lat[:-1]) / 2
    if (np.abs(lat[-1]-lat[0])-180.) < 1.e-4:
        #print('polepoint')
        yw[0] = lat[0]
        yw[-1] = lat[-1]
    else:
        #print('nopole point')
        yw[0] = lat[0]-(lat[1]-lat[0])/2.
        yw[-1] = lat[-1]+(lat[1]-lat[0])/2.

    if alt == True:
        #print('alt True')
        dy = np.diff(yw*pio180)*np.cos(lat*pio180)
    else:
        # wall in sin(latitude in radians)
        yw = np.sin(yw*pio180)
        dy = np.diff(yw)
    dx = float(lon[1]-lon[0])*pio180
    dxa = np.full([len(lat),len(lon)],dx)
    area = dxa*dy[:,np.newaxis]

    return area
    
def s3load_zarr(s3path):
    """ return an xarray DataSet given an s3path """
    return xr.open_dataset(fsspec.get_mapper(s3path),engine = 'zarr')

def center_time(DS1):
    """center_time(DS1)                                                                                                  
                                                                                                                         
    correct the time coordinate in DS1 to represent the center of the time bounds
    if the bounds variable is missing this only works for monthly data
                                                                                                                         
    """
    DS = DS1.copy()
    # the time coord is often registered at the end of the time averaging interval                                       
    # reset to the midpoint of the interval                                                                              
    time = DS['time'].copy()
    #print('center_time time',time)
    #print('xxx',time.values)                                                                                            
    bndname = time.attrs['bounds']
    #print('bndname',bndname)
    # check whether bndname is present
    #print('yyy',DS1)
    if bndname not in DS1:
        #print('need to improvise')
        DS['time'] = DS.indexes['time'].shift(-1,"MS")
        month_length = DS.time.dt.days_in_month.values
        month_half = month_length/2
        dm = month_half.astype(int)
        hm = ((month_half - dm)*24.).astype(int)
        DS['time'] = DS.indexes['time'].shift(dm,"D")
        DS['time'] = DS.indexes['time'].shift(hm,"h")
    else:
        time_bnds = DS1[bndname]
        #print('time_bnds',time_bnds)                                                                                        
        tbdims = time_bnds.dims
        #print('tbdims',tbdims)                                                                                              
        tbd_name = ''
        # find the bounds name (the dim that isn't time)                                                                     
        for tbd in tbdims:
            if tbd != 'time':
                tbd_name = tbd
        #print('tbd_name',tbd_name)                                                                                          
        #print('tbdims',tbdims)                                                                                              
        # if no bounds, then do nothing                                                                                      
        if tbd_name != '':
            #tb = time_bnds.values                                                                                           
            #print('time_bnds',time_bnds)                                                                                    
            tbm = time_bnds.mean(dim=tbd_name).values
            #print('yyy',tbm)
            #1./0.
            DS.coords["time"] = tbm
            DS['time'].attrs['long_name'] = 'time'
            DS['time'].attrs['bounds'] = 'time_bnds'
    #time = DS['time'].values
    #print('center_time returning DS', time[0],time[1],time[-2],time[-1])
    return DS

def fix_area(DS1):
    """fix_area(DS1)                                                                                                  
                                                                                                                         
    add an area variable to DS1 if absent  
    only works on E3SM grids for ne30pg2 grid
                                                                                                                         
    """
    DS = DS1.copy()
# check for area variable. If it is missing add it (differs for E3SM and CESM grids)
    if 'area' not in DS1:
        if 'ncol' in DS1.dims:
            #print ('get CS data')
            #areafile = '~/NetCDF_Files/F2010_PJR1.eam.h0.0001-01.nc'
            areafile = '~jupyter-adminphil/NetCDF_Files/ne30pg2.nc'
            DSA = xr.open_mfdataset(areafile)
            if len(DSA['ncol']) != len(DS1['ncol']):
                raise ValueError('area file mismatch')
            area = DSA.area
            #lon = DSA.lon
            #lat = DSA.lat
        else:
            #print('calculating fv area weights')
            lat = DS1['lat'].values
            lon = DS1['lon'].values
            aread = make_fvarea(lon,lat)
            area = xr.DataArray(aread, dims=['lat','lon'], coords={'lon':lon,'lat':lat})
            area.attrs['units']='steradians'
        DS['area'] = area
        #DS['lon'] = lon
        #DS['lat'] = lat

    #print('area',area)

    #print('Var1',Var1)
    #1./0.
    return DS

def fix_DS(DS):
    """fix_DS(DS)                                                                                                  
    center the time coordinate within the averaging interval                                                                                                                     
    add an area variable to DS if absent                                      
                                                                                                                         
    """
    DS = fix_area(DS)
    #print('after area ', DS)
    DS = center_time(DS)
    #print('after time ', DS)
    return DS

```

```python
var = 'CLDLOW'
experiment = "20220930.v2.LR.F2010.E1_CNTL"  # control
s3_path = f's3://mcb-zarr/E3SMv2/{experiment}/atm/proc/tseries/month_1/{experiment}.eam.h0.{var}.zarr'
#E3SMv2/20220930.v2.LR.F2010.E1_CNTL/atm/proc/tseries/month_1/20220930.v2.LR.F2010.E1_CNTL.eam.h0.CLDLOW
#E3SMv2/20230330.v2.LR.F2010.MCB-SSLT-EM.R1-3.test01/atm/proc/tseries/month_1/20230330.v2.LR.F2010.MCB-SSLT-EM.R1-3.test01.eam.h0.CLDLOW
EDSc = s3load_zarr(s3_path)
#print('xxx',DS1)
#print('time', DS1.time)
EDSc = fix_DS(EDSc)
#print('control times',EDSc.time)
experiment = '20230426.v2.LR.F2010.MCB-SSLT-EM.R1-3.test01' # experiment
s3_path = f's3://mcb-zarr/E3SMv2/{experiment}/atm/proc/tseries/month_1/{experiment}.eam.h0.{var}.zarr'
EDSe = s3load_zarr(s3_path)
EDSe = fix_DS(EDSe)
#print('experiment times', EDSe.time)
```

```python
#experiment = "b.e21.BSSP245smbb_MCB600cm_R1R2R3.f09_g17.LE2-1011.001"
#CESM2/F2010climo.ss_NEP_SEP_SEA.1.5Tg/atm/proc/tseries/month_1/F2010climo.ss_NEP_SEP_SEA.1.5Tg.cam.h0.CLDLOW.1-25
experiment = "F2010climo.ss_NEP_SEP_SEA.1.5Tg"
s3_path = f's3://mcb-zarr/CESM2/{experiment}/atm/proc/tseries/month_1/{experiment}.cam.h0.{var}.1-25.zarr'
#print('s3_path',s3_path)
CDSe = s3load_zarr(s3_path)
CDSe = fix_DS(CDSe)
print('CESM DSe', CDSe)
# CESM2/F2010climo/atm/proc/tseries/month_1/F2010climo.cam.h0.CLDLOW.1-19
experiment = "F2010climo"
s3_path = f's3://mcb-zarr/CESM2/{experiment}/atm/proc/tseries/month_1/{experiment}.cam.h0.{var}.1-19.zarr'
#print('s3_path',s3_path)
CDSc = s3load_zarr(s3_path)
CDSc = fix_DS(CDSc)
print('CESM DSc', CDSc)
```

```python
def make_AA_tser(Var):
    ''' make Annual Average from tseries files
        Var: Xarray object
        returns xarray object with annual averages
    '''
    month_length = Var.time.dt.days_in_month
    #twgts = month_length.groupby("time.year") / month_length.groupby("time.year").sum()
    V1AY = (Var*month_length).groupby("time.year").sum()/month_length.groupby("time.year").sum()
    #print('V1AY.values', V1AY.values)
    # for some reason, the groupby operator does not propogate the variable name
    V1AY = V1AY.rename(Var.name)
    return V1AY

```

```python
CCLDLOWc = make_AA_tser(CDSc[var]).mean(dim='year')
CCLDLOWe = make_AA_tser(CDSe[var]).mean(dim='year')
DCCLDLOW = CCLDLOWe-CCLDLOWc
DCCLDLOW = DCCLDLOW*100.
DCCLDLOW.attrs['units'] = '%'
```

```python
ECLDLOWc = make_AA_tser(EDSc[var]).mean(dim='year')
ECLDLOWe = make_AA_tser(EDSe[var]).mean(dim='year')
DECLDLOW = ECLDLOWe-ECLDLOWc

# convert from cubed sphere to the CESM lat/lon grid
xoutm,youtm=np.meshgrid(DCCLDLOW.lon.values,DCCLDLOW.lat.values)
Vnew1 = xr.DataArray(interp_ap(xoutm, youtm, DECLDLOW.values,EDSc.lat.values,EDSc.lon.values), 
                    coords={'lat': DCCLDLOW.lat.values,'lon': DCCLDLOW.lon.values}, #,'lev': Var_e3sm.lev.values},
                    attrs=DECLDLOW.attrs,
                    dims=["lat", "lon"]) #,"lev"])
DECLDLOW = Vnew1
DECLDLOW = DECLDLOW*100.
DECLDLOW.attrs['units'] = '%'
```

```python
def xr_llhplot2 (xrVar, cbar='default', plotproj=None, ax=None, cax=None,
                 ylabels=None, clevs=None, cmap=None, title=None, cbartitle=None,regmark=False):
    """xr_llhplot xarray lat lon horizontal plot
       returns a "mappable" containing the artist info that is needed to plot a colorbar
       that is
       mpbl = xr_llhplot2()
       plt.colorbar(mpbl, orientation='horizontal',cax=cax,...)
    """
    #print(' entering xr_llhplot', xrVar)
    
    lon=xrVar['lon'].values
    lat=xrVar['lat'].values
    xv,yv=np.meshgrid(lon,lat)
    data_regridded = xrVar.values
    #print('aaa',data_regridded.shape, xv.shape, yv.shape)
    df = data_regridded.flatten()
    dsub = df[np.isfinite(df)] # ignore NaN
    zmax = dsub.max()
    zmin = dsub.min()
    #print('masked interpolated range',zmin,zmax)
    dataproj=ccrs.PlateCarree()    # data is always assumed to be lat/lon
    if ylabels is None: ylabels = True
    if clevs is None:
        clevs = findNiceContours(np.array([zmin,zmax]),nlevs=10)
    #print('clevs',clevs)
    if cmap is None:
        #print('aaa, grabbing cmap default')
        #cmap = mpl.cm.get_cmap()
        cmap = plt.get_cmap()
        #print('bbb',cmap.N)
    #print('cmap',cmap)
    extend = 'both'
    norm = mpl.colors.BoundaryNorm(clevs,cmap.N,extend=extend)
    #print('norm',norm(clevs))
    clat = (lat.min()+lat.max())/2.
    clon = (lon.min()+lon.max())/2.
    if plotproj is None:
        plotproj = ccrs.PlateCarree()
        plotproj = ccrs.Mollweide()
 
    # if no ax argument, could get current axis, or create it
    if ax is None:
        #print('grab current axis')
        #ax = plt.gca()
        ax = plt.axes(projection=plotproj)
        
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,
                      linewidth=2, color='gray', alpha=0.5)
    pl = ax.contourf(xv, yv, data_regridded, levels=clevs, # vmin=zmin, vmax=zmax,
                     norm=norm, cmap=cmap,
                     extend=extend, transform=ccrs.PlateCarree())
    
    gl.left_labels=ylabels
    gl.right_labels=ylabels
    ax.coastlines(linewidth=1,color='blue')
 
    ## Find the location of the main plot axes
    ## has to be done after some plotting is done in projection space
    posn = ax.get_position()
    
    fig = plt.gcf()
    ax2 = fig.add_axes([0,0,0.1,0.1])
    ax2.set_position([posn.x0-0.005, posn.y0-0.005, posn.width+0.01, posn.height+0.01])
    ax2.patch.set_alpha(0.0)
    ax2.set_axis_off()
    ax2.set_xlim([0,1])
    ax2.set_ylim([0,1])
    if regmark:
        # print some registration marks to help in lining up figures
        ax2.scatter([0,0,1,1], [0,1,0,1], c="r", s=100)

    if not title is None:
        #print('title is ', title)
        #ax2.set_title(title)
        ax2.text(0.01,0.88,title,fontsize=6)

        
    if cbar == 'default':
        # Add colorbar to plot
        if cbartitle is None:
            cbartitle = xrVar.long_name
            
        if cax is not None:
            cax = ax
        else:
            # create an colorbar axis
            cax = fig.add_axes([0,0,0.1,0.1])
            ## Adjust the positioning and orientation of the colorbar
            #ax.set_position([posn.x0, posn.y0-0.06, posn.width, 0.04])
            cax.set_position([posn.x0, posn.y0-0.02, posn.width, 0.015])
        
        units = xrVar.units
        if units == 'mm/day':
            units = 'mm day$^{-1}$'
            
        cb = plt.colorbar(
             pl, orientation='horizontal',ticks=clevs,cax=cax,
             label='%s (%s)'%(cbartitle, units)
             )
        cb.ax.tick_params(labelsize=7)
    
        
    return pl
```

```python
# plot all models

def pltllbox(xri, yri,ax=None):
    if ax is None:
        ax = plt.gca()
    if xri[1] < xri[0]:
        xri[1] += 360.
    regcx = [xri[0],xri[1],xri[1],xri[0],xri[0]]
    regcy = [yri[0],yri[0],yri[1],yri[1],yri[0]]
    ax.plot(regcx,regcy,color='red',transform=ccrs.PlateCarree())
    
def setfign ():
    """
    return fig and axes for a single panel figure
    """
    plotproj = ccrs.Mollweide(central_longitude=-80)
    plotproj._threshold /= 100.
    fig, axes = plt.subplots(ncols=3,nrows=4,
                             #gridspec_kw={'width_ratios': [1]},
                             subplot_kw={'projection': plotproj},
                             figsize=(8.,9.0),
                            )
    plt.subplots_adjust(left=None, bottom=0.2, right=None, top=None, wspace=0.02, hspace=0.05)

    fig.set_dpi(300.0)
    return fig, axes;

def plotall (fig, ax, row, coltitles, UD, ED, CD, dmap=None):# (DUKESM, DSUKESM, DE3SM, DSE3SM, DCESM, DSCESM,dmap=None):
       
    if row == 0:
        ax[0,0].set_title('UKESM1')
        ax[0,1].set_title('E3SMv2')
        ax[0,2].set_title('CESM2')

               
    lat = UD['lat'].values
    lon = UD['lon'].values
    area = make_fvarea(lon,lat)
    weights = UD.copy()
    weights.data =area
    weights.attrs['units']='1'
    wdims = weights.dims
    weights = weights/(weights.sum(dim=wdims))
    sDSU=' GA:%5.2f' % UD.weighted(weights).mean().values

    drmin = np.min([UD.min(),ED.min(),CD.min()]) #, DUKESM.min(),DSCESM.min(),DCESM.min(),DSE3SM.min(),DE3SM.min()])
    drmax = np.max([UD.max(),ED.max(),CD.max()]) #, DUKESM.max(),DSCESM.max(),DCESM.max(),DSE3SM.max(),DE3SM.max()])
    #print('plotting range', drmin, drmax)
    factor = 0.7
    dlevs = findNiceContours(np.array([drmin,drmax])*factor,nlevs = 15,rmClev=0.,sym=True)
    #print('dlevs',dlevs)

    if dmap is None:
        dmap = diverge_map()

    xr_llhplot2(UD, ax=ax[row,0],clevs=dlevs,cmap=dmap,title=sDSU, ylabels=False,cbar=None)

    pltllbox([-150.,-110.],[0.,30.],ax=ax[row,0])
    pltllbox([-110.,-70.],[-30.,0.],ax=ax[row,0])
    pltllbox([-25.,15.],[-30.,0.],ax=ax[row,0])

    lat = ED['lat'].values
    lon = ED['lon'].values
    area = make_fvarea(lon,lat)
    weights = ED.copy()
    weights.data =area
    weights.attrs['units']='1'
    wdims = weights.dims
    weights = weights/(weights.sum(dim=wdims))

    sDSE='GA:%5.2f' % ED.weighted(weights).mean().values
    #print('sDSE', sDSE)

    pl = xr_llhplot2(ED, ax=ax[row,1],clevs=dlevs,cmap=dmap,title=sDSE, ylabels=False,cbar=None)
    pltllbox([-150.,-110.],[0.,30.],ax=ax[row,1])
    pltllbox([-110.,-70.],[-30.,0.],ax=ax[row,1])
    pltllbox([-25.,15.],[-30.,0.],ax=ax[row,1])

    lat = CD['lat'].values
    lon = CD['lon'].values
    area = make_fvarea(lon,lat)
    weights = CD.copy()
    weights.data =area
    weights.attrs['units']='1'
    wdims = weights.dims
    weights = weights/(weights.sum(dim=wdims))

    sDSC='GA:%5.2f' % CD.weighted(weights).mean().values
    #print('sDSE', sDSC)

    xr_llhplot2(CD, ax=ax[row,2],clevs=dlevs,cmap=dmap,title=sDSC, ylabels=False,cbar=None)
    pltllbox([-150.,-110.],[0.,30.],ax=ax[row,2])
    pltllbox([-110.,-70.],[-30.,0.],ax=ax[row,2])
    pltllbox([-25.,15.],[-30.,0.],ax=ax[row,2])

    posn = ax[row,1].get_position()
    ax2 = fig.add_axes([0,0,0.1,0.1])
    ax2.set_position([posn.x0-0.005, posn.y0-0.005, posn.width+0.01, posn.height+0.01])
    ax2.patch.set_alpha(0.0)
    ax2.set_axis_off()
    ax2.set_xlim([0,1])
    ax2.set_ylim([0,1])
    if True:
        ## Find the location of the main plot axes
        ## has to be done after some plotting is done in projection space
        posn = ax[row,1].get_position()
        # create an colorbar axis
        cax = fig.add_axes([0,0,0.1,0.1])
        ## Adjust the positioning and orientation of the colorbar
        xof = 0.2
        cax.set_position([posn.x0-xof, posn.y0-0.015, posn.width+2*xof, 0.011])
        #cax.set_position([posn.x0, posn.y0-0.02, posn.width, 0.015])
        #cax.set_ylim([-1,1])
        if True:
            cb = plt.colorbar(pl, orientation='horizontal',ticks=dlevs,cax=cax)
                 #label='%s (%s)'%(cbartitle, xrVar.units)
            cb.ax.tick_params(labelsize=7)
        units = ED.units
        if units == 'mm/day':
            units = 'mm day$^{-1}$'
        long_name = ED.long_name # 'longname'
        cax.text(0.5,-2.6,'%s (%s)'%(long_name, units),va='center',ha='center')
    
tcmap = diverge_map()
pcmap = plt.get_cmap('BrBG')

fig, axes = setfign()
axf = axes.flatten()
    
plotall (fig, axes, 0, 'ROW FSNT', UDFSNT, EDFSNT,CDFSNT) # (D1245_UKESM_Pr, DS1245_UKESM_Pr, D1245_E3SM_Pr, DS1245_E3SM_Pr, D1245_CESM_Pr, DS1245_CESM_Pr,dmap=pcmap)
plotall (fig, axes, 1, 'ROW FSNT', UDFSNTC, EDFSNTC,CDFSNTC) # (D1245_UKESM_Pr, DS1245_UKESM_Pr, D1245_E3SM_Pr, DS1245_E3SM_Pr, D1245_CESM_Pr, DS1245_CESM_Pr,dmap=pcmap)
#plotall (fig, axes, 2, 'ROW FSNT', UDCLDLOW, EDCLDTOT,CDCLDTOT) # (D1245_UKESM_Pr, DS1245_UKESM_Pr, D1245_E3SM_Pr, DS1245_E3SM_Pr, D1245_CESM_Pr, DS1245_CESM_Pr,dmap=pcmap)
plotall (fig, axes, 2, 'ROW FSNT', UDCLDLOW, DECLDLOW, DCCLDLOW) # (D1245_UKESM_Pr, DS1245_UKESM_Pr, D1245_E3SM_Pr, DS1245_E3SM_Pr, D1245_CESM_Pr, DS1245_CESM_Pr,dmap=pcmap)
plotall (fig, axes, 3, 'ROW FSNT', UDLWP, EDLWP,CDLWP) # (D1245_UKESM_Pr, DS1245_UKESM_Pr, D1245_E3SM_Pr, DS1245_E3SM_Pr, D1245_CESM_Pr, DS1245_CESM_Pr,dmap=pcmap)
plt.savefig('fig5_rev2'+'.pdf',format='pdf',dpi=300,transparent=True)#,facecolor='xkcd:mint green')
plt.show()

#
```

```python

```
