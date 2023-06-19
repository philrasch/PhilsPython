---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.14.5
  kernelspec:
    display_name: Python [conda env:pjrpy3] *
    language: python
    name: conda-env-pjrpy3-py
---

**compare two cases over the globe assuming they are on lat/lon grid at same resolution**

*configured for the GMD paper*

```python
import sys
print(sys.version)
%matplotlib inline
%run -i ~/Python/pjr3
```

```python
def setfig3b1x1 ():
    """
    return fig and axes for a single panel figure
    """
    plotproj = ccrs.Mollweide()
    plotproj._threshold /= 100.
    fig, axes = plt.subplots(ncols=1,
                             gridspec_kw={'width_ratios': [1]},
                             subplot_kw={'projection': plotproj},
                             figsize=(6,3),
                            )
    fig.set_dpi(300.0)
    return fig, axes;


def pltllbox(xri, yri):
    if xri[1] < xri[0]:
        xri[1] += 360.
    regcx = [xri[0],xri[1],xri[1],xri[0],xri[0]]
    regcy = [yri[0],yri[0],yri[1],yri[1],yri[0]]
    plt.plot(regcx,regcy,color='red',transform=ccrs.PlateCarree())
    
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
            ds = ds.rename({Vname:'CFPJR'})
        
        if (('T_surface_K' in filename) & (Vname == 'surface_temperature')):
            ds[Vname].attrs['UKESM name'] = Vname
            ds[Vname].attrs['long_name'] = ds[Vname].attrs['standard_name']
            ds = ds.rename({Vname:'TS'})
            #print('fixed TS')
            
        if (('_SW_W' in filename) & (Vname == 'unknown')):
            ds[Vname].attrs['UKESM name'] = Vname
            ds[Vname].attrs['units'] = 'W/m2'
            ds[Vname].attrs['long_name'] = 'Net TOA Shortwave'
            ds = ds.rename({Vname:'FSNT'})
            #print('fixed FSNT')
            
        if (('_SW_Clear' in filename) & (Vname == 'toa_outgoing_shortwave_flux_assuming_clear_sky')):
            ds[Vname] = -ds[Vname]
            ds[Vname].attrs['UKESM name'] = Vname
            ds[Vname].attrs['units'] = 'W/m2'
            ds[Vname].attrs['long_name'] = 'outgoing SW assuming clearsky (upward -ive)'
            ds = ds.rename({Vname:'MFSUTC'})
              

    return ds


```

```python
def xr_getvar_sl(VN, DS1, method='surface', verbose=False):
    """ get a field from a netCDF file.
    If it is a multi-level field, do something useful to provide single level information
    This function assumes that the x-xoordinate increases monotonically

    """
    Var1 = xr_getvar(VN,DS1)
    #V1 = Var1.mean(dim='time',keep_attrs=True)
    dimlist = Var1.dims
    #print('dimlist',dimlist)
    if 'model_level_number' in dimlist:
        level_height = xr_getvar('level_height',DS1)
        sigma = xr_getvar('sigma',DS1)
        print('sigma',sigma.values)
        surface_altitude = xr_getvar('surface_altitude',DS1)
        altitude = level_height + (sigma * surface_altitude)
        altitude.attrs['long_name'] = 'altitude above mean sea-level'
        #print('altitude',altitude)
        if method == 'surface':
            print('method:surface')
            V1 = Var1.isel(model_level_number=0)
            V1.attrs['long_name'] = V1.attrs['long_name'] + ' (surface level)'
        elif method == 'maxb850':
            i = 1
            j = 1
            print('method:max below 850')
            dz1 = altitude - surface_altitude   # height above surface
            print('altitude',altitude[:,i,j].values,altitude)
            print('surface_altitude',surface_altitude[i,j].values,surface_altitude)
            print('dz1(i,j)',dz1[:,i,j].values)
            dz2 = (sigma - 1)*surface_altitude + level_height  # height above sea level
            print('dz2(i,j)',dz2[:,i,j].values)
            pmb = 850.
            psmb = 1000.
            scaleheight = 8.4e3
            altmb = -np.log(pmb/psmb)*scaleheight
            print('altmb is ', altmb)
            V1 = Var1.copy()
            V1.values = dz1.values ### for the moment, overwrite field with height above surface height
            print('V1',V1[:,i,j].values)
            V2 = V1.where(dz1 <= altmb+50.)
            V3 = V2.max(dim='model_level_number')
            V3.attrs['long_name'] = 'max value below pmb of '+V2.attrs['long_name']
            V1 = V3
            
    else:
        V1 = Var1
    return V1

weights = None

Varlist = np.array(['LWP_kg_m2'])
#Varlist = np.array(['p_surface_Pa'])
#Varlist = np.array(['T_surface_K'])
#Varlist = np.array(['AOD_550nm','LWP_kg_m2','net_ToA_LW_W_m2','net_ToA_SW_W_m2'])
#Varlist = np.array(['Outgoing_SW_Clear_W_m2','p_surface_Pa','T_surface_K','precip_rate_kg_m2_sec','PBL_depth_metres','cloudtop_r_e_microns','AOD_550nm','LWP_kg_m2','net_ToA_LW_W_m2','net_ToA_SW_W_m2'])
Varlist = np.array(['Outgoing_SW_Clear_W_m2','precip_rate_kg_m2_sec','PBL_depth_metres','cloudtop_r_e_microns','AOD_550nm','LWP_kg_m2','net_ToA_LW_W_m2','net_ToA_SW_W_m2'])
Varlist = np.array(['cloud_fraction'])
Varlist = np.array(['Outgoing_SW_Clear_W_m2','net_ToA_SW_W_m2'])
                   
Vdict = {'net_ToA_LW_W_m2':'toa_outgoing_longwave_flux'
        ,'net_ToA_SW_W_m2':'FSNT'
        ,'AOD_550nm':'AOD'
        ,'LWP_kg_m2':'TGCLDLWP'
        ,'p_surface_Pa':'PS'
        ,'T_surface_K':'TS'
        ,'precip_rate_kg_m2_sec':'PRECT'
        ,'PBL_depth_metres':'PBLH'
        ,'cloudtop_r_e_microns':'REPJR'
        ,'cloud_fraction':'CFPJR'
        ,'Outgoing_SW_Clear_W_m2':'MFSUTC'
        }


# specify regions (assume lon always specified as west, then east limit)
xreg = np.array([[-150.,-110.],[-110,-70],[-25.,15.],[170.,-120.],[-170.,-90.]])%360.
yreg = np.array([[0.,30.],     [-30.,0.], [-30.,0.], [30.,50.],   [-50.,-30.] ])
namereg = ['NEP','SEP','SEA','NP','SP']
#xreg = [[0.,360.]]
#yreg = [[-90.,91.]]


REG_ID = 'R1_NEP'
REG_ID = 'R2_SEP'
REG_ID = 'R3_SEA'

# coupled simulations
case_start1 = '~/NetCDF_Files/UKESM1_data/'+REG_ID+'_20450101_20490101_mean_'
case_end1 = ".nc"
fstring1 ='%s%s%s' 
pref1=REG_ID+'_UKESM1_50Tgpy_Cpld'

case_start2 = '~/NetCDF_Files/UKESM1_data/CTL_20450101_20490101_mean_'
case_end2 = ".nc"
pref2='UKESM1_control'
fstring2 ='%s%s%s' 

if True:
    # fixed SST simulations
    case_start1 = '~/NetCDF_Files/UKESM1_data/'+REG_ID+'_AtmosOnly_19840101_19880101_mean_'
    case_end1 = ".nc"
    fstring1 ='%s%s%s' 
    pref1=REG_ID+'_50Tgpy_FixedSST'

    case_start2 = '~/NetCDF_Files/UKESM1_data/CTL_AtmosOnly_19840101_19880101_mean_'
    case_end2 = ".nc"
    pref2='Control'
    fstring2 ='%s%s%s' 

Varname='<Varname>'
ind1 = fstring1 % (case_start1,Varname,case_end1)
print('example string used for file open',ind1)


for Varname in Varlist:
    print()
    print('-------------------------------'+Varname)    
    ind1 = fstring1 % (case_start1,Varname,case_end1)
    #print('xxx opening',ind1)
    DS1 = xr.open_mfdataset(ind1)
    DS1 = fix_UKMO_ds(ind1, DS1)
    #print('DS1.lon',DS1.lon.values)
    #DS1 = center_time(DS1)
    VN = Vdict[Varname]
    print('VN is ',VN)
    V1 = xr_getvar_sl(VN,DS1,method='maxb850')
    #print('V1',V1)
    #V1 = Var1.mean(dim='time',keep_attrs=True)
       
    ind2 = fstring2 % (case_start2,Varname,case_end2)
    #print('opening ind2',ind2)
    #DS2 = xr.open_mfdataset(ind2)
    DS2 = xr.open_mfdataset(ind2)
    DS2 = fix_UKMO_ds(ind1, DS2)
    V2 = xr_getvar_sl(VN,DS2,method='maxb850')


    DV = V1-V2
    print('DV range', DV.min().values, DV.max().values)
    if 'area' in DS1:
        area = DS1['area']
    elif 'area' in DS2:
        area = DS2['area']
    else:
        print('calculating weights')
        lat = V1['lat'].values
        lon = V1['lon'].values
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
    sDVA = ' (%5.2f)' % DVA
    print('area avgs '+pref1+' %5.2f' % (V1A.values),' '+pref2+' %5.2f' % (V2A.values),' Delta %5.2f' % (DVA.values))
    
    cbartitle = None
    if VN == 'unknown':
        cbartitle = Varname

    filefmt = 'pdf'
    
    if V1.min().values == V1.max().values:
        print('constant field skipping plot ')
    else:
        clev_rng = {'CDNUMC':np.array([0.,3.e11]),'FSNT':np.array([40.,360]),
                    'TGCLDLWP':np.array([0.,280.]),'PRECL':np.array([0.,10.]),
                    'PRECT':np.array([0.,16.]),'SWCF':np.array([-140.,0.]),
                    'CLDLOW':np.array([0.,90.]),'XXX':np.array([-45.,45.]),
                   }
        dlev_rng = {'CDNUMC':np.array([0.,3.e11])/2.,'FSNT':np.array([-45.,45.]),
                   'TGCLDLWP':np.array([-80.,80.]),'PRECL':np.array([-1.,1.]),
                    'PRECT':np.array([-5.,5.]),'SWCF':np.array([-45.,45.]),
                    'CLDLOW':np.array([-10.,10.]),'XXX':np.array([-45.,45.]),
                   }
        if V1.name in clev_rng:
            clevs = findNiceContours(clev_rng[V1.name],nlevs = 12)
        else:
            clevs = findNiceContours(np.array([V1.values,V2.values]),nlevs = 12)
        if V1.name in dlev_rng:
            dlevs = findNiceContours(dlev_rng[V1.name],nlevs = 15,rmClev=0.,sym=True)
        else:
            dlevs = findNiceContours(np.array([DV.min().values,DV.max().values]),nlevs = 15, rmClev=0.,sym=True)
        #dlevs = [-5.,-2.,-1.,-0.5,-0.2,-0.1,0.1,0.2,0.5,1.,2.,5.]
        #print('xxx',dlevs)
        dmap = diverge_map()

        plconf = '3-1x1'
        #plconf = '1x3'
        # good setup for 1 row of 3 columns
        if plconf == '1x3':
            fig, axes = plt.subplots(ncols=3
                                     ,gridspec_kw={'width_ratios': [1, 1, 1]}
                                     ,subplot_kw={'projection': ccrs.Mollweide()}
                                     ,figsize=(16,5)
                                    )

            xr_llhplot(V1, ax=axes[0],clevs=clevs,title=pref1+sV1A, cbartitle=cbartitle)
            xr_llhplot(V2, ax=axes[1],clevs=clevs,ylabels=False,title=pref2+sV2A, cbartitle=cbartitle)
            xr_llhplot(DV, ax=axes[2],clevs=dlevs,cmap=dmap,title=pref1+'-'+pref2+sDVA, cbartitle=cbartitle)
            plt.savefig(pref1+'_'+Varname+'.'+filefmt,format=filefmt,dpi=300)
            plt.show()
            
        # good setup for 3 rows of 1 columns
        if plconf == '3-1x1':
            
            fig, axes = setfig3b1x1()
            xr_llhplot(V1, ax=axes,clevs=clevs,title=pref1+sV1A, cbartitle=cbartitle)
            plt.savefig(pref1+'_'+Varname+'.'+filefmt,format=filefmt,dpi=300)
            plt.show()
            
            fig, axes = setfig3b1x1()
            xr_llhplot(V2, ax=axes,clevs=clevs,ylabels=False,title=pref2+sV2A, cbartitle=cbartitle)
            plt.savefig(pref2+'_'+Varname+'.'+filefmt,format=filefmt,dpi=300)
            plt.show()
            
            fig, axes = setfig3b1x1()
            xr_llhplot(DV, ax=axes,clevs=dlevs,cmap=dmap,title=pref1+'-'+pref2+sDVA, cbartitle=cbartitle)
            #plt.savefig(pref1+'_'+Varname+'-D.jpg',format='jpg',dpi=300)
            plt.savefig(pref1+'_'+Varname+'-D.'+filefmt,format=filefmt,dpi=300)
            pltllbox([-150.,-110.],[0.,30.])
            pltllbox([-110.,-70.],[-30.,0.])
            pltllbox([-25.,15.],[-30.,0.])
            plt.show()


        
    print('field processing complete')

```

```python
def xr_getvar_UK(case_start1,Varname,case_end1):
    """ get a field from a netCDF file.
    If it is a multi-level field, do something useful to provide single level information
    This function assumes that the x-xoordinate increases monotonically

    """
    
    Vdict = {'net_ToA_LW_W_m2':'toa_outgoing_longwave_flux'
        ,'net_ToA_SW_W_m2':'FSNT'
        ,'AOD_550nm':'AOD'
        ,'LWP_kg_m2':'TGCLDLWP'
        ,'p_surface_Pa':'PS'
        ,'T_surface_K':'TS'
        ,'precip_rate_kg_m2_sec':'PRECT'
        ,'PBL_depth_metres':'PBLH'
        ,'cloudtop_r_e_microns':'REPJR'
        ,'cloud_fraction':'CFPJR'
        ,'Outgoing_SW_Clear_W_m2':'MFSUTC'
        }

    ind1 = fstring1 % (case_start1,Varname,case_end1)
    #print('xxx opening',ind1)
    DS1 = xr.open_mfdataset(ind1)
    DS1 = fix_UKMO_ds(ind1, DS1)

    VN = Vdict[Varname]
    # print('VN is ',VN)
    Var1 = xr_getvar_sl(VN,DS1,method='maxb850')
    V1 = Var1
    return V1

Varlist = np.array(['LWP_kg_m2'])
#Varlist = np.array(['p_surface_Pa'])
#Varlist = np.array(['T_surface_K'])
#Varlist = np.array(['AOD_550nm','LWP_kg_m2','net_ToA_LW_W_m2','net_ToA_SW_W_m2'])
#Varlist = np.array(['Outgoing_SW_Clear_W_m2','p_surface_Pa','T_surface_K','precip_rate_kg_m2_sec','PBL_depth_metres','cloudtop_r_e_microns','AOD_550nm','LWP_kg_m2','net_ToA_LW_W_m2','net_ToA_SW_W_m2'])
Varlist = np.array(['Outgoing_SW_Clear_W_m2','precip_rate_kg_m2_sec','PBL_depth_metres','cloudtop_r_e_microns','AOD_550nm','LWP_kg_m2','net_ToA_LW_W_m2','net_ToA_SW_W_m2'])
Varlist = np.array(['cloud_fraction'])
Varlist = np.array(['Outgoing_SW_Clear_W_m2','net_ToA_SW_W_m2'])

REG_ID = 'R1_NEP'
REG_ID = 'R2_SEP'
REG_ID = 'R3_SEA'

# coupled simulations
case_start1 = '~/NetCDF_Files/UKESM1_data/'+REG_ID+'_20450101_20490101_mean_'
case_end1 = ".nc"
fstring1 ='%s%s%s' 
pref1=REG_ID+'_UKESM1_50Tgpy_Cpld'

case_start2 = '~/NetCDF_Files/UKESM1_data/CTL_20450101_20490101_mean_'
case_end2 = ".nc"
pref2='UKESM1_control'
fstring2 ='%s%s%s' 

if True:
    # fixed SST simulations
    case_start1 = '~/NetCDF_Files/UKESM1_data/'+REG_ID+'_AtmosOnly_19840101_19880101_mean_'
    case_end1 = ".nc"
    fstring1 ='%s%s%s' 
    pref1=REG_ID+'_50Tgpy_FixedSST'

    case_start2 = '~/NetCDF_Files/UKESM1_data/CTL_AtmosOnly_19840101_19880101_mean_'
    case_end2 = ".nc"
    pref2='Control'
    fstring2 ='%s%s%s' 

Varname = 'Outgoing_SW_Clear_W_m2'
MFSUTC1 = xr_getvar_UK(case_start1,Varname,case_end1)
MFSUTC2 = xr_getvar_UK(case_start2,Varname,case_end2)

Varname = 'net_ToA_SW_W_m2'
FSNT1 = xr_getvar_UK(case_start1,Varname,case_end1)
FSNT2 = xr_getvar_UK(case_start2,Varname,case_end2)

SWCF1 = FSNT1+MFSUTC1
SWCF2 = FSNT2+MFSUTC2
print('xxx',SWCF1)

V1 = FSNT1
V2 = MFSUTC1
DV = SWCF1
clevs = None;
dlevs = None;

print('DV range', DV.min().values, DV.max().values)
if 'area' in DS1:
    area = DS1['area']
elif 'area' in DS2:
    area = DS2['area']
else:
    print('calculating weights')
    lat = V1['lat'].values
    lon = V1['lon'].values
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
sDVA = ' (%5.2f)' % DVA
print('area avgs '+pref1+' %5.2f' % (V1A.values),' '+pref2+' %5.2f' % (V2A.values),' Delta %5.2f' % (DVA.values))

cbartitle = None
if VN == 'unknown':
    cbartitle = Varname

filefmt = 'pdf'


if plconf == '3-1x1':

    fig, axes = setfig3b1x1()
    xr_llhplot(V1, ax=axes,clevs=clevs,title=pref1+sV1A, cbartitle=cbartitle)
    plt.savefig(pref1+'_'+Varname+'.'+filefmt,format=filefmt,dpi=300)
    plt.show()

    fig, axes = setfig3b1x1()
    xr_llhplot(V2, ax=axes,clevs=clevs,ylabels=False,title=pref2+sV2A, cbartitle=cbartitle)
    plt.savefig(pref2+'_'+Varname+'.'+filefmt,format=filefmt,dpi=300)
    plt.show()

    fig, axes = setfig3b1x1()
    xr_llhplot(DV, ax=axes,clevs=dlevs,cmap=dmap,title=pref1+'-'+pref2+sDVA, cbartitle=cbartitle)
    #plt.savefig(pref1+'_'+Varname+'-D.jpg',format='jpg',dpi=300)
    plt.savefig(pref1+'_'+Varname+'-D.'+filefmt,format=filefmt,dpi=300)
    pltllbox([-150.,-110.],[0.,30.])
    pltllbox([-110.,-70.],[-30.,0.])
    pltllbox([-25.,15.],[-30.,0.])
    plt.show()

print('field processing complete')

```

```python
1./0.
```

```python
print('level_height',level_height.values)
print('sigma',sigma.values)
print('surface_altitude',surface_altitude)
print('alt',altitude)
i = 1
j = 1
#print('sij',surface_altitude[i,j].values)
#print('aij',altitude[:,i,j].values)
dz1 = altitude - surface_altitude
dz2 = (sigma - 1)*surface_altitude + level_height
#print('dz1',dz1[:,i,j].values)
#print('dz2',dz2[:,i,j].values)
pmb = 850.
psmb = 1000.
scaleheight = 8.4e3
altmb = -np.log(pmb/psmb)*scaleheight
print('alt850, i,j,pmb, alt near pmb',i,j,pmb, altmb) # approximate height of pressure surface at pmb hPa


# identify the vertical level where dz1 is nearest altmb
llmin1 =  np.abs(dz1.values-altmb).argmin(axis=0)
print('llmin',np.shape(llmin1))
lind = llmin1.flatten()
print('lind',np.shape(lind),lind[0:3])
indi = np.arange(0,len(lind))
print('indi',np.shape(indi))
ni,nj,nk = np.shape(dz1.values)
print('ni,nj,nk',ni,nj,nk)
VAR = dz1.copy()
# now extract data for each column
data = VAR.values
datas = data.reshape((ni,nj*nk))
#print('datas',np.shape(datas))
# select the level nearest the index
dataz = datas[lind,indi]
#print('dataz',np.shape(dataz))
datazr = dataz.reshape((nj,nk))
#print('datazr shape ij',datazr.shape, datazr[i,j])
#datam = data[llmin1,:,:]
#print('datam',datam.shape)
VARs = VAR.sel(model_level_number=1,method='nearest')
VARs.attrs['long_name'] = VARs.attrs['long_name']+' nearest altmb'
VARs.data = datazr

print('Vars.data',i,j,VARs[i,j].values)
```

```python
print('level_height',level_height.values)
print('sigma',sigma.values)
print('surface_altitude',surface_altitude)
print('alt',altitude)
i = 100
j = 100
#print('sij',surface_altitude[i,j].values)
#print('aij',altitude[:,i,j].values)
dz1 = altitude - surface_altitude
dz2 = (sigma - 1)*surface_altitude + level_height
#print('dz1',dz1[:,i,j].values)
#print('dz2',dz2[:,i,j].values)
pmb = 850.
psmb = 1000.
scaleheight = 8.4e3
altmb = -np.log(pmb/psmb)*scaleheight
print('alt850, i,j,pmb, alt near pmb',i,j,pmb, altmb) # approximate height of pressure surface at pmb hPa

# mask values where dz1 > altmb
VAR = dz1.copy()
V2 = VAR.where(dz1 <= altmb+50.)
print('V2',V2)
print('V2vals',V2[:,i,j].values)
V3 = V2.max(dim='model_level_number')
V3.attrs['long_name'] = 'max value below pmb of '+V2.attrs['long_name']
print('V3',V3)
print('V3.values',V3[i,j].values)
if False:
    # identify the vertical level where dz1 is nearest altmb
    llmin1 =  np.abs(dz1.values-altmb).argmin(axis=0)
    print('llmin',np.shape(llmin1))
    lind = llmin1.flatten()
    print('lind',np.shape(lind),lind[0:3])
    indi = np.arange(0,len(lind))
    print('indi',np.shape(indi))
    ni,nj,nk = np.shape(dz1.values)
    print('ni,nj,nk',ni,nj,nk)
    VAR = dz1.copy()
    # now extract data for each column
    data = VAR.values
    datas = data.reshape((ni,nj*nk))
    #print('datas',np.shape(datas))
    # select the level nearest the index
    dataz = datas[lind,indi]
    #print('dataz',np.shape(dataz))
    datazr = dataz.reshape((nj,nk))
    #print('datazr shape ij',datazr.shape, datazr[i,j])
    #datam = data[llmin1,:,:]
    #print('datam',datam.shape)
    VARs = VAR.sel(model_level_number=1,method='nearest')
    VARs.attrs['long_name'] = VARs.attrs['long_name']+' nearest altmb'
    VARs.data = datazr

    print('Vars.data',i,j,VARs[i,j].values)
```

$$
   alt(i,j) = levelheight(l) + (\sigma(l) * surf(i,j))
$$

$$
\begin{aligned}
  dz(i,j) &= alt(i,j) - surf(i,j) \cr
  &= (levelheight(l) + (\sigma(l)*surf(i,j)) - surf(i,j) \cr
  &= (levelheight(l) + (\sigma(l)-1)*surf(i,j) \cr
\end{aligned}
$$


```python
# sample multilevel field on hybrid levels
fig, axes = setfig3b1x1()

def getvarDSM(Varname,fstring,case_start,case_end):
    """getvar DSM
       get variable from file specifying the formatting
    """
    ind = fstring % (case_start,Varname,case_end)
    #print('opening',ind1)
    DS = xr.open_mfdataset(ind)
    Var = xr_getvar(Varname,DS)
    VM = Var.mean(dim='time',keep_attrs=True)
    return VM;

mylev=850.
FREQL = getvarDSM('FREQL', fstring1, case_start1, case_end1).sel(lev=mylev,method='nearest')
CLOUD = getvarDSM('CLOUD', fstring1, case_start1, case_end1).sel(lev=mylev,method='nearest')/100.

NUMLIQ = getvarDSM('NUMLIQ', fstring1, case_start1, case_end1).sel(lev=mylev,method='nearest')*1.e-6
ICNUMLIQ = NUMLIQ/(CLOUD+1.e-2)
ICNUMLIQ = ICNUMLIQ.rename('ICNUMLIQ')
ICNUMLIQ.attrs['long_name'] = 'approx in-cl number @850hPa hybsfc'
ICNUMLIQ.attrs['units'] = '#/cc'

xr_llhplot(ICNUMLIQ, ax=axes)#,clevs=clevs,title=pref1+sV1A)
#plt.savefig(pref1+'_'+Varname+'.jpg',format='jpg',dpi=300)
plt.show()

```

```python
# sample multilevel field at level nearest PBLH
fig, axes = setfig3b1x1()

PBLH = getvarDSM('PBLH', fstring1, case_start1, case_end1)
Z3 = getvarDSM('Z3', fstring1, case_start1, case_end1)
#P3 = getvarDSM('P3', fstring1, case_start1, case_end1)
PHIS = getvarDSM('PHIS', fstring1, case_start1, case_end1)

# find the 1D array of indices closest to the PBLH
# I think Z3 is height above sea-level, and PBLH is height above surface
Z3D = Z3 - PHIS/9.8 - PBLH
llmin1 =  np.abs(Z3D.values).argmin(axis=0)
#print('llmin',np.shape(llmin1))
lind = llmin1.flatten()
#print('lind',np.shape(lind),lind[0:3])
indi = np.arange(0,len(lind))
#print('indi',np.shape(indi))
ni,nj,nk = np.shape(Z3.values)
#print('ni,nj,nk',ni,nj,nk)

NUMLIQ = getvarDSM('NUMLIQ', fstring1, case_start1, case_end1)*1.e-6
CLOUD = getvarDSM('CLOUD', fstring1, case_start1, case_end1)/100.
ICNUMLIQ = NUMLIQ/(CLOUD+1.e-2)
ICNUMLIQ = ICNUMLIQ.rename('ICNUMLIQ')
ICNUMLIQ.attrs['long_name'] = 'approx in-cloud number conc'
ICNUMLIQ.attrs['units'] = '#/cc'

VAR = ICNUMLIQ.copy()
# now extract data for each column
data = VAR.values
datas = data.reshape((ni,nj*nk))
#print('datas',np.shape(datas))
dataz = datas[lind,indi]
#print('dataz',np.shape(dataz))
datazr = dataz.reshape((nj,nk))
#print('datazr shape ij',datazr.shape, datazr[i,j])
#datam = data[llmin1,:,:]
#print('datam',datam.shape)
VARs = VAR.sel(lev=1000.,method='nearest')
VARs.attrs['long_name'] = VARs.attrs['long_name']+' near PBLH'
VARs.data = datazr

xr_llhplot(VARs, ax=axes)#,clevs=clevs,title=pref1+sV1A)
#plt.savefig(pref1+'_'+Varname+'.jpg',format='jpg',dpi=300)
plt.show()

```

```python

lon = xr_getvar('lon',DS1)
lat = xr_getvar('lat',DS1)
DVcum = 0.
#DV = V1
print('analyzing variable DV.name=',DV.name)
for i in range(len(xreg)):
    xri = xreg[i]
    yri = yreg[i]
    print('namereg, xreg',namereg[i],xri)
    #print('lon',lon.values)
    #R2 = DV.sel(lon=slice(220,250),lat=slice(15,35))
    #%print('R2',R2)
    # allow longitudes to cross greenwich
    print('DV range', DV.min().values, DV.max().values)

    if xri[1]>xri[0]: 
        DVI = DV.sel(lon=(lon < xri[1]) & (lon >= xri[0]),
                     lat=(lat >=yri[0]) & (lat < yri[1])
                    ).load()
        areai = area.sel(lon=(lon < xri[1]) & (lon >= xri[0]),
                     lat=(lat >=yri[0]) & (lat < yri[1])
                    ).load()
    else:
        DVI = DV.sel(lon=(lon >= xri[0]) | (lon < xri[1]),
                     lat=(lat >=yri[0]) & (lat < yri[1])
                    ).load()
        areai = area.sel(lon=(lon >= xri[0]) | (lon < xri[1]),
                     lat=(lat >=yri[0]) & (lat < yri[1])
                    ).load()

    #print('DVI',DVI)
    #print('DVI range', DVI.min().values, DVI.max().values)

    DVIG = (DVI*areai).sum()/area.sum()
    DVcum = DVcum + DVIG.values
    print('i, areai, area, fglob, focean ', i, areai.sum().values, area.sum().values,
         (areai*of).sum().values/area.sum().values,
         (areai*of).sum().values/(area*of).sum().values)
    #print('i, DVI.sum, DVI*areai,', DVI.sum().values, (DVI*areai).sum().values)
    print('i, DVIG and cum',i, DVIG.values, DVcum)
    print(' ')
    #break
```

```python
xmesh, ymesh = np.meshgrid(lon, lat)
indlf = fstring1 % (case_start1,'LANDFRAC',case_end1)
DSLF = xr.open_mfdataset(indlf)
lf = xr_getvar('LANDFRAC',DSLF).squeeze()
indof = fstring1 % (case_start1,'OCNFRAC',case_end1)
DSOF = xr.open_mfdataset(indof)
of = xr_getvar('OCNFRAC',DSOF).squeeze()
indif = fstring1 % (case_start1,'ICEFRAC',case_end1)
DSIF = xr.open_mfdataset(indif)
ifr = xr_getvar('ICEFRAC',DSIF).squeeze()

DVcum = 0.
i = 2
if True:
#for i in range(len(xreg)):
    xri = xreg[i].copy()
    yri = yreg[i]
    # find indices in region of interest
    if xri[1]>=xri[0]:
        inds = ((xmesh>xri[0]) & (xmesh <= xri[1]) & (ymesh > yri[0]) & (ymesh <= yri[1]))
    else:
        inds = (((xmesh >= xri[0]) | (xmesh < xri[1])) &
                (ymesh >=yri[0]) & (ymesh < yri[1]))
    # remove points with land from the mask
    inds = inds & (lf < 0.02)
    #print('DISABLED land mask')
    #DVx = DV.where((xmesh>xri[0]) & (xmesh <= xri[1]) & (ymesh > yri[0]) & (ymesh <= yri[1]))
    DVx = DV.where(inds)
    areax = area.where(inds)
    DVIG = (DVx*areax).sum()/area.sum()
    DVIL = (DVx*areax).sum()/areax.sum()
    DVcum = DVcum + DVIG.values
    fglob = (areax.sum()/area.sum()).values
    focn = (areax.sum()/(area*of).sum()).values
    print('areax',areax.sum().values)
    print('area',area.sum().values)
    print('areao',(area*of).sum().values)
    print('areai',(area*ifr).sum().values)
    print('areal',(area*lf).sum().values)
    print(i," DVIG = {:.2f}".format(DVIG.values), 
          "DVcum = {:.2f}".format(DVcum),
          "DVIL = {:.2f}".format(DVIL.values), 
          "fglob = {:.2%}".format(fglob),
          "focn = {:.2%}".format(focn)
         )
    if xri[1] < xri[0]:
        xri[1] += 360.
    regcx = [xri[0],xri[1],xri[1],xri[0],xri[0]]
    regcy = [yri[0],yri[0],yri[1],yri[1],yri[0]]
    #fig, axes = setfig3b1x1()
    map_proj = ccrs.Mollweide()
    map_proj._threshold /= 100.
    fig, axes = plt.subplots(ncols=1,
                             gridspec_kw={'width_ratios': [1]},
                             subplot_kw={'projection':map_proj},
                             figsize=(6,3),
                            )
    fig.set_dpi(300.0)
    xr_llhplot(DVx, ax=axes)
    plt.plot(regcx,regcy,color='green',transform=ccrs.PlateCarree())
    plt.show()
    #break
```

```python
xmesh, ymesh = np.meshgrid(lon, lat)
indlf = fstring1 % (case_start1,'LANDFRAC',case_end1)
DSLF = xr.open_mfdataset(indlf)

# add new lon coordinate that runs from -180 to 180. for use when region crosses the gm
lf = xr_getvar('LANDFRAC',DSLF).squeeze()
lf = lf.assign_coords(longm=(((lf.lon +180.) % 360) - 180))

area = xr_getvar('area',DSLF).squeeze()
area = area.assign_coords(longm=(((area.lon + 180) % 360) - 180))

DV = DV.assign_coords(longm=(((DV.lon + 180) % 360) - 180))
#print('DV',DV)

DVcum = 0.
i = 2
if True:
#for i in range(len(xreg)):
    print(' ')
    print('xreg i ', xreg[i])
    xri = xreg[i].copy()
    xrigm = ((xreg[i]+180)%360)-180.
    print('xrigm',xrigm)
    yri = yreg[i]
    # find indices in region of interest
    indlat = np.where((DV.lat.values>yri[0]) & (DV.lat.values <= yri[1]))[0]
    print('indlat',indlat)
    if xri[1]>=xri[0]:
        indlon = np.where((DV.lon.values>xri[0]) & (DV.lon.values<= xri[1]))[0]
        DVxx = DV.isel(lon=indlon,lat=indlat).copy()
        print('DVxx',DVxx)
        areaxx = area.isel(lon=indlon,lat=indlat)
        lfxx = lf.isel(lon=indlon,lat=indlat)
    else:
        indlon = np.where((DV.longm.values>xrigm[0]) & (DV.longm.values<= xrigm[1]))[0]
        DVxx = DV.isel(lon=indlon,lat=indlat).copy()
        print('DVxx',DVxx)
        areaxx = area.isel(lon=indlon,lat=indlat)
        lfxx = lf.isel(lon=indlon,lat=indlat)
    print('indlon',indlon)
    print('DVxx.longm',DVxx.longm.values)
    print('DVnew', DVxx.sortby("longm"))

    """
    # remove points with land from the mask
    inds = inds & (lf < 0.02)
    #print('DISABLED land mask')
    #DVx = DV.where((xmesh>xri[0]) & (xmesh <= xri[1]) & (ymesh > yri[0]) & (ymesh <= yri[1]))
    DVx = DV.where(inds)
    areax = area.where(inds)
    DVIG = (DVx*areax).sum()/area.sum()
    DVIL = (DVx*areax).sum()/areax.sum()
    DVcum = DVcum + DVIG.values
    fglob = (areax.sum()/area.sum()).values
    focn = (areax.sum()/(area*of).sum()).values
    print('areax',areax.sum().values)
    print('area',area.sum().values)
    print('areao',(area*of).sum().values)
    print('areai',(area*ifr).sum().values)
    print('areal',(area*lf).sum().values)
    print(i," DVIG = {:.2f}".format(DVIG.values), 
          "DVcum = {:.2f}".format(DVcum),
          "DVIL = {:.2f}".format(DVIL.values), 
          "fglob = {:.2%}".format(fglob),
          "focn = {:.2%}".format(focn)
         )
    if xri[1] < xri[0]:
        xri[1] += 360.
    regcx = [xri[0],xri[1],xri[1],xri[0],xri[0]]
    regcy = [yri[0],yri[0],yri[1],yri[1],yri[0]]
    #fig, axes = setfig3b1x1()
    map_proj = ccrs.Mollweide()
    map_proj._threshold /= 100.
    print('xri',xri)
    print('2nd xreg i', xreg[i])
    fig, axes = plt.subplots(ncols=1,
                             gridspec_kw={'width_ratios': [1]},
                             subplot_kw={'projection':map_proj},
                             figsize=(6,3),
                            )
    fig.set_dpi(300.0)
    xr_llhplot(DVx, ax=axes)
    plt.plot(regcx,regcy,color='green',transform=ccrs.PlateCarree())
    plt.show()
    #break
    """
```

```python

if True:
    map_proj = ccrs.Mollweide()
    #map_proj = ccrs.PlateCarree()

    map_proj._threshold /= 100.
    fig, axes = plt.subplots(ncols=1,
                             gridspec_kw={'width_ratios': [1]},
                             subplot_kw={'projection':map_proj},
                             figsize=(6,3),
                            )
    fig.set_dpi(300.0)
    #xr_llhplot(DVxx, ax=axes)
    DVnew = DVxx.sortby("longm")
    #DVnew = DVxx.sortby("lon")
    xmesh, ymesh = np.meshgrid(DVnew.longm.values, DVnew.lat.values)
    axes.coastlines()
    axes.contourf(xmesh, ymesh, DVnew.values, transform=ccrs.PlateCarree())
    ylabels = True 
    gl = axes.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,
                      linewidth=2, color='gray', alpha=0.5)
    gl.left_labels=ylabels
    gl.right_labels=ylabels
    gl.bottom_labels= True
    #axes.set_extent([-180., 180., -90., 90.], crs=ccrs.PlateCarree())
    #axes.set_extent([-180., 180., -90., 90.], crs=map_proj)
    axes.set_global()
    plt.show()
```

```python
lonsub = DS1['lon'+regtag]#.where(pmask)#.isel(time=0)
latsub = DS1['lat'+regtag]#.where(pmask)#.isel(time=0)
lat = latsub.values
print('subreg',lonsub.min().values,lonsub.max().values, latsub.min().values,latsub.max().values)
#print('subreg size',lonsub.shape)
print('shape and size of variables',lonsub.shape, lonsub.size,' number of unmasked cells ',np.count_nonzero(lonsub.notnull().values))

if 'area' in DS1:
    print('area present in DS1')
    area = xr_getvar('area',DS1,regtag=regtag)#.where(pmask)
    wtsh = area.fillna(0).values
else:
    print('area not found')
    wtsh = make_fvarea(lonsub.values,latsub.values)

#print('weights shape',wtsh.shape,wtsh.sum())

Varlist = np.array(['T','AREL','CLOUD'])
#Varlist = np.array(['AREL'])

#print('wtsh',wtsh)

for Vname in Varlist:
    print()
    print('-------------------------------')
    ind1 = fstring1 % (case_start1,Vname,case_end1)
    #print('ind1 is ', ind1)
    DS1 = xr.open_mfdataset(ind1).transpose(...,'lev') # make sure lev is last index
    DS1 = center_time(DS1)
    ind2 = fstring2 % (case_start2,Vname,case_end2)  # make sure lev is last index
    DS2 = xr.open_mfdataset(ind2).transpose(...,'lev') 
    DS2 = center_time(DS2)
    
    V1 = xr_getvar(Vname,DS1,regtag=regtag)#.where(pmask)
    if V1.min().values == V1.max().values:
        print('constant field skipping plot ')
    else:
        V1 = xr_getvar(Vname, DS1, regtag).mean(dim='time',keep_attrs=True)#.where(pmask).squeeze()
        #print('V1xxx ', V1)
        #print('V1',V1.min().values,V1.max().values)
        #print('V1 shape, size, realsize', V1.shape, np.size(V1.values), np.count_nonzero(V1.notnull().values) )
        V2 = xr_getvar(Vname, DS2, regtag).mean(dim='time',keep_attrs=True)#.where(pmask).squeeze()
        DV = V1-V2
        print(Vname, V1.attrs['long_name'],'Range V1 and V2 ',V1.min().values, V1.max().values, V2.min().values, V2.max().values)

        # """

        DPOG1 = xr_getvar('DPOG',DS1,regtag=regtag).mean(dim='time',keep_attrs=True)
        weights1 = wtsh[...,np.newaxis]*DPOG1
        DPOG2 = xr_getvar('DPOG',DS2,regtag=regtag).mean(dim='time',keep_attrs=True)
        weights2 = wtsh[...,np.newaxis]*DPOG2
        V1A = V1.weighted(weights1).mean()
        V2A = V2.weighted(weights2).mean()
        print('Area/Mass weighted avg '+pref1+' %5.3f' % (V1A.values),' '+pref2+' %5.3f' % (V2A.values))
        
        # """
        
        data1 = V1.mean(dim='lon',keep_attrs=True).transpose().values
        #print('data1 shape',data1.shape)
        data2 = V2.mean(dim='lon',keep_attrs=True).transpose().values
        datad = data1-data2
        

        lev = V1['lev'].values
        #print('lev shape',lev.shape)
    
#        plotZMf(data1, lat_h, lev)
        fig, axes = plt.subplots(ncols=3
                                 ,gridspec_kw={'width_ratios': [1, 1, 1]}
#                                 ,subplot_kw={'projection': ccrs.PlateCarree()}
                                 ,figsize=(16,5)
                                )
        ytop = 20.
        plotZMf(data1, lat, lev,axesa=axes[0],plotOpt={'colorbar':"botnd",'units':V1.units,'ltitle':pref1,'ytop':ytop})
        plotZMf(data2, lat, lev,axesa=axes[1],plotOpt={'colorbar':"bot",'units':V2.units,'ltitle':pref2,'rtitle':V2.long_name,'ytop':ytop})
        dlevs = findNiceContours(np.array([datad.min(),datad.max()]),nlevs = 10, rmClev=0.,sym=True)
        dmap = diverge_map()
        plotZMf(datad, lat, lev,axesa=axes[2],plotOpt={'clevs':dlevs,'cmap':dmap,'colorbar':"bot",'units':V2.units,'ytop':ytop,'ltitle':pref1+'-'+pref2})

        #print('attribute check on xarray',hasattr(V1,'units'))
        #plt.savefig('test_'+Vname+'.jpg',format='jpg')
        #plt.tight_layout()
        plt.show()
      

```

```python
# next cells not used
help(xr_getvar)
1./0.
```

```python
Varlist = np.array(['T','Q','CLOUD','CLDLIQ','ICWMR','CLDICE','RELHUM','NUMICE','NUMLIQ','Mass_bc'])
#                    RESTOM','FLNT','FSNT','TS','TMQ','PRECT','AEROD_v','CLDLOW','CLDTOT','LWCF','SWCF','TGCLDIWP','TGCLDLWP','SHFLX','LHFLX','PBLH','PCONVT','PRECC','PRECS'])
#Varlist = np.array(['TS','TMQ','PRECT'])
#Varlist = np.array(['RESTOM','LWCF','SWCF','FLNT','FSNT'])

Varname = 'TS' # Varlist[0]

case_start1 = "/home/jupyter-haruki/work/CESM_MCB/Fixed_SST/"
case_end1 = "Fixed_SST.cam.h0."+Varname+".y1-19.nc"
ind1 = case_start1+case_end1
DS1 = xr.open_mfdataset(ind1)


Var1a = DS1[Varname]
Var1am = Var1a.mean(dim='time')
DS1 = center_time(DS1)
Var1 = DS1[Varname]
#print('xxx',Var1)

case_start2 = "/home/jupyter-haruki/work/CESM_MCB/MCB_R1R2R3_CN375cm/" 
case_end2 = "MCB_R1R2R3_CN375cm.cam.h0."+Varname+".y1-10.nc"
ind2 = case_start2+case_end2
DS2 = xr.open_mfdataset(ind2)
DS2 = center_time(DS2)
Var2 = DS2[Varname]
#print('yyy',Var2)

```
