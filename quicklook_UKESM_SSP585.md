---
jupyter:
  jupytext:
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.15.2
  kernelspec:
    display_name: pjrpy3
    language: python
    name: pjrpy3
---

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

```

```python
Varlist = np.array(['ts'])
# coupled simulations
case_start1 = '~/NetCDF_Files/UKESM1_data_v2/Coupled_SSP585/'
case_end1 = '_Amon_UKESM1-0-LL_ssp585_r1i1p1f2_gn_205001-210012.nc'
case_end1 = '_Amon_UKESM1-0-LL_ssp585_r1i1p1f2_gn_201501-204912.nc'
case_end1 = '_Amon_UKESM1-0-LL_ssp585_r1i1p1f2_gn_*.nc'
case_start1 = '~/NetCDF_Files/UKESM1_data_v2/Coupled_SSP245/'
case_end1 = '_Amon_UKESM1-0-LL_ssp245_r1i1p1f2_gn_*.nc'
fstring1 ='%s%s%s' 
REG_ID=''
pref1=REG_ID+'_UKESM1_50Tgpy_Cpld'
Varname=Varlist[0]
#print('xxx',Varname)
ind1 = fstring1 % (case_start1,Varname,case_end1)
print('example string used for file open',ind1)
DS1 = xr.open_mfdataset(ind1)
print(DS1)
time = DS1['time']
time_bnds = DS1['time_bnds']
Var = DS1[Varname]
#print('time',time[0])
#print('time_bnds',time_bnds[0].values)
#print('Var',Var[0])
```

```python
def fix_UKMOds_ts(filename, dsin: xr.Dataset) -> xr.Dataset:
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
    
    # create an area variable if it isn't present
    if not 'area' in ds:
        #print('calculating weights')
        lat = ds['lat'].values
        lon = ds['lon'].values
        areavals = make_fvarea(lon,lat)
        #print('area',area.shape,area.sum()-4.*np.pi)
        area = xr.DataArray(areavals,coords=[lat, lon],dims=['lat','lon'])
        #area['units'] = 'steradian'
        #area['long_name'] = 'cell_area'
        ds['area'] = area
        ds['area'].attrs['units'] = 'steradian'
        ds['area'].attrs['long_name'] = 'areacella'
        
    for Vname in ds:
        #print('Vname',Vname)
        if (('PBL' in filename) & ('height_after' in Vname)):
            ds[Vname].attrs['UKESM name'] = Vname
            ds[Vname].attrs['units'] = 'kg/m2'
            ds[Vname].attrs['long_name'] = 'PBL height'
            ds = ds.rename({Vname:'PBLH'})
            #print('fixed PBLH')
            

    #print('ds',ds['area'].shape)
    return ds


```

```python
decade = 4
ybeg = 2000+10*decade+1
yend = 2000+10*(decade+1)+1
#yend = 2000+10*(decade+1)
fstring1 ='%d%s' 
tbeg = fstring1 % (ybeg,"-01-01")
tend = fstring1 % (yend,"-01-01")
DS2 = fix_UKMOds_ts(ind1, DS1)
DS2 = DS2.sel(time=slice(tbeg, tend))
DS2.time
```

```python
DS2 = fix_UKMOds_ts(ind1, DS1)
print(DS2)
aa = DS2['area']
print('aa',aa.sum().values,4*np.pi)
DS2 = DS2.sel(time=slice(tbeg, tend)).mean('time')
V1 = DS2['ts']
#print('V1',V1)
#DS2
weights = DS2.area
V1A = V1.weighted(weights).mean()
print('V1A',V1A.values)
if False:
    # write the decadal average to a netCDF file
    DS2.to_netcdf("saved_on_disk.nc")
```

```python

fig, axes = setfig3b1x1()
#print('V1XXX',V1)
clevs = findNiceContours(V1.values,nlevs = 12)
xr_llhplot(V1, ax=axes,clevs=clevs,title='title')
#pltllbox([-150.,-110.],[0.,30.])
#pltllbox([-110.,-70.],[-30.,0.])
#pltllbox([-25.,15.],[-30.,0.])
#plt.savefig(pref1+'_'+Varname+'.pdf',format='pdf',dpi=300)
plt.show()

```

```python

```
