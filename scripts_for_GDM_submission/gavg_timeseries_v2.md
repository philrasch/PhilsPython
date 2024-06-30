---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.16.2
  kernelspec:
    display_name: pjrpy3
    language: python
    name: pjrpy3
---

**plot timeseries of Surface Temperature**

```python
### script rewrite not completed because Haruki already built a nice plot update
1./0.
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
%run -i ~/Python/pjr3
```

```python
# open a file that will hold all the filenames used by the program
flname = '/tmp/flname'
file = open(flname, 'w')
file.write('list of files used by gavg_xxx\n')
file.close()
file = open(flname, 'a')
```

```python
fnprefix='GAVG_tsers_E3SM_and_CESM'
```

```python


var = 'TS'

experiment = "20230724.v2.LR.WCYCLSSP245.MCB-SSLT-EM.R1-3.test01"
s3_path = f's3://mcb-zarr/E3SMv2/{experiment}/atm/proc/tseries/month_1/{experiment}.eam.h0.{var}.zarr'
V3 = make_GA_ts2(var, s3_path )
print('V3',V3)

experiment = "20221014.v2.LR.WCYCLSSP245.E2_CNTL_01"
s3_path = f's3://mcb-zarr/E3SMv2/{experiment}/atm/proc/tseries/month_1/{experiment}.eam.h0.{var}.zarr'
#DS1 = s3load_zarr(s3_path)
#print('xxx',DS1)
#DS1 = xr.open_mfdataset(ind1)
#DS1 = center_time(DS1)
V2 = make_GA_ts2(var, s3_path )
print('V2',V2)

experiment = '20221128.v2.LR.WCYCLSSP245.E2_CNTL_01.R1-3_CDNC2000'
s3_path = f's3://mcb-zarr/E3SMv2/{experiment}/atm/proc/tseries/month_1/{experiment}.eam.h0.{var}.zarr'
V1 = make_GA_ts2(var, s3_path )
print('V1',V1)


```

```python

tit1 = "E3SM CDNC2000, NEP,SEP,SEA"
case1 = '20221128.v2.LR.WCYCLSSP245.E2_CNTL_01.R1-3_CDNC2000'
var = 'TS'
s3_path = f's3://mcb-zarr/E3SMv2/{case1}/atm/proc/tseries/month_1/{case1}.eam.h0.{var}.zarr'
print('yyy',s3_path)
DS1 = s3load_zarr(s3_path)
#print('DS1',DS1)

tit2 = "E3SM Control"
experiment = "20221014.v2.LR.WCYCLSSP245.E2_CNTL_01"
s3_path = f's3://mcb-zarr/E3SMv2/{experiment}/atm/proc/tseries/month_1/{experiment}.eam.h0.{var}.zarr'
DS2 = s3load_zarr(s3_path)
#print('DS2',DS2)

tit3 = "E3SM 49Tg/yr, NEP,SEP,SEA"
experiment = "20230724.v2.LR.WCYCLSSP245.MCB-SSLT-EM.R1-3.test01"
s3_path = f's3://mcb-zarr/E3SMv2/{experiment}/atm/proc/tseries/month_1/{experiment}.eam.h0.{var}.zarr'
DS3 = s3load_zarr(s3_path)
#print('DS3',DS3)

tit4 = "CESM CDNC600, NEP,SEP,SEA"
experiment = "b.e21.BSSP245smbb_MCB600cm_R1R2R3.f09_g17.LE2-1011.001"
s3_path = f's3://mcb-zarr/CESM2/{experiment}/atm/proc/tseries/month_1/{experiment}.cam.h0.{var}.2015-2064.zarr'
print('s3_path',s3_path)
DS4 = s3load_zarr(s3_path)
#print('DS4',DS4)

experiment = "b.e21.BSSP245smbb.f09_g17.001"
s3_path = f's3://mcb-zarr/CESM2/{experiment}/atm/proc/tseries/month_1/{experiment}.cam.h0.{var}.zarr'
print('s3_path',s3_path)
DS5 = s3load_zarr(s3_path)
print('DS5',DS5)
1./0.

ext5 = ".cam.h0"
ext5f = ".201501-206412.nc"
tit5 = "CESM Control(e1)"
string5 = "{dir:s}/{case:s}{ext:s}.{var:s}"+ext5f

dir6 = "/e3sm_prod/phil/timeseries/cesm2-mcb-reshaped/CESM2_SSP245/ens002"
case6 = "b.e21.BSSP245smbb.f09_g17.002"
ext6 = ".cam.h0"
ext6f = ".201501-206412.nc"
tit6 = "CESM Control(e2)"
string6 = "{dir:s}/{case:s}{ext:s}.{var:s}"+ext6f

dir7 = "/e3sm_prod/phil/timeseries/cesm2-mcb-reshaped/CESM2_SSP245/ens003"
case7 = "b.e21.BSSP245smbb.f09_g17.003"
ext7 = ".cam.h0"
ext7f = ".201501-206412.nc"
tit7 = "CESM Control(e3)"
string7 = "{dir:s}/{case:s}{ext:s}.{var:s}"+ext7f

dir8 = "/e3sm_prod/phil/timeseries/cesm2-mcb-reshaped"
case8 = "b.e21.BSSP245smbb_MCBss7TgYr_R1R2R3.f09_g17.LE2-1011.001"
ext8 = ".cam.h0.2015-2065"
ext8f = ".nc"
tit8 = "CESM SSE 7 Tg/yr"
string8 = "{dir:s}/{case:s}/{case:s}{ext:s}.{var:s}"+ext8f

filespec=string8.format(dir=dir8,case=case8,ext=ext8,var='TS')
print(filespec)
file.write(filespec+'\n')
ind = xr.open_mfdataset(filespec)

Varlist = np.array(['RESTOM','FLNT','FSNT','TS','TMQ','PRECT','AEROD_v','CLDLOW','CLDTOT','LWCF','SWCF','TGCLDIWP','TGCLDLWP',
                    'SHFLX','LHFLX','PBLH','PCONVT','PRECC','PRECS'])
Varlist = np.array(['FLNT'])
Varlist = np.array(['TS'])

for Varname in Varlist:
    print()
    print('-------------------------------'+Varname)
    
    V8AY = make_GA_tser(Varname, string8, dir8, case8, ext8)
    V8AY.coords['year'] = V8AY.coords['year'] - V8AY.coords['year'][0]
    
    V1AY = make_GA_tser(Varname, string1, dir1, case1, ext1)
    #time = V1AY.coords['time']
    V1AY.coords['year'] = V1AY.coords['year'] - V1AY.coords['year'][0]
    V1AY
    print('time',V1AY.coords['year'].values)
    #V1AY.plot()
    
    V2AY = make_GA_tser(Varname, string2, dir2, case2, ext2)
    V2AY.coords['year'] = V2AY.coords['year'] - V2AY.coords['year'][0]
    print('time2',V2AY.coords['year'].values)
    #V2AY.plot()

    if True:
        V3AY = make_GA_tser(Varname, string3, dir3, case3, ext3)
        V3AY.coords['year'] = V3AY.coords['year'] - V3AY.coords['year'][0]
        print('time3',V3AY.coords['year'].values)
        
    V4AY = make_GA_tser(Varname, string4, dir4, case4, ext4)
    #time = V1AY.coords['time']
    V4AY.coords['year'] = V4AY.coords['year'] - V4AY.coords['year'][0]
    #print('time',V4AY.coords['year'].values)

    V5AY = make_GA_tser(Varname, string5, dir5, case5, ext5)
    V5AY.coords['year'] = V5AY.coords['year'] - V5AY.coords['year'][0]
    
    V6AY = make_GA_tser(Varname, string6, dir6, case6, ext6)
    V6AY.coords['year'] = V6AY.coords['year'] - V6AY.coords['year'][0]
    
    V7AY = make_GA_tser(Varname, string7, dir7, case7, ext7)
    V7AY.coords['year'] = V7AY.coords['year'] - V7AY.coords['year'][0]
    

    #V2AY.plot()
    
    fig, axes = plt.subplots(nrows=1,
                             #gridspec_kw={'width_ratios': [1]},
                             #subplot_kw={'projection': plotproj},
                             figsize=(6,4.1),
                            )
    fig.set_dpi(300.0)
    
    # two color tables (one for cesm, one for e3sm. Stay away from the smallest indices for icolors as it will be whitish)
    ecolors = plt.cm.viridis(np.linspace(0, 1, 8))
    ccolors = plt.cm.inferno_r(np.linspace(0, 1, 8))
    ecolors = plt.cm.Purples(np.linspace(0, 1, 8))
    ccolors = plt.cm.Oranges(np.linspace(0, 1, 8))
    
    # pert and control linestyes
    pstyle = 'dotted'
    cstyle = 'solid'
    
    V1AY = V1AY - V1AY[0]
    V2AY = V2AY - V2AY[0]
    V3AY = V3AY - V3AY[0]
    V4AY = V4AY - V4AY[0]
    V5AY = V5AY - V5AY[0]
    V6AY = V6AY - V6AY[0]

    V7AY = V7AY - V7AY[0]
    V8AY = V8AY - V8AY[0]


    #axes.text(0.5,1.01,"ax title")
    V1AY.plot.line(ax=axes,label=tit1,color=ecolors[1,:],linestyle=pstyle)
    V2AY.plot.line(ax=axes,label=tit2,color=ecolors[2,:],linestyle=cstyle)
    V3AY.plot.line(ax=axes,label=tit3,color=ecolors[3,:],linestyle=pstyle)
    V4AY.plot.line(ax=axes,label=tit4,color=ccolors[2,:],linestyle=pstyle)
    V5AY.plot.line(ax=axes,label=tit5,color=ccolors[3,:],linestyle=cstyle)
    V6AY.plot.line(ax=axes,label=tit6,color=ccolors[4,:],linestyle=cstyle)
    V7AY.plot.line(ax=axes,label=tit7,color=ccolors[5,:],linestyle=cstyle)
    V8AY.plot.line(ax=axes,label=tit8,color=ccolors[6,:],linestyle=pstyle)


    Vmin = np.minimum(np.minimum(V5AY,V6AY),V7AY)
    Vmax = np.maximum(np.maximum(V5AY,V6AY),V7AY)
    axes.fill_between(Vmin.year.values, Vmin.values, Vmax.values, color=ccolors[3,:], alpha=0.2)

    if Varname == 'TS':
        axlabel = "T$_S$"
    else:
        axlabel = V1AY.long_name
        
    axes.legend(fontsize=8, framealpha=0.)
    axes.set_xlabel('Year beyond 2020')
    axes.set_ylabel(axlabel+' Change from 2020 IC ('+V1AY.units+')')
    plt.savefig(fnprefix+'_'+Varname+'.pdf',format='pdf',dpi=300)
    plt.show()
    
    print('field processing complete')
```

```python

import fsspec

def s3load_zarr(s3path):
    return xr.open_dataset(fsspec.get_mapper(s3path),engine = 'zarr')

var = 'TS'
experiment = "20221014.v2.LR.WCYCLSSP245.E2_CNTL_01"
s3_path = f's3://mcb-zarr/E3SMv2/{experiment}/atm/proc/tseries/month_1/{experiment}.eam.h0.{var}.zarr'
DS1 = s3load_zarr(s3_path)
print('xxx',DS1)
DS1f = center_time(DS1)
DS1f.time
```

```python
DS1f.time_bnds
```

```python
file.close()
```

```python

```

```python

import fsspec

def s3load_zarr(s3path):
    return xr.open_dataset(fsspec.get_mapper(s3path),engine = 'zarr')

experiment = 'b.e21.BSSP245smbb_MCBsse2.5Tg_R1.LE2-1031.002'
var = 'TREFHT'
s3_path = f's3://mcb-zarr/CESM2/{experiment}/atm/proc/tseries/month_1/{experiment}.cam.h0.{var}.zarr'

print('s3_path',s3_path)

## data is a lazy-loaded xarray dataset
data = s3load_zarr(s3_path)
time = data.time
lat = data.lat
print(data.TREFHT)
data
```

```python
def make_GA_tser(Varname, filespec, dirname, casename, extname):
    ''' make global avg from tseries files
        Varname: Variable name
        filespec: the string used to format the filename for a model
        dirname: directory where netcdf file is found
        casename: name of case
        extname: string used to specify years, etc
        returns DataArray with timeseries of global averages
    '''
    ind1 = filespec.format(dir=dirname,case=casename,ext=extname,var=Varname)
    print('opening',ind1)
    file.write(ind1+'\n')
    DS1 = xr.open_mfdataset(ind1)
    DS1 = center_time(DS1)
    month_length = DS1.time.dt.days_in_month
    twgts = month_length.groupby("time.year") / month_length.groupby("time.year").sum()
    #print('twgts',twgts)
    #print('xxx',DS1.time)
    Var1 = xr_getvar(Varname,DS1)


    if 'area' in DS1:
        area = DS1['area']
        print('area on DS1')
    else: 
        if 'ncol' in DS1.dims:
            print ('get CS data')
            #areafile = '~/NetCDF_Files/F2010_PJR1.eam.h0.0001-01.nc'
            areafile = '~/NetCDF_Files/ne30pg2.nc'
            DSA = xr.open_mfdataset(areafile)
            if len(DSA['ncol']) != len(DS1['ncol']):
                raise ValueError('area file mismatch')
            area = DSA.area
        else:
            print('calculating fv area weights')
            lat = Var1['lat'].values
            lon = Var1['lon'].values
            aread = make_fvarea(lon,lat)
            area = xr.DataArray(aread, dims=['lat','lon'], coords={'lon':lon,'lat':lat})
            area.attrs['units']='steradians'
        #print('area',area)
        
    wdims = area.dims
    #print('wdims',wdims)
    weights = area/(area.sum(dim=wdims))
    #print('weights sum',weights.shape,weights.sum(dim=wdims).values)

    #print(Varname, Var1)
    # Global Avg by month
    V1A = Var1.weighted(weights).mean(dim=wdims)
    # calculate ann avgs accounting for number of days in month
    V1AY = (V1A*month_length).groupby("time.year").sum()/month_length.groupby("time.year").sum()
    #print('V1AY.values', V1AY.values)
    return V1AY

regtag = ""
weights = None



```

```python
def center_timex(DS1):
    """center_time(DS1)                                                                                                  
                                                                                                                         
    correct the time coordinate in DS1 to represent the center of the time bounds                                        
                                                                                                                         
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
    time = DS['time'].values
    print('center_time returning DS', time[0],time[1],time[-2],time[-1])
    return DS

def make_GA_ts2(Varname, s3_path ):
    ''' make global avg from tseries files
        Varname: Variable name
        filespec: the string used to format the filename for a model
        dirname: directory where netcdf file is found
        casename: name of case
        extname: string used to specify years, etc
        returns DataArray with timeseries of global averages
    '''
    #print('xx s3_path',s3_path)
    DS1 = s3load_zarr(s3_path)
    #print('xxx DS1',DS1)
    #DS1 = xr.open_mfdataset(ind1)
    DS1 = center_timex(DS1)
    month_length = DS1.time.dt.days_in_month
    twgts = month_length.groupby("time.year") / month_length.groupby("time.year").sum()
    #print('twgts',twgts)
    #print('xxx',DS1.time)
    Var1 = xr_getvar(Varname,DS1)
    area = DS1['area']
    #print('Var1',Var1)
    #1./0.    
    wdims = area.dims
    #print('wdims',wdims)
    weights = area/(area.sum(dim=wdims))
    #print('weights sum',weights.shape,weights.sum(dim=wdims).values)

    #print(Varname, Var1)
    # Global Avg by month
    V1A = Var1.weighted(weights).mean(dim=wdims)
    # calculate ann avgs accounting for number of days in month
    V1AY = (V1A*month_length).groupby("time.year").sum()/month_length.groupby("time.year").sum()
    #print('V1AY.values', V1AY.values)
    return V1AY
```
