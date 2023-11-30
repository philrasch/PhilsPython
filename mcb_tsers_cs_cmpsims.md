---
jupyter:
  jupytext:
    formats: ipynb,md
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

**compare two simulations on the ne30 native grid
does zonal averaging, and can focus on a small region**

```python
import sys
print(sys.version)
%matplotlib inline
#from xhistogram.xarray import histogram
%run -i ~/Python/pjr3

ne30area = '~/NetCDF_Files/F2010_PJR1.eam.h0.0001-01.nc'
DSA = xr.open_mfdataset(ne30area)
lon = xr_getvar('lon',DSA)
lat = xr_getvar('lat',DSA)
area = xr_getvar('area',DSA)
```

```python
def getvarDS(Varname,fstring,case_start,case_end,regtag=''):
    """getvar DS
       get variable from file specifying the formatting of the dataset file names
    """
    ind = fstring % (case_start,Varname,case_end)
    print('opening',ind)
    DS = xr.open_mfdataset(ind).chunk({'time': 12}).transpose('ncol'+regtag,...) 
    DS = center_time(DS)
    #DS.coords['lon'] = (DS.coords['lon'] + 180) % 360 - 180
    #DS = DS.sortby(DS.lon)
    
    Var = xr_getvar(Varname,DS)
    #VM = Var.mean(dim='time',keep_attrs=True)
    return Var;

def derfldDS(VN, fstring1, case_start1, case_end1):
    """ calculate some derived fields specifying formatting of file names
    RESTOM is Top of Model Residual energy imbalance (In - Out)
    PRECT is total precip
    """
    if VN == 'RESTOM':    
        FSNT = getvarDS('FSNT', fstring1, case_start1, case_end1)
        FLNT = getvarDS('FLNT', fstring1, case_start1, case_end1)
        RESTOM = FSNT - FLNT
        RESTOM = RESTOM.rename('RESTOM')
        RESTOM.attrs['long_name'] = 'TOA imbalance'
        return RESTOM   
    elif VN == 'PRECT':
        PRECL = getvarDS('PRECL', fstring1, case_start1, case_end1)
        PRECC = getvarDS('PRECC', fstring1, case_start1, case_end1)
        PRECT = PRECL+PRECC
        PRECT = PRECT.rename('PRECT')
        PRECT.attrs['long_name'] = 'total precipitation (liq + ice)'
        return PRECT
    else:
        return getvarDS(VN, fstring1, case_start1, case_end1)
    

```

```python

case_start2 = "/e3sm_prod/phil/climo/cesm/Fixed_SST/ne30pg2/Fixed_SST.cam.h0.1-20."
case_end2 = ".nc"
pref2='CESMcontrol'
fstring2 ='%s%s%s' 

case_start1 = "/e3sm_prod/phil/climo/e3sm/20220930.v2.LR.F2010.E1_CNTL/ne30pg2/20220930.v2.LR.F2010.E1_CNTL.eam.h0.1-14."
case_end1 = ".nc"
fstring1 ='%s%.0s%.0s' 
fstring1 ='%s%s%s' 
pref1='E3SMcontrol'

case_start2 = '/e3sm_prod/e3sm-reshaped/20220930.v2.LR.F2010.E1_CNTL/atm/20220930.v2.LR.F2010.E1_CNTL.eam.h0.1-14.'
case_end2 = '.nc'
pref2='0Tg/yr,R1-3'
fstring2 ='%s%s%s' 

# working from the original files
case_start2 = '/e3sm_prod/mingxuan/archive/20230330.v2.LR.F2010.MCB-SSLT-EM.R1-3.test02/run/*.eam.h0*.nc'
case_end2 = ''
pref2='18Tg/yr,R1-3'
fstring2 ='%s%s%s'
fstring2 ='%s%.0s%.0s'


# working from the time series files
#case_start2 = '/e3sm_prod/mingxuan/archive/20230405.v2.LR.F2010.MCB-SSLT-EM.R1-3.test01/reshaped/20230405.v2.LR.F2010.MCB-SSLT-EM.R1-3.test01.eam.h0.1-6.'
#case_end2 = '.nc'
case_start2 = '/e3sm_prod/phil/tseries/e3sm/20230405.v2.LR.F2010.MCB-SSLT-EM.R1-3.test01/'
case_end2 = '_20230405.v2.LR.F2010.MCB-SSLT-EM.R1-3.test01_000101_002112.nc'
pref2='28Tg/yr,R1-3'
fstring2 ='%s%s%s'

case_start2 = '/e3sm_prod/phil/tseries/e3sm/20230423.v2.LR.F2010.MCB-SSLT-EM.R1-3.test01/'
case_end2 = '_20230423.v2.LR.F2010.MCB-SSLT-EM.R1-3.test01_000101_001112.nc'
pref2='42Tg/yr,R1-3'
fstring2 ='%s%s%s'

if False:
    # working from the time series files
    #case_start2 = '/e3sm_prod/mingxuan/archive/20230330.v2.LR.F2010.MCB-SSLT-EM.R1-3.test02/reshaped/20230330.v2.LR.F2010.MCB-SSLT-EM.R1-3.test02.eam.h0.1-11.'
    #case_end2 = '.nc'
    case_start2 = '/e3sm_prod/phil/tseries/e3sm/20230330.v2.LR.F2010.MCB-SSLT-EM.R1-3.test02/'
    case_end2 = '_20230330.v2.LR.F2010.MCB-SSLT-EM.R1-3.test02_000101_001112.nc'
    pref2='18Tg/yr,R1-3'
    fstring2 ='%s%s%s'

#
```

```python
regtag = ""
weights = None

Varlist = np.array(['RESTOM','FLNT','FSNT','TS','TMQ','PRECT','AEROD_v','CLDLOW','CLDTOT','LWCF','SWCF','TGCLDIWP','TGCLDLWP',
                    'SHFLX','LHFLX','PBLH','PCONVT','PRECC','PRECS'])
Varlist = np.array(['RESTOM','TMQ','AODVIS','FLNT','FSNT','FSNTC','FLNTC','TS','PRECC','PRECL','CLDLOW','CLDHGH','CLDTOT','LWCF','SWCF','TGCLDIWP','TGCLDLWP',
                    'SHFLX','LHFLX','PBLH'])
#Varlist = np.array(['TS','TMQ','PRECT'])
#Varlist = np.array(['RESTOM','LWCF','SWCF','FLNT','FSNT'])
#Varlist = np.array(['AEROD_v'])
#Varlist = np.array(['CLDLOW','TGCLDLWP','PRECL','PRECC','SWCF'])
#Varlist = np.array(['FSNT','FSNTC','FLNT','FLNTC'])
#Varlist = np.array(['SWCF','TGCLDLWP','CLDLOW'])
#Varlist = np.array(['FSUTOA','FSNT','FSNTC'])
#Varlist = np.array(['CLDHGH','ICEFRAC'])

# specify regions (assume lon always specified as west, then east limit)
xreg = np.array([[-150.,-110.],[-110,-70],[-25.,15.],[170.,-120.],[-170.,-90.]])%360.
yreg = np.array([[0.,30.],     [-30.,0.], [-30.,0.], [30.,50.],   [-50.,-30.] ])
namereg = ['NEP','SEP','SEA','NP','SP']
#xreg = [[0.,360.]]
#yreg = [[-90.,91.]]

if False:
    ind = fstring2 % (case_start2,'FLNT',case_end2)
    print('opening',ind)
    DS = xr.open_mfdataset(ind).chunk({'time': 12}).transpose('ncol'+regtag,...) 
    #DS
    #print('DS1',DS)
    #print('DS1.time',DS.time)
    #print('DDD',DS.indexes['time'])
    ##XX0 = DS['FLNT']
    print('xx0',XX0)
    print('xx0a',XX0.time)
    DS = center_time(DS)
    print('DS2',DS)
    print('DS2.time',DS.time)

for Varname in Varlist:
    print()
    print('-------------------------------'+Varname)
    Var1 = derfldDS(Varname, fstring1, case_start1, case_end1)
    Var1y = tavg_mon_wt(Var1)
    Var1yga = Var1y.weighted(area).mean('ncol',keep_attrs=True)
    V1 = Var1yga

    Var2 = derfldDS(Varname, fstring2, case_start2, case_end2)
    #rint('xx1',Var2.time.units)
    Var2y = tavg_mon_wt(Var2)
    time = Var2y.time.values
    #print('time',Var2y.time)
    #print('Var2y time',time)
    Var2yga = Var2y.weighted(area).mean('ncol',keep_attrs=True)
    V2 = Var2yga
    
    DV = V2-V1.values
    DVM = DV.mean().values
    DVM2 = DV[1:6].mean().values
    #print('xxx',DVM2, time[1:6])
    DVM3 = DV[1:11].mean().values

    print('DV range, mean', DV.min().values, DV.max().values, DVM)
    print(Varname, V1.attrs['long_name'],'Range V1 and V2 ',V1.min().values, V1.max().values, V2.min().values, V2.max().values)
    DV.plot()
    print('time yrs', time)
    plt.plot([time[0],time[-1]],[DVM,DVM])
    plt.plot([time[1],time[6]],[DVM2,DVM2])
    plt.plot([time[1],time[-1]],[DVM3,DVM3])

    #print('xxx',time[0])
    plt.show()
    
    print('field processing complete')

```

```python
Varname = 'FSNT'
Var2 = derfldDS(Varname, fstring2, case_start2, case_end2)
Var2y = tavg_mon_wt(Var2)# .isel(time=slice(0,9))
Var2yga = Var2y.weighted(area).mean('ncol')
FSNT = Var2yga
FSNT.plot()
Varname = 'FLNT'
Var2 = derfldDS(Varname, fstring2, case_start2, case_end2)
Var2y = tavg_mon_wt(Var2)
Var2yga = Var2y.weighted(area).mean('ncol')
FLNT = Var2yga
FLNT.plot()
plt.show()
```

```python
1./0.
```

```python
# process file holding data only over a region at every timestep of a model run

pref1 = 'exp_'
ind1 = '/global/cscratch1/sd/pjr/E3SMv2/v2.LR.histAMIP_x4/tests/M_1x10_ndays/run/v2.LR.histAMIP_x4.eam.h1.2015-07-*0.nc'
pref2 = 'con_'
ind2 = '/global/cscratch1/sd/pjr/E3SMv2/v2.LR.histAMIP_e4/tests/M_1x10_ndays/run/v2.LR.histAMIP_e4.eam.h1.2015-07-*0.nc'
regtag = '_190e_to_250e_0n_to_35n'
pref1 = 'exp_'
ind1 = '/global/cscratch1/sd/pjr/E3SMv2/v2.LR.histAMIP_e5/tests/L_1x1_nmonths/run/v2.LR.histAMIP_e5.eam.h0.2015-07.nc'
pref2 = 'con_'
ind2 = '/global/cscratch1/sd/pjr/E3SMv2/v2.LR.histAMIP_c5/tests/L_1x1_nmonths/run/v2.LR.histAMIP_c5.eam.h0.2015-07.nc'
pref1 = 'exp_'
ind1 = '/global/cscratch1/sd/pjr/E3SMv2/v2.LR.histAMIP_e6/run/v2.LR.histAMIP_e6.eam.h0.2000-07.nc'
pref2 = 'con_'
ind2 = '/global/cscratch1/sd/pjr/E3SMv2/v2.LR.histAMIP_c6/run/v2.LR.histAMIP_c6.eam.h0.2000-07.nc'
pref1 = 'exp_'
ind1 = '/global/cscratch1/sd/pjr/E3SMv2/v2.LR.histAMIP_e6/climo/v2.LR.histAMIP_e6_ANN_200001_200412_climo.nc'
pref2 = 'con_'
ind2 = '/global/cscratch1/sd/pjr/E3SMv2/v2.LR.histAMIP_c6/climo/v2.LR.histAMIP_c6_ANN_200001_200412_climo.nc'
regtag = ''
ind1 = '/e3sm_prod/mingxuan/archive/20230330.v2.LR.F2010.MCB-SSLT-EM.R1-3.test02/run/20230330.v2.LR.F2010.MCB-SSLT-EM.R1-3.test02.eam.h0.00*.nc'
Varname = 'SWCF'
ind2 = '/e3sm_prod/e3sm-reshaped/20220930.v2.LR.F2010.E1_CNTL/atm/20220930.v2.LR.F2010.E1_CNTL.eam.h0.1-14.'+Varname+'.nc'

xr.set_options(keep_attrs=True)
# reorder coords so ncol is alway first dim 
# to allow lat/lon interpolation across multiple dimensions
DS1 = xr.open_mfdataset(ind1).chunk({'time': 12}).transpose('ncol'+regtag,...) 
DS1 = center_time(DS1)
DS2 = xr.open_mfdataset(ind2).chunk({'time': 12}).transpose('ncol'+regtag,...) 
DS2 = center_time(DS2)
print(DS2)


lon = DS1['lon'+regtag]#.isel(time=0)#.squeeze()
print('lon',lon.shape,lon.min().values,lon.max().values)
lat = DS1['lat'+regtag]#.isel(time=0)
print('lat',lat.shape,lat.min().values,lat.max().values)



```

```python
Varname='FSNT'
regtag=''
ind2 = '/e3sm_prod/e3sm-reshaped/20220930.v2.LR.F2010.E1_CNTL/atm/20220930.v2.LR.F2010.E1_CNTL.eam.h0.1-14.'+Varname+'.nc'
DS2 = xr.open_mfdataset(ind2).chunk({'time': 12}).transpose('ncol'+regtag,...) 
DS2 = center_time(DS2)
Var2 = xr_getvar(Varname,DS2)
Var2y = tavg_mon_wt(Var2)
area = xr_getvar('area',DS1).isel(time=0)
print('area',area)
ntimes = len(Var2y['time'].values)
print('ntimes',ntimes)
print('Var2y',Var2y)
for yr in np.arange(ntimes):
    V2 = Var2y[:,yr]
    V2a = V2.weighted(area).mean()
    #print(yr, V2a.values)
#Var2yga = Var2y.mean('ncol')
Var2yga = Var2y.weighted(area).mean('ncol')
Var2yga.plot()
plt.show()
```

```python
DS2 = center_time(DS2)
tlast = DS2.time[-1]
tbnds = DS2.time_bnds[-1]
print('tlast','tbnds',tlast.values,tbnds.values)
```

```python
SWCF1 = xr_getvar('SWCF',DS1)
SWCF1 = tavg_mon_wt(SWCF1)
SWCF1
area = xr_getvar('area',DS1)
for yr in np.arange(6):
    print('yr',yr)
    V1 = SWCF1[:,yr]
    V1a = V1.weighted(area).mean()
    print(V1a.values)
```

```python
pmask = ((lon > 220) & (lon < 250) & (lat > 15) & (lat < 35))#[0] # select a subregion
pmask = (lon > -999) # select all points
xxx = pmask.load()
colinds = np.where(xxx.values)[0]
#print('xxx',type(colinds), colinds)
```

```python
lonsub = DS1['lon'+regtag].isel(ncol=colinds)
latsub = DS1['lat'+regtag].isel(ncol=colinds)
print('subreg',lonsub.min().values,lonsub.max().values, latsub.min().values,latsub.max().values)
#print('subreg size',lonsub.shape)
print('shape and size of variables',lonsub.shape, lonsub.size,' number of unmasked cells ',np.count_nonzero(lonsub.notnull().values))
dinc = 1.  # increment of mesh in degrees
lon_h=np.arange(np.floor(lonsub.min().values),np.ceil(lonsub.max().values+dinc), dinc)
if (np.abs(lon_h[0]-lon_h[-1]%360) < 0.01): # delete wrap lon for creating zonal average
    print('removing wrap lon')
    lon_h = lon_h[0:-1]
    
lat_h=np.arange(np.floor(latsub.min().values),np.ceil(latsub.max().values+dinc), dinc)
xoutm,youtm=np.meshgrid(lon_h,lat_h)
print('xxx',xoutm.shape,xoutm.min(),xoutm.max(),youtm.min(),youtm.max())
area = xr_getvar('area',DS1,regtag=regtag).isel(ncol=colinds)
wtsh = area.fillna(0)
#print('wtsh',wtsh)
DPOG1 = xr_getvar('DPOG',DS1,regtag=regtag).isel(ncol=colinds)
weights1 = wtsh*DPOG1
weights1 = weights1.fillna(0)
DPOG2 = xr_getvar('DPOG',DS2,regtag=regtag).isel(ncol=colinds)
weights2 = wtsh*DPOG2
weights2 = weights2.fillna(0)
#print('weights',weights2)

Varlist = np.array(['T','Q','CLOUD','CLDLIQ','ICWMR','CLDICE','RELHUM','NUMICE','NUMLIQ','Mass_bc'])
#                    RESTOM','FLNT','FSNT','TS','TMQ','PRECT','AEROD_v','CLDLOW','CLDTOT','LWCF','SWCF','TGCLDIWP','TGCLDLWP','SHFLX','LHFLX','PBLH','PCONVT','PRECC','PRECS'])
Varlist = np.array(['T'])
#Varlist = np.array(['RESTOM','LWCF','SWCF','FLNT','FSNT'])

for Vname in Varlist:
    print()
    print('-------------------------------')
    V1 = xr_getvar(Vname,DS1,regtag=regtag).isel(ncol=colinds)
    if V1.min().values == V1.max().values:
        print('constant field skipping plot ')
    else:
        #V1 = xr_getvar(Vname, DS1, regtag).where(pmask).squeeze()
        V1 = xr_getvar(Vname, DS1, regtag).isel(ncol=colinds).squeeze()
        V2 = xr_getvar(Vname, DS2, regtag).isel(ncol=colinds).squeeze()
        DV = V1-V2
        print(Vname, V1.attrs['long_name'],'Range V1 and V2 ',V1.min().values, V1.max().values, V2.min().values, V2.max().values)
        V1A = V1.weighted(weights1).mean()
        V2A = V2.weighted(weights2).mean()
        print('mass weight average: V1A %5.3f' % (V1A.values),' V2A %5.3f' % (V2A.values))
        # create regular lat/lon gridding to make zonal averages. Use dataarray to make NaNs easier to process
        Vnew1 = xr.DataArray(interp_ap(xoutm, youtm, V1.values,latsub.values,lonsub.values), 
                            coords={'lat': lat_h,'lon': lon_h,'lev': V1.lev.values}, 
                            dims=["lat", "lon","lev"])
        Vnew1_xa = Vnew1.mean(dim='lon')
        data1 = Vnew1_xa.values.transpose()
        
        Vnew2 = xr.DataArray(interp_ap(xoutm, youtm, V2.values,latsub.values,lonsub.values), 
                            coords={'lat': lat_h,'lon': lon_h,'lev': V2.lev.values}, 
                            dims=["lat", "lon","lev"])
        Vnew2_xa = Vnew2.mean(dim='lon')
        data2 = Vnew2_xa.values.transpose()
        
        Vnewd_xa = Vnew1_xa - Vnew2_xa
        dmin = Vnewd_xa.min().values
        dmax = Vnewd_xa.max().values
        print('dmin,dmax',dmin,dmax)

        datad = data1-data2
#       data1 = data1.mean(axis=1).transpose()
        lev = V1['lev'].values
#        plotZMf(data1, lat_h, lev)
        fig, axes = plt.subplots(ncols=3
                                 ,gridspec_kw={'width_ratios': [1, 1, 1]}
#                                 ,subplot_kw={'projection': ccrs.PlateCarree()}
                                 ,figsize=(16,5)
                                )
        ytop = 1.
        plotZMf(data1, lat_h, lev,axesa=axes[0],plotOpt={'colorbar':"botnd",'units':V1.units,'ltitle':pref1,'ytop':ytop})
        plotZMf(data2, lat_h, lev,axesa=axes[1],plotOpt={'colorbar':"bot",'units':V2.units,'ltitle':pref2,'rtitle':V2.long_name,'ytop':ytop})
        dlevs = findNiceContours(np.array([dmin,dmax]),nlevs = 10, rmClev=0.,sym=True)
        dmap = diverge_map()
        plotZMf(datad, lat_h, lev,axesa=axes[2],plotOpt={'clevs':dlevs,'cmap':dmap,'colorbar':"bot",'units':V2.units,'ytop':ytop,'ltitle':pref1+'-'+pref2})

        #print('attribute check on xarray',hasattr(V1,'units'))
        #plt.savefig('test_'+Vname+'.jpg',format='jpg')
        #plt.tight_layout()
        plt.show()
      

```

```python

print('subreg',lonsub.min().values,lonsub.max().values, latsub.min().values,latsub.max().values)
#print('subreg size',lonsub.shape)
print('shape and size of variables',lonsub.shape, lonsub.size,' number of unmasked cells ',np.count_nonzero(lonsub.notnull().values))

area = xr_getvar('area', DS1, regtag).isel(ncol=colinds).squeeze()

weights = area.fillna(0)

Varlist = np.array(['RESTOM','FLNTC','FLNT','FSNTC','FSNT','TS','TMQ','PRECT','AEROD_v','CLDLOW','CLDTOT','LWCF','SWCF','TGCLDIWP','TGCLDLWP',
                    'SHFLX','LHFLX','PBLH','PCONVT','PRECC','PRECS'])
#Varlist = np.array(['TS','TMQ','PRECT'])
#Varlist = np.array(['RESTOM','LWCF','SWCF','FLNT','FSNT'])
Varlist = np.array(['AEROD_v'])


for Vname in Varlist:
    print()
    print('-------------------------------')
    V1 = xr_getvar(Vname,DS1,regtag=regtag).where(pmask)
    if V1.min().values == V1.max().values:
        print('constant field skipping plot ')
    else:
        V1 = xr_getvar(Vname, DS1, regtag).isel(ncol=colinds).squeeze()
        V2 = xr_getvar(Vname, DS2, regtag).isel(ncol=colinds).squeeze()
        DV = V1-V2
        print(Vname, V1.attrs['long_name'],'Range V1 and V2 ',V1.min().values, V1.max().values, V2.min().values, V2.max().values)
        V1A = V1.weighted(weights).mean()
        V2A = V2.weighted(weights).mean()
        print('V1A %5.3f' % (V1A.values),' V2A %5.3f' % (V2A.values))
        fig, axes = plt.subplots(ncols=3
                                 ,gridspec_kw={'width_ratios': [1, 1, 1]}
                                 ,subplot_kw={'projection': ccrs.PlateCarree()}
                                 ,figsize=(16,5)
                                )

        clevs = findNiceContours(np.array([V1.values,V2.values]),nlevs = 10)
        dlevs = findNiceContours(np.array([DV.min().values,DV.max().values]),nlevs = 20, rmClev=0.,sym=True)
        #dlevs = [-5.,-2.,-1.,-0.5,-0.2,-0.1,0.1,0.2,0.5,1.,2.,5.]
        #print('xxx',dlevs)
        dmap = diverge_map()

        xr_cshplot(V1, lonsub, latsub,ax=axes[0],clevs=clevs,title=pref1)
        xr_cshplot(V2, lonsub, latsub,ax=axes[1],clevs=clevs,ylabels=False,title=pref2)
        xr_cshplot(DV, lonsub, latsub,ax=axes[2],clevs=dlevs,cmap=dmap,title=pref1+'-'+pref2)
        plt.savefig('test_'+Vname+'.jpg',format='jpg')
        plt.show()

```

```python
# next cells not used
help(xr_getvar)
1./0.
```

```python
# example of locating max of a multidimensional array
a = np.array([[1,2,3],[4,3,1]])
a = np.array([[1,4,3,np.nan],[np.nan,4,3,1]])
af = a.flatten()
afd = af[np.isfinite(af)]
print('a',a)
print('len a, af, afd',len(a), len(af), len(afd))
print('afd max',afd.max())
i,j = np.where(a==a.max())
print('i,j',i,j)
i,j = np.where(a==afd.max())
print('i,j',i,j)
a[i,j]
```

```python
from matplotlib import pylab as plt
import numpy

fig = plt.figure()
ax = fig.add_subplot(111)
ax.grid()
# set labels and font size
ax.set_xlabel('X axis', fontsize = 12)
ax.set_ylabel('Y axis', fontsize = 12)

ax.plot(np.random.random(100))

# change font size for x axis
#ax.xaxis.get_label().set_fontsize(20)

plt.show()
ax.xaxis.get_label().get_fontsize()
```
