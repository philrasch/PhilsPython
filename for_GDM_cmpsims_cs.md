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
    display_name: Python [conda env:.conda-pjrpy3] *
    language: python
    name: conda-env-.conda-pjrpy3-py
---

**compare two simulations on the ne30 native grid
does zonal averaging, and can focus on a small region**

```python
import sys
print(sys.version)
%matplotlib inline
%run -i ~/Python/pjr3
from jupytext.config import find_jupytext_configuration_file
print('jupytext config file is ',find_jupytext_configuration_file('.'))


ne30area = '~/NetCDF_Files/F2010_PJR1.eam.h0.0001-01.nc'
DSA = xr.open_mfdataset(ne30area)
lon = xr_getvar('lon',DSA)
lat = xr_getvar('lat',DSA)
area = xr_getvar('area',DSA)
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
pref2='28Tgpyr,R1-3'
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
Varlist = np.array(['RESTOM','FLNTC','FLNT','FSNTC','FSNT','TS','TMQ','PRECT','AODVIS','CLDLOW','CLDTOT','LWCF','SWCF','TGCLDIWP','TGCLDLWP',
                    'SHFLX','LHFLX','PBLH'])
#Varlist = np.array(['TS','TMQ','PRECT'])
#Varlist = np.array(['RESTOM','LWCF','SWCF','FLNT','FSNT'])
#Varlist = np.array(['AODVIS'])


for Varname in Varlist:
    print()
    print('-------------------------------'+Varname)
    Var1 = derfldDS(Varname, fstring1, case_start1, case_end1)
    Var1y = tavg_mon_wt(Var1)
    V1 = Var1y.mean('time')
    Var1yga = V1.weighted(area).mean('ncol',keep_attrs=True)
    Var2 = derfldDS(Varname, fstring2, case_start2, case_end2)
    #rint('xx1',Var2.time.units)
    Var2y = tavg_mon_wt(Var2)
    V2 = Var2y.mean('time')
    Var2yga = V2.weighted(area).mean('ncol',keep_attrs=True)

    DV = V2-V1
    print('DV range, mean', DV.min().values, DV.max().values, (Var2yga-Var1yga).values)
    print()
    print(Varname, V1.attrs['long_name'],'Range V1 and V2 ',V1.min().values, V1.max().values, V2.min().values, V2.max().values)

    print('V1A %5.3f' % (Var1yga.values),' V2A %5.3f' % (Var2yga.values))


    clevs = findNiceContours(np.array([V1.values,V2.values]),nlevs = 10)
    dlevs = findNiceContours(np.array([DV.min().values,DV.max().values]),nlevs = 20, rmClev=0.,sym=True)
    #dlevs = [-5.,-2.,-1.,-0.5,-0.2,-0.1,0.1,0.2,0.5,1.,2.,5.]
    #print('xxx',dlevs)
    dmap = diverge_map()
    sV1A = ' (%5.2f)' % Var1yga.values
    sV2A = ' (%5.2f)' % Var2yga.values
    sDVA = ' (%5.2f)' % (Var2yga-Var1yga).values

    if False:
        fig, axes = plt.subplots(ncols=3
                         ,gridspec_kw={'width_ratios': [1, 1, 1]}
                         ,subplot_kw={'projection': ccrs.PlateCarree()}
                         ,figsize=(16,5)
                        )
        plt.savefig('test_'+Vname+'.jpg',format='jpg')
        plt.show()

    plconf = '3-1x1'
    plconf = '1x3'

    # good setup for 1 row of 3 columns
    if plconf == '1x3':
        plotproj = ccrs.Mollweide()
        plotproj._threshold /= 100.
        fig, axes = plt.subplots(ncols=3
                                 ,gridspec_kw={'width_ratios': [1, 1, 1]}
                                 ,subplot_kw={'projection': plotproj}
                                 ,figsize=(16,5)
                                )

#        xr_llhplot(V1, ax=axes[0],clevs=clevs,title=pref1+sV1A)
        xr_cshplot(V1, lon, lat,ax=axes[1],clevs=clevs,ylabels=False,title=pref1+sV1A)
#        xr_llhplot(V2, ax=axes[1],clevs=clevs,ylabels=False,title=pref2+sV2A)
        xr_cshplot(V2, lon, lat,ax=axes[0],clevs=clevs,title=pref2+sV2A)
#        xr_llhplot(DV, ax=axes[2],clevs=dlevs,cmap=dmap,title=pref1+'-'+pref2+sDVA)        
        xr_cshplot(DV, lon, lat,ax=axes[2],clevs=dlevs,cmap=dmap,title=pref2+'-'+pref1+sDVA)
        pltllbox([-150.,-110.],[0.,30.])
        pltllbox([-110.,-70.],[-30.,0.])
        pltllbox([-25.,15.],[-30.,0.])

        #plt.savefig(pref1+'_'+Varname+'.pdf',format='pdf',dpi=300)
        plt.show()

    # good setup for 3 rows of 1 columns
    if plconf == '3-1x1':

        fig, axes = setfig3b1x1()
        xr_cshplot(V2, lon, lat,ax=axes,clevs=clevs,ylabels=False,title=pref2+sV2A)

        pltllbox([-150.,-110.],[0.,30.])
        pltllbox([-110.,-70.],[-30.,0.])
        pltllbox([-25.,15.],[-30.,0.])
        plt.savefig(pref2+'_'+Varname+'.pdf',format='pdf',dpi=300)
        plt.show()
        
        fig, axes = setfig3b1x1()
        print('V1XXX',V1)
        xr_cshplot(V1, lon, lat,ax=axes,clevs=clevs,title=pref1+sV1A)

        pltllbox([-150.,-110.],[0.,30.])
        pltllbox([-110.,-70.],[-30.,0.])
        pltllbox([-25.,15.],[-30.,0.])
        plt.savefig(pref1+'_'+Varname+'.pdf',format='pdf',dpi=300)
        plt.show()

        fig, axes = setfig3b1x1()
        xr_cshplot(DV, lon, lat,ax=axes,clevs=dlevs,cmap=dmap,title=pref1+'-'+pref2+sDVA)

        pltllbox([-150.,-110.],[0.,30.])
        pltllbox([-110.,-70.],[-30.,0.])
        pltllbox([-25.,15.],[-30.,0.])

        plt.savefig(pref1+'_'+Varname+'-D.pdf',format='pdf',dpi=300)
        plt.show()


        
    print('field processing complete')


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
