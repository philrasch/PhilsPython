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
    display_name: Python [conda env:.conda-pjrpy3] *
    language: python
    name: conda-env-.conda-pjrpy3-py
---

**compare two simulations from data on the ne30 native grid
does zonal averaging, and horizontal plots**

```python
import sys
print(sys.version)
%matplotlib inline
%run -i ~/Python/pjr3
```

```python
# experiment is case1
# control is case2
pref1 = 'exp_'
pref2 = 'con_'
regtag = ''
ind1 = '/scratch2/ec2-user/PJR/E3SM/20221123.v2.LR.F2010.E1_R1-3_CDNC2000/ne30pg2/20221123.v2.LR.F2010.E1_R1-3_CDNC2000.eam.h0.1-11.T.nc' 
ind2 = '/scratch2/ec2-user/PJR/E3SM/20220930.v2.LR.F2010.E1_CNTL/ne30pg2/20220930.v2.LR.F2010.E1_CNTL.eam.h0.1-14.T.nc'

xr.set_options(keep_attrs=True)
# make sure horizontal dimensions ordered first to facilitate interpolations
DS1 = xr.open_mfdataset(ind1).transpose('ncol'+regtag,...) 
DS2 = xr.open_mfdataset(ind2).transpose('ncol'+regtag,...) 


```

**scan globe for subregion mask**

```python

lon = DS1['lon'+regtag]#.isel(time=0)#.squeeze()
print('lon',lon.shape,lon.min().values,lon.max().values)
lat = DS1['lat'+regtag]#.isel(time=0)
print('lat',lat.shape,lat.min().values,lat.max().values)

pmask = ((lon > 220) & (lon < 250) & (lat > 15) & (lat < 35))#[0] # select a subregion
# pmask = (lon > -999) # select all points
xxx = pmask.load()
colinds = np.where(xxx.values)[0]
#print('xxx',type(colinds), colinds)
```

**select subregion, calculate weights**

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
# weights for horizontal (area) averages
wtsh = area.fillna(0)
#print('wtsh',wtsh)
# weights for 3D fields (mass and area weights) DPOG is delta-p over g
DPOG1 = xr_getvar('DPOG',DS1,regtag=regtag).isel(ncol=colinds)
weights1 = wtsh*DPOG1
weights1 = weights1.fillna(0)
DPOG2 = xr_getvar('DPOG',DS2,regtag=regtag).isel(ncol=colinds)
weights2 = wtsh*DPOG2
weights2 = weights2.fillna(0)
#print('weights',weights2)
```

**lat/lon plots on model surfaces**

```python

print('subreg',lonsub.min().values,lonsub.max().values, latsub.min().values,latsub.max().values)
#print('subreg size',lonsub.shape)
print('shape and size of variables',lonsub.shape, lonsub.size,' number of unmasked cells ',np.count_nonzero(lonsub.notnull().values))

area = xr_getvar('area', DS1, regtag).isel(ncol=colinds).squeeze()

weights = area.fillna(0)

Varlist = np.array(['RESTOM','FLNTC','FLNT','FSNTC','FSNT','TS','TMQ','PRECT','AEROD_v','CLDLOW','CLDTOT','LWCF','SWCF','TGCLDIWP','TGCLDLWP',
                    'SHFLX','LHFLX','PBLH','PCONVT','PRECC','PRECS','PS'])
#Varlist = np.array(['TS','TMQ','PRECT'])
Varlist = np.array(['T'])
mylev = 850.
#Varlist = np.array(['PS'])
Varlist = np.sort(Varlist)

for Vname in Varlist:
    print()
    print('-------------------------------')
    V1 = xr_getvar(Vname, DS1, regtag).isel(ncol=colinds).squeeze()
    V2 = xr_getvar(Vname, DS2, regtag).isel(ncol=colinds).squeeze()
    dlist = list(V1.dims)
    if 'lev' in dlist:
        V1 = V1.sel(lev=mylev,method='nearest')
        V2 = V2.sel(lev=mylev,method='nearest')
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

**latitude eta (lon is averaged)**

```python
Varlist = np.array(['T','Q','CLOUD','CLDLIQ','ICWMR','CLDICE','RELHUM','NUMICE','NUMLIQ','Mass_bc'])
#                    RESTOM','FLNT','FSNT','TS','TMQ','PRECT','AEROD_v','CLDLOW','CLDTOT','LWCF','SWCF','TGCLDIWP','TGCLDLWP','SHFLX','LHFLX','PBLH','PCONVT','PRECC','PRECS'])
Varlist = np.array(['T'])

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
