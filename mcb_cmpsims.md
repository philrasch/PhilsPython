---
jupyter:
  jupytext:
    formats: ipynb,md
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
from xhistogram.xarray import histogram
%run -i ~/Python/pjr3
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
ind1 = '/global/cscratch1/sd/pjr/E3SMv2/v2.LR.histAMIP_e6/climo/v2.LR.histAMIP_e6_ANN_200001_200012_climo.nc'
pref2 = 'con_'
ind2 = '/global/cscratch1/sd/pjr/E3SMv2/v2.LR.histAMIP_c6/climo/v2.LR.histAMIP_c6_ANN_200001_200012_climo.nc'
regtag = ''
xr.set_options(keep_attrs=True)
# reorder coords so ncol is alway first dim 
# to allow lat/lon interpolation across multiple dimensions
DS1 = xr.open_mfdataset(ind1).transpose('ncol'+regtag,...) 
DS2 = xr.open_mfdataset(ind2).transpose('ncol'+regtag,...) 

tlast = DS1.time[-1]
print('tlast',tlast.values)

lon = DS1['lon'+regtag]#.isel(time=0)#.squeeze()
print('lon',lon.shape,lon.min().values,lon.max().values)
lat = DS1['lat'+regtag]#.isel(time=0)
print('lat',lat.shape,lat.min().values,lat.max().values)



```

```python
def xr_cshplot(xrVar, xrLon, xrLat, plotproj=None, ax=None, cax=None,ylabels=None,clevs=None, cmap=None, title=None):
    """xr_cshplot xarray cubed sphere horizontal plot
    """

    dinc = 1.  # increment of mesh in degrees
    lon_h=np.arange(np.floor(xrLon.min().values),np.ceil(xrLon.max().values+dinc), dinc)
    lat_h=np.arange(np.floor(xrLat.min().values),np.ceil(xrLat.max().values+dinc), dinc)
    xv,yv=np.meshgrid(lon_h,lat_h)
    data_regridded = interp_ap(xv, yv, xrVar.values,xrLat.values,xrLon.values)
    df = data_regridded.flatten()
    dsub = df[np.isfinite(df)] # ignore NaN
    zmax = dsub.max()
    zmin = dsub.min()
    #print('masked interpolated range',zmin,zmax)
    dataproj=ccrs.PlateCarree()    # data is always assumed to be lat/lon
    #plotproj=ccrs.Orthographic(central_latitude=0,central_longitude=55)   # any projections should work
    #print('plotproj is ',plotproj)
    if plotproj is None: plotproj = ccrs.Mercator
    if ax is None: ax = plt.gca()
    if cax is None: cax = ax
    if ylabels is None: ylabels = True
    if clevs is None:
        clevs = findNiceContours(np.array([zmin,zmax]),nlevs=10)
    #print('clevs',clevs)
    if cmap is None:
        #print('aaa, grabbing cmap default')
        cmap = mpl.cm.get_cmap()
        #print('bbb',cmap.N)
    #print('cmap',cmap)
    extend = 'both'
    norm = mpl.colors.BoundaryNorm(clevs,cmap.N,extend=extend)
    #print('norm',norm(clevs))

    plotproj=ccrs.Mollweide(central_longitude=200)   # any projections should work 
    clat = (lat.values.min()+lat.values.max())/2.
    clon = (lon.values.min()+lon.values.max())/2.
    #plotproj=ccrs.NearsidePerspective(central_longitude=clon, central_latitude=clat)
    plotproj=ccrs.Mercator()
    #ax = plt.axes(projection=plotproj)
    #ax.set_extent([lon.values.min(), 260., lat.values.min(), lat.values.max()])
    #ax.set_global()
    pl = ax.contourf(xv, yv, data_regridded, levels=clevs, # vmin=zmin, vmax=zmax,
                     norm=norm, cmap=cmap,
                     extend=extend, transform=ccrs.PlateCarree())
    # Add colorbar to plot
    cb = plt.colorbar(
        pl, orientation='horizontal',ticks=clevs,ax=ax,
        label='%s (%s)'%(xrVar.long_name, xrVar.units), pad=0.1
    )
    if not title is None:
        ax.set_title(title)
        
    cb.ax.tick_params(labelsize=8)
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=2, color='gray', alpha=0.5)
    gl.left_labels=ylabels
    gl.right_labels=ylabels
    ax.coastlines(linewidth=1,color='blue')
    return

    
```

```python
# create a mask to isolate a region of interest
pmask = ((lon > 220) & (lon < 250) & (lat > 15) & (lat < 35))#[0] # select a subregion
pmask = (lon > -999) # select all points
#print('pmask',pmask)
lonsub = DS1['lon'+regtag].where(pmask)#.isel(time=0)
latsub = DS1['lat'+regtag].where(pmask)#.isel(time=0)
print('subreg',lonsub.min().values,lonsub.max().values, latsub.min().values,latsub.max().values)
#print('subreg size',lonsub.shape)
print('shape and size of variables',lonsub.shape, lonsub.size,' number of unmasked cells ',np.count_nonzero(lonsub.notnull().values))
dinc = 1.  # increment of mesh in degrees
lon_h=np.arange(np.floor(lonsub.values.min()),np.ceil(lonsub.values.max()+dinc), dinc)
if (np.abs(lon_h[0]-lon_h[-1]%360) < 0.01): # delete wrap lon for creating zonal average
    print('removing wrap lon')
    lon_h = lon_h[0:-1]
    
lat_h=np.arange(np.floor(latsub.values.min()),np.ceil(latsub.values.max()+dinc), dinc)
xoutm,youtm=np.meshgrid(lon_h,lat_h)

area = xr_getvar('area',DS1,regtag=regtag).where(pmask)
weights = area.fillna(0)

Varlist = np.array(['T','Q','CLOUD','CLDLIQ','ICWMR','CLDICE','RELHUM','NUMICE','NUMLIQ','Mass_bc'])
#                    RESTOM','FLNT','FSNT','TS','TMQ','PRECT','AEROD_v','CLDLOW','CLDTOT','LWCF','SWCF','TGCLDIWP','TGCLDLWP','SHFLX','LHFLX','PBLH','PCONVT','PRECC','PRECS'])
#Varlist = np.array(['TS','TMQ','PRECT'])
#Varlist = np.array(['RESTOM','LWCF','SWCF','FLNT','FSNT'])

for Vname in Varlist:
    print()
    print('-------------------------------')
    V1 = xr_getvar(Vname,DS1,regtag=regtag).where(pmask)
    if V1.min().values == V1.max().values:
        print('constant field skipping plot ')
    else:
        V1 = xr_getvar(Vname, DS1, regtag).where(pmask).squeeze()
        #print('V1xxx ', V1)
        #print('V1',V1.min().values,V1.max().values)
        #print('V1 shape, size, realsize', V1.shape, np.size(V1.values), np.count_nonzero(V1.notnull().values) )
        V2 = xr_getvar(Vname, DS2, regtag).where(pmask).squeeze()
        DV = V1-V2
        print(Vname, V1.attrs['long_name'],'Range V1 and V2 ',V1.min().values, V1.max().values, V2.min().values, V2.max().values)
        V1A = V1.weighted(weights).mean()
        V2A = V2.weighted(weights).mean()
        print('V1A %5.3f' % (V1A.values),' V2A %5.3f' % (V2A.values))
        
        data1 = interp_ap(xoutm, youtm, V1.values,latsub.values,lonsub.values).mean(axis=1).transpose()
        
        data2 = interp_ap(xoutm, youtm, V2.values,latsub.values,lonsub.values).mean(axis=1).transpose()
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
        dlevs = findNiceContours(np.array([datad.min(),datad.max()]),nlevs = 10, rmClev=0.,sym=True)
        dmap = diverge_map()
        plotZMf(datad, lat_h, lev,axesa=axes[2],plotOpt={'clevs':dlevs,'cmap':dmap,'colorbar':"bot",'units':V2.units,'ytop':ytop,'ltitle':pref1+'-'+pref2})

        #print('attribute check on xarray',hasattr(V1,'units'))
        #plt.savefig('test_'+Vname+'.jpg',format='jpg')
        #plt.tight_layout()
        plt.show()
      

```

```python
# create a mask to isolate a region of interest
pmask = ((lon > 220) & (lon < 250) & (lat > 15) & (lat < 35))#[0] # select a subregion
pmask = (lon > -999) # select all points
#print('pmask',pmask)
lonsub = DS1['lon'+regtag].where(pmask)#.isel(time=0)
latsub = DS1['lat'+regtag].where(pmask)#.isel(time=0)
print('subreg',lonsub.min().values,lonsub.max().values, latsub.min().values,latsub.max().values)
#print('subreg size',lonsub.shape)
print('shape and size of variables',lonsub.shape, lonsub.size,' number of unmasked cells ',np.count_nonzero(lonsub.notnull().values))


area = xr_getvar('area',DS1,regtag=regtag).where(pmask)
weights = area.fillna(0)

Varlist = np.array(['RESTOM','FLNT','FSNT','TS','TMQ','PRECT','AEROD_v','CLDLOW','CLDTOT','LWCF','SWCF','TGCLDIWP','TGCLDLWP',
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
        V1 = xr_getvar(Vname, DS1, regtag).where(pmask)
        #print('V1',V1.min().values,V1.max().values)
        #print('V1 shape, size, realsize', V1.shape, np.size(V1.values), np.count_nonzero(V1.notnull().values) )
        V2 = xr_getvar(Vname, DS2, regtag).where(pmask)
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
        M1 = V1.mean(axis=1)
        M2 = V2.mean(axis=1)
        DM = M1-M2

        clevs = findNiceContours(np.array([M1.values,M2.values]),nlevs = 10)
        dlevs = findNiceContours(np.array([DM.min().values,DM.max().values]),nlevs = 20, rmClev=0.,sym=True)
        #dlevs = [-5.,-2.,-1.,-0.5,-0.2,-0.1,0.1,0.2,0.5,1.,2.,5.]
        #print('xxx',dlevs)
        dmap = diverge_map()

        xr_cshplot(M1, lonsub, latsub,ax=axes[0],clevs=clevs,title=pref1)
        xr_cshplot(M2, lonsub, latsub,ax=axes[1],clevs=clevs,ylabels=False,title=pref2)
        xr_cshplot(DM, lonsub, latsub,ax=axes[2],clevs=dlevs,cmap=dmap,title=pref1+'-'+pref2)
        plt.savefig('test_'+Vname+'.jpg',format='jpg')
        plt.show()

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

```python
# next cells not used
help(xr_getvar)
1./0.
```