---
jupyter:
  jupytext:
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.16.4
  kernelspec:
    display_name: pjrpy3
    language: python
    name: pjrpy3
---

```python
import sys
print(sys.version)
%matplotlib inline
%run -i ./pjr3.py
```

```python
def findNiceContlog(cmax,ndec=2):
    """ find a reasonable set of contours spaced logarithmically"""
    decade = int(np.log10(clevs[-1])+1)
    mylist = []
    for i in range(decade-ndec,decade):
        cdec = np.array([1.,2.,5.])*10.**i
        mylist.extend(cdec)
    ind = np.array(np.where(mylist > np.abs(cmax)))
    #print('ind',ind,ind[0],len(ind[0]))
    if len(ind[0]) == 1: 
        ind = ind.item()
        mylist = mylist[0:ind]
    #print('indf',ind)

    #if mylist[-1] > np.abs(cmax):
    #   mylist = mylist[0:-2] 
    flist = np.array([-np.array(mylist[::-1]),np.array(mylist)]).flatten()
    return flist
```

```python

case_start1 = "/global/cfs/projectdirs/m1199/huoyilin/forPhilRasch/20230821.F2010.arcticx4v1pg2_ARRM10to60E2r1.pm-cpu.eam.h0.*.nc"        
case_end1 = ""

fstring1 ='%s%.0s%.0s' 
#fstring1 ='%s%s%s' 
Varname = ''

ind1 = fstring1 % (case_start1,Varname,case_end1)
print('opening ind1 ',ind1)

DS1 = xr.open_mfdataset(ind1)
DS1 = fix_DS(DS1)
time = DS1.time

area = DS1.area
lon = DS1.lon
lat = DS1.lat
if len(time) > 0:
    area = area.isel(time=0).squeeze()
    lon = lon.isel(time=0).squeeze()
    lat = lat.isel(time=0).squeeze()
    
weights = area/area.sum(dim='ncol')

print(area.shape)

```

```python
# figure out the distance between points near the NP
northind = lat.argmax()
dinc = gcd(lon.values,lat.values,lon[northind].values,lat[northind].values)
dincs = np.sort(dinc)[1].item() # distance from point near NP to next closest point (km)
distnpsp = np.pi*6371
dincr = 180.*dincs/distnpsp # radial distance between nearest points
dincr
```

```python
dinc = dincr
#dinc = 1.
```

```python

month_length = DS1.time.dt.days_in_month
    #twgts = month_length.groupby("time.year") / month_length.groupby("time.year").sum()
    # Calculate the weights by grouping by 'time.season'.
twgts = (
        month_length.groupby("time.season") / month_length.groupby("time.season").sum()
        )


def get_NHDJF(Varname, DS):
    """Get the variable called Varname from DS
       construct DJF average and mask values north of 60N
    """
    V1 = xr_getvar(Varname,DS)
    V1s = (V1 * twgts).groupby("time.season").sum(dim="time")
    V1s = V1s.where(lat > 60.)
    V1s = V1s.sel(season='DJF')
    #print('yyy',V1s)
    return V1s

```

```python

case_start2 = "/pscratch/sd/s/seefeldt/_forPhil/R2302aRBRcaaa01i/multi_yr/R2302aRBRcaaa01i.wrf.h.201909_202008-DJF.nc"
case_end2 = ""

fstring2 ='%s%.0s%.0s' 
#fstring1 ='%s%s%s' 
Varname = 'ps'

ind2 = fstring2 % (case_start2,Varname,case_end2)
print('opening ind2 ',ind2)

DS2 = xr.open_mfdataset(ind2)
time = DS2.time
Var2 = DS2[Varname].squeeze()
```

```python
Var2.plot()
```

```python
def xr_pshplot(xrVar, xrLon, xrLat, dinc=None, plotproj=None, ax=None, cax=None,ylabels=None,clevs=None, cmap=None, title=None,cbar='default'):
    """xr_pshplot xarray WRF polar stereohorizontal plot
    """

    if dinc is None: dinc = 1.  # increment of mesh in degrees
    #lon_h=np.arange(np.floor(xrLon.min().values),np.ceil(xrLon.max().values+dinc), dinc)
    #lat_h=np.arange(np.floor(xrLat.min().values),np.ceil(xrLat.max().values+dinc), dinc)
    lon_h = np.arange(0.,360+dinc,dinc)-180.
    #lat_h = np.arange(-90.,90+dinc,dinc)
    latn = np.floor((180.+dinc)/dinc).astype(int)
    dincns = 180./latn
    lat_h = np.linspace(-90.,90.,latn)
    #print('lon_h',lon_h)
    #print('lat_h',lat_h)
    xv,yv=np.meshgrid(lon_h,lat_h)
    xrvf = xrVar.values.flatten()
    xrlat = xrLat.values.flatten()
    xrlon = xrLon.values.flatten()
    #for i in range(0, 11):
        #print(i,xrvf[i],xrlat[i],xrlon[i])
    data_regridded = interp_ap(xv, yv, xrvf, xrlat, xrlon)
    #inds = np.where(yv < 60.)
    #data_regridded[inds] = np.nan
    #print('data_regridded',data_regridded.shape)
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
        #cmap = mpl.cm.get_cmap()
        cmap = plt.get_cmap()
        #print('bbb',cmap.N)
    #print('cmap',cmap)
    extend = 'both'
    norm = mpl.colors.BoundaryNorm(clevs,cmap.N,extend=extend)
    #print('norm',norm(clevs))

    pl = ax.contourf(xv, yv, data_regridded, levels=clevs, # vmin=zmin, vmax=zmax,
                     norm=norm, cmap=cmap,
                     extend=extend, transform=ccrs.PlateCarree())
    if cbar == 'default':
        # Add colorbar to plot
        cb = plt.colorbar(
            pl, orientation='horizontal',ticks=clevs,ax=ax,
            label='%s (%s)'%(xrVar.long_name, xrVar.units), pad=0.1
        )
        cb.ax.tick_params(labelsize=8)
        
    if not title is None:
        ax.set_title(title)
        
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,
                      linewidth=2, color='gray', alpha=0.5)
    gl.left_labels=ylabels
    gl.right_labels=ylabels
    ax.coastlines(linewidth=1,color='blue')
    return pl

```

```python
def xr_intum2ll(xrVar, lon_o, lat_o, xrLon=None, xrLat=None):
    """xr_intum2ll: interpolate xarray data on unstructured mesh to lat/lon
    """

    if xrLon is None and xrLat is None:
        #print('assuming lon and lat are part of xarray DataArray')
        #print('xrVar', xrVar)
        xrlon = xrVar.lon.values.flatten()
        xrlat = xrVar.lat.values.flatten()
    else:
        xrlon = xrLon.values
        xrlat = xrLat.values
    #print('xrlon',xrlon.shape)
    #print('xrlat',xrlat.shape)
    xv,yv=np.meshgrid(lon_o,lat_o)
    xrvf = xrVar.values.flatten()
    #xrlat = xrLat.values.flatten()
    #xrlon = xrLon.values.flatten()
    #for i in range(0, 11):
        #print(i,xrvf[i],xrlat[i],xrlon[i])
    data_regridded = interp_ap(xv, yv, xrvf, xrlat, xrlon)
    #print('data_regridded', data_regridded.shape)
    #print('XrVar',xrVar)
    Xrout = xr.DataArray(data_regridded, 
                        coords={'lat': lat_o,'lon': lon_o}, 
                        dims=["lat", "lon"])
    Xrout.attrs = xrVar.attrs
    #print('Xrout',Xrout)
    return Xrout
    
```

```python
xrlon = xr_getvar('lon',DS1).isel(time=0)
xrlat = xr_getvar('lat',DS1).isel(time=0)

if dinc is None: dinc = 1.  # increment of mesh in degrees
#dinc = 1.
lon_o = np.arange(0.,360+dinc,dinc)-180.
latn = np.floor((180.+dinc)/dinc).astype(int)
dincns = 180./latn
lat_o = np.linspace(-90.,90.,latn)

WRFnmdict = {'PS':'ps','TS':'ts','TGCLDLWP':'xxx'}
Varlist = np.array(['TS'])
for Varname in Varlist:
    # E3SM
    Var1 = get_NHDJF(Varname,DS1)    
    Var1ll = xr_intum2ll(Var1, lon_o, lat_o, xrLon=xrlon,xrLat=xrlat)
    fig, axes = plt.subplots(ncols=3,
                             gridspec_kw={'width_ratios': [1,1,1]},
                             subplot_kw={
                                 #'projection': ccrs.Orthographic(central_longitude=0.0, central_latitude=90.0)
                                 'projection': ccrs.NorthPolarStereo(central_longitude=0.0)
                                 },
                             figsize=(6,3),
                            )
    fig.set_dpi(300.0)
    #ax = axes[0]
    axes[0].set_extent([-180.,180., 60., 90], ccrs.PlateCarree())
    xr_llhplot(Var1ll,ax=axes[0])

    # WRF
    WRFname = WRFnmdict.get(Varname) 
    ind2 = fstring2 % (case_start2,WRFname,case_end2)
    print('opening ind2 ',ind2)
    
    DS2 = xr.open_mfdataset(ind2)
    time = DS2.time
    Var2 = DS2[WRFname].squeeze()
    Var2 = Var2.where(Var2.lat > 60.)
    Var2ll = xr_intum2ll(Var2, lon_o, lat_o)
    axes[1].set_extent([-180.,180., 60., 90], ccrs.PlateCarree())
    xr_llhplot(Var2ll,ax=axes[1])

    # Difference variable
    DV = Var1ll - Var2ll
    axes[2].set_extent([-180.,180., 60., 90], ccrs.PlateCarree())
    zmin = DV.min().values
    zmax = DV.max().values
    clevs = findNiceContours(np.array([zmin,zmax]),nlevs=10,rmClev=0.,sym=True)
    cmax = np.array([abs(zmin),abs(zmax)]).max()
    clevs = findNiceContlog(cmax,ndec=2)
    xr_llhplot(DV,ax=axes[2],clevs=clevs)

```
