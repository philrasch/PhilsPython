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
    decade = int(np.log10(cmax)+1)
    mylist = []
    for i in range(decade-ndec,decade):
        cdec = np.array([1.,2.,5.])*10.**i
        mylist.extend(cdec)
    ind = np.array(np.where(mylist > np.abs(cmax)))
    #print('cut1',mylist)
    #print('ind',ind,ind[0],len(ind[0]))
    if len(ind[0]) > 0: 
        ind = ind[0][0].item()
        mylist = mylist[0:ind]
    #print('cut2',mylist)

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
def fixWRF(Var):
    if Var.name == 'LWP':
        Var = Var*1000.
        Var.attrs['units'] = 'g m-2'
    if Var.name == 'IWP':
        Var = Var*1000.
        Var.attrs['units'] = 'g m-2'
    return Var
    
xrlon = xr_getvar('lon',DS1).isel(time=0)
xrlat = xr_getvar('lat',DS1).isel(time=0)

dmap = diverge_map()

if dinc is None: dinc = 1.  # increment of mesh in degrees
#dinc = 1.
lon_o = np.arange(0.,360+dinc,dinc)-180.
latn = np.floor((180.+dinc)/dinc).astype(int)
dincns = 180./latn
lat_o = np.linspace(-90.,90.,latn)

WRFnmdict = {'PS':'ps','TS':'ts','TGCLDLWP':'LWP',
             'TGCLDIWP':'IWP',
#             'FLDS':'rlds','FLUS':'rlus',
             'FLDS':'LW_d','FLUS':'LW_u',
             'LHFLX':'LH','SHFLX':'SH'}
Varlist = np.array(['PS','TS','TGCLDLWP','TGCLDIWP','FLDS','FLUS','LHFLX','SHFLX'])
# Varlist = np.array(['TGCLDIWP'])

for Varname in Varlist:
    # E3SM
    Var1 = get_NHDJF(Varname,DS1)
    print('range of Var1',Var1.min().values, Var1.max().values)
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

    # WRF
    WRFname = WRFnmdict.get(Varname) 
    ind2 = fstring2 % (case_start2,WRFname,case_end2)
    print('opening ind2 ',ind2)
    
    DS2 = xr.open_mfdataset(ind2)
    time = DS2.time
    Var2 = DS2[WRFname].squeeze()
    Var2 = fixWRF(Var2)
    Var2 = Var2.where(Var2.lat > 60.)
    print('range of Var2', Var2.min().values, Var2.max().values)
    Var2ll = xr_intum2ll(Var2, lon_o, lat_o)
    axes[1].set_extent([-180.,180., 60., 90], ccrs.PlateCarree())

    zmin = np.min([Var1.min().values,Var2.min().values])
    zmax = np.max([Var1.max().values,Var2.max().values])
    clevs = findNiceContours(np.array([zmin,zmax]),nlevs=12)
    print('clevs',clevs)

    xr_llhplot(Var1ll,ax=axes[0],clevs=clevs,cax=False,title='E3SM')
    pl = xr_llhplot(Var2ll,ax=axes[1],clevs=clevs,cax=False,title='RASM')

    # Difference variable
    DV = Var1ll - Var2ll
    axes[2].set_extent([-180.,180., 60., 90], ccrs.PlateCarree())
    zmin = DV.min().values
    zmax = DV.max().values
    print('DV range', zmin, zmax)
    #dlevs = findNiceContours(np.array([zmin,zmax]),nlevs=10,rmClev=0.,sym=True)
    dmax = np.array([abs(zmin),abs(zmax)]).max()
    dlevs = findNiceContlog(dmax,ndec=2)
    print('dlevs', dlevs)
    pl2 = xr_llhplot(DV,ax=axes[2],clevs=dlevs,cax=False,cmap=dmap,title='E3SM-RASM')

    posn = axes[1].get_position()
    # create an colorbar axis
    cax = fig.add_axes([0.,0.,0.8,0.1])
    cax2 = fig.add_axes([0.,0.,0.8,0.1])

    ## Adjust the positioning and orientation of the colorbar
    xoff = 0.5
    yoff = 0.06
    yoff2 = 0.18
    cax.set_position([posn.x0-0.5*xoff, posn.y0-yoff, posn.width+xoff, 0.03])
    cax2.set_position([posn.x0-0.5*xoff, posn.y0-yoff2, posn.width+xoff, 0.03])
    #cax.set_position([posn.x0, posn.y0-0.02, posn.width, 0.015])
    cbartitle = Var1ll.long_name
    cb = plt.colorbar(
         pl, orientation='horizontal',ticks=clevs,ax=cax,cax=cax,cmap=dmap,
         #label='%s (%s)'%(cbartitle, Var1ll.units)
         )
    cb.ax.tick_params(labelsize=7)
    cb2 = plt.colorbar(
         pl2, orientation='horizontal',ticks=dlevs,ax=cax2,cax=cax2,cmap=dmap,
         label='%s (%s)'%(cbartitle, Var1ll.units)
         )
    cb2.ax.tick_params(labelsize=7)
    
    plt.show()
```
