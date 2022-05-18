---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.13.8
  kernelspec:
    display_name: pjrpy3
    language: python
    name: pjrpy3
---

```python
# paired with markdown file
# added a comment in markdown version and annotated it
import sys
print(sys.version)
%matplotlib inline
%run -i ~/Python/pjr3
```

```python
indir = os.path.expanduser('~/NetCDF_Files/*F2010*.nc')
#indir = os.path.expanduser('/lustre/choi040/20210920.F2010.1Nudg.ne30pg2_r05_oECv3/run/20210920.F2010.1Nudg.ne30pg2_r05_oECv3.eam.h2.2015-01-01-00000.nc')
indir = os.path.expanduser('~/NetCDF_Files/F2010*01.nc')
#indir = os.path.expanduser('~/NetCDF_Files/vd05_ANN_climo.nc')
indir = os.path.expanduser('/global/cscratch1/sd/pjr/E3SMv2/v2.LR.histAMIP_e1/tests/M_1x1_nmonths/run/v2.LR.histAMIP_e1.eam.h4.2015-07-01-00000.nc')
print(indir)
#print('exists',os.path.exists(indir))
DS = xr.open_mfdataset(indir)#.chunk({'time': 20})
xr.set_options(keep_attrs=True)
#print('DS',DS)
weights = DS.area
weights.name = 'weights'
Var = DS.PRECT
Var = Var*8.64e7
Var.attrs['units'] = 'mm/day'

Var = DS.PCONVT
Var = Var/100.
Var.attrs['units'] = 'hPa'
print('Var range ',Var.min().values, Var.max().values)
Var.attrs['units'] = 'm'
Var = -7.1e3*np.log(Var/1000.)  # rough conversion to m using 7km scale height
Var.attrs['long_name'] = 'convection top height'

#Var = DS.PRECC/DS.PRECT
#Var.attrs['long_name'] = 'frac of precip that is convective'
#Var.attrs['units'] = 1.
#Var.attrs['name'] = 'FRACC'

Var = Var.mean(dim='time')

print('Var',Var)
Varm2 = Var.weighted(weights).mean() # calculate global average
#zavg = -7.1e3*np.log(Varm2.values/1000.)
#print('area weighted mean', Varm2.values, zavg)
```

```python
def great_circle(lon1, lat1, lon2, lat2,n=None):
    """see http://www.edwilliams.org/avform147.htm#Intermediate"""
    from math import radians, degrees, sin, cos, asin, acos, sqrt, atan2

    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
    #print('xxx',lon1,lat1,lon2,lat2)
    d = (acos(sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(lon1 - lon2)))
    #print('d',d)
    f = 0.5
    if n == None:
        f = np.array(np.linspace(0.,1.,3))

    else:
        f = np.array(np.linspace(0.,1.,n))
        
    #print ('f',f.dtype,f.shape)
    a = np.sin((1-f)*d)/np.sin(d)
    b = np.sin(f*d)/np.sin(d)
    clat1 = np.cos(lat1)
    clon1 = np.cos(lon1)
    clat2 = np.cos(lat2)
    clon2 = np.cos(lon2)
    slon1 = np.sin(lon1)
    slon2 = np.sin(lon2)
    slat1 = np.sin(lat1)
    slat2 = np.sin(lat2)
    x = a*clat1*clon1 + b*clat2*clon2
    y = a*clat1*slon1 + b*clat2*slon2
    z = a*slat1       + b*slat2
    lat = np.arctan2(z,np.sqrt(x**2+y**2))
    lon = np.arctan2(y,x)
    lon,lat = np.degrees([lon,lat])

    rad = 6371. 
    return lon,lat


#lon,lat = great_circle(-0.08,51.53,132.,43.17,n=10)
xc,yc = great_circle(0.,50,110.,-40,n=10)
xc,yc = great_circle(-30.,50,110.,-40,n=10)
xc,yc = great_circle(190.,0,230.,35,n=10)


#yc = -yc
#xc =np.concatenate((xc,xc))
#yc =np.concatenate((yc,-yc))

print('xc,yc',xc,yc)

dataproj=ccrs.PlateCarree()    # data is always assumed to be lat/lon
plotproj=ccrs.Orthographic(central_latitude=-50)   # any projections should work 
plotproj=ccrs.Mollweide(central_longitude=55) 
ax = plt.axes(projection=plotproj)
ax.set_global()
ax.coastlines(linewidth=0.2)
plt.plot(xc,yc,marker='*',color='green',transform=ccrs.PlateCarree())

plt.show
```

```python
# create a regular lat/lon grid
from scipy.interpolate import griddata
import numpy

xi = numpy.linspace(0, 360, 361)  # to regrid to 1/2 degree
yi = numpy.linspace(-90, 90, 181)  # to regrid to 1/2 degree

data2d = Var
print(data2d)
print('d',data2d.shape)
lon = DS['lon']
print('lon',lon.shape)
lat = DS['lat']
print('lat',lat.shape)

#data_regridded = griddata((lon, lat), data2d, (xi[None,:], yi[:,None]), method='linear')
data_regridded = interp_to_latlon(data2d.values,lat.values,lon.values,yi,xi)

df = data_regridded.flatten()
dsub = df[np.isfinite(df)] # ignore NaN
zmax = dsub.max()
zmin = dsub.min()

dataproj=ccrs.PlateCarree()    # data is always assumed to be lat/lon
#plotproj=ccrs.Orthographic(central_latitude=0,central_longitude=55)   # any projections should work 
plotproj=ccrs.Robinson(central_longitude=240)   # any projections should work 
ax = plt.axes(projection=plotproj)
ax.set_global()
ax.coastlines(linewidth=0.2,color='white')
clevs = findNiceContours(np.array([zmin,zmax]),nlevs=10)
pl = ax.contourf(xi, yi, data_regridded, clevs, vmin=zmin, vmax=zmax,  extend='both', transform=ccrs.PlateCarree())
print('xxx',data2d.long_name, data2d.units)
# Add colorbar to plot
cb = plt.colorbar(
    pl, orientation='horizontal',ticks=clevs,
    label='%s (%s)'%(data2d.long_name, data2d.units), pad=0.05
)
plt.plot(xc,yc,marker='*',color='green',transform=ccrs.PlateCarree())

plt.show
```

```python
inds = np.where(
    #((lat.values < 35.)&(lat.values > 0.)&(lon.values > 190)&(lon.values < 230))
    ((lat <= 35.)&(lat >= 0.)&(lon >= 190)&(lon <= 240))
    )[0]
Varsub = Var[inds]
print('Varsub',np.shape(Varsub))
print('eee',Varsub)
latsub = lat[inds]
lonsub = lon[inds]
dinc = 1.  # increment of mesh in degrees
lon_h=np.arange(np.floor(lonsub.values.min()),np.ceil(lonsub.values.max()+dinc), dinc)
lat_h=np.arange(np.floor(latsub.values.min()),np.ceil(latsub.values.max()+dinc), dinc)
print('ddd',latsub.values.min(),latsub.values.max())
print('fff',lonsub.values.min(),lonsub.values.max())

print('bbb',lon_h.min(),lon_h.max())
print('ccc',lat_h.min(),lat_h.max())
xv,yv=np.meshgrid(lon_h,lat_h)
print('aaa',xv.shape)

```

```python
#%%timeit
data_regridded = interp_ap(xv, yv, Varsub.values,latsub.values,lonsub.values)


```

```python
df = data_regridded.flatten()
dsub = df[np.isfinite(df)] # ignore NaN
zmax = dsub.max()
zmin = dsub.min()
#zmax = 1000.
#zmin = 700.

dataproj=ccrs.PlateCarree()    # data is always assumed to be lat/lon
#plotproj=ccrs.Orthographic(central_latitude=0,central_longitude=55)   # any projections should work 
plotproj=ccrs.Mollweide(central_longitude=200)   # any projections should work 
clat = (latsub.values.min()+latsub.values.max())/2.
clon = (lonsub.values.min()+lonsub.values.max())/2.
plotproj=ccrs.NearsidePerspective(central_longitude=clon, central_latitude=clat)
ax = plt.axes(projection=plotproj)
ax.set_global()
ax.coastlines(linewidth=0.2)
clevs = findNiceContours(np.array([zmin,zmax]),nlevs=10)
pl = ax.contourf(xv, yv, data_regridded, clevs, vmin=zmin, vmax=zmax,  extend='both', transform=ccrs.PlateCarree())
# Add colorbar to plot
cb = plt.colorbar(
    pl, orientation='horizontal',ticks=clevs,
    label='%s (%s)'%(data2d.long_name, data2d.units), pad=0.05
)
plt.plot(xc,yc,marker='*',color='green',transform=ccrs.PlateCarree())

plt.show
```

```python
indir = '/global/cscratch1/sd/pjr/E3SMv2/v2.LR.histAMIP_e1/tests/M_1x2_ndays/run/v2.LR.histAMIP_e1.eam.h1.2015-07-*.nc'
regtag = '_190e_to_250e_0n_to_35n'
xr.set_options(keep_attrs=True)
DS0 = xr.open_mfdataset(indir)
DS = DS0
# reorder coords so ncol is alway first dim 
# to allow interpolation across multiple dimensions
DS = DS0.transpose('ncol'+regtag,...) 

#print('DS',DS)
print('Z3v1',DS['Z3'+regtag].shape)
weights = np.cos(DS['lat'+regtag]/57.296) # approximate area weight
weights.attrs['name'] = 'weights'
weights.attrs['units'] = '1'

lon = DS['lon'+regtag].isel(time=0)#.squeeze()
print('lon0',lon)
print('lon',lon.shape,lon.min().values,lon.max().values)
lat = DS['lat'+regtag].isel(time=0)
print('lat',lat.shape,lat.min().values,lat.max().values)

Varname = 'ZCONVT'
if Varname == "PRECT":
    Var = DS['PRECT'+regtag]
    Var = Var*8.64e7
    Var.attrs['units'] = 'mm/day'
elif Varname == 'PCONVT':
    Var = DS['PCONVT'+regtag]
    Var = Var/100.
    Var.attrs['units'] = 'hPa'
elif Varname == 'ZCONVT':
    Var = DS['PCONVT'+regtag]
    Var = Var/100.
    Var.attrs['units'] = 'm'
    Var = -7.1e3*np.log(Var/1000.)  # rough conversion to m using 7km scale height
    Var.attrs['long_name'] = 'convection top height'
    Var = Var.rename(Varname)
else:
    print('Varname not found: ', Varname)


print('Var',Var)
#Var = Var.mean(dim='time')
#Var = Var.isel(time=-1)
print('shape',Var.shape)

dinc = 1.  # increment of mesh in degrees
lon_h=np.arange(np.floor(lon.values.min()),np.ceil(lon.values.max()+dinc), dinc)
lat_h=np.arange(np.floor(lat.values.min()),np.ceil(lat.values.max()+dinc), dinc)

xv,yv=np.meshgrid(lon_h,lat_h)
#print('yyy',xv.shape)
#data_regridded = interp_ap(xv, yv, Var.values,xv,lon.values)
data_regridded = interp_ap(xv, yv, Var.values,lat.values,lon.values)
data_regridded = data_regridded[...,-1]
print('xxx',data_regridded.shape)

df = data_regridded.flatten()
dsub = df[np.isfinite(df)] # ignore NaN
zmax = dsub.max()
zmin = dsub.min()
print('zmin,zmax', zmin, zmax)
#zmax = 1000.
#zmin = 700.

dataproj=ccrs.PlateCarree()    # data is always assumed to be lat/lon
#plotproj=ccrs.Orthographic(central_latitude=0,central_longitude=55)   # any projections should work 
plotproj=ccrs.Mollweide(central_longitude=200)   # any projections should work 
clat = (lat.values.min()+lat.values.max())/2.
clon = (lon.values.min()+lon.values.max())/2.
plotproj=ccrs.NearsidePerspective(central_longitude=clon, central_latitude=clat)
ax = plt.axes(projection=plotproj)
ax.set_global()
ax.coastlines(linewidth=0.2)
clevs = findNiceContours(np.array([zmin,zmax]),nlevs=10)
#clevs = [0.,0.1,0.2,0.5,1.,2.,5.,10.,20.]
pl = ax.contourf(xv, yv, data_regridded, clevs, vmin=zmin, vmax=zmax,  extend='both', transform=ccrs.PlateCarree())
# Add colorbar to plot
cb = plt.colorbar(
    pl, orientation='horizontal',ticks=clevs,
    label='%s (%s)'%(Var.long_name, Var.units), pad=0.05
)
#plt.plot(xc,yc,marker='*',color='green',transform=ccrs.PlateCarree())

plt.show
      
```

```python
# interpolate xarray datavariable on eta hybrid surfaces to 20 pressure levels
Tin = DS['CLDLIQ'+regtag]
#print("Tin",Tin)
#print('Tin.shape',Tin.shape)
#print('Tin.dims',Tin.dims)

PS = DS['PS'+regtag]
Pin = (DS.hyam*DS.P0 + DS.hybm*PS)/100.
Pin.attrs["units"] = 'hPa'
Pin.attrs["long_name"] = 'Pressure'
# make sure the two quantities have the same coordinate order
ldims = list(Tin.dims)
Pin = Pin.transpose(*ldims)
#print('newPin.dims', Pin.dims)
#print('newPin.shape', Pin.shape)

#Tin = Tin.isel(time=0)
#Pin = Pin.isel(time=0)
# specify the pressures to interpolate to
nzout = 40
pout = np.linspace(1.,1000.,nzout)


Tout = hy2plev(Tin, Pin, pout)
#print('Tout',Tout)
Tinc = Tin.isel(time=0,ncol_190e_to_250e_0n_to_35n=100).values
Pinc = Pin.isel(time=0,ncol_190e_to_250e_0n_to_35n=100).values
Toutc = Tout.isel(time=0,ncol_190e_to_250e_0n_to_35n=100).values
plt.plot(Tinc,Pinc)
plt.plot(Toutc,pout)
plt.show()
```

```python
xc,yc = great_circle(190.,0,230.,35,n=10)
data_regridded = interp_ap(xc, yc, Tout,lat.values,lon.values)
#transect coord, time, pressure)
print('xxx',data_regridded.shape)
```

```python
def gcd(lon1, lat1, lon2, lat2):
    """calculate great circle distance in km
       given input longitude and latitudes in degrees
    """
    d2r = 57.296
    lonr1 = lon1/d2r
    latr1 = lat1/d2r
    lonr2 = lon2/d2r
    latr2 = lat2/d2r

    d = 6371.*(np.arccos(np.sin(latr1) * np.sin(latr2) + np.cos(latr1) * np.cos(latr2) * np.cos(lonr1 - lonr2)))
    return d

lat = DS['lat'].values#[0,:]
lon = DS['lon'].values#[0,:]
print('xxx',lat[10],lon[10])
# calculate the GCD between tlon,tlat, and all the lons and lats from the model grid
tlon = 321.6
tlat = -35.8
dd = gcd(lon,lat, tlon, tlat)

# identify the index for the nearest lat and lon and print it
ind = np.where(dd == np.min(dd))
ind = int(ind[0])
print('ind',ind,lat[ind],lon[ind])
```

```python
data2d = DS['PS'].squeeze()
print(data2d)
#print('d',data2d.shape)
lon = DS['lon']
#print('lon',lon.shape)
lat = DS['lat']
#print('lat',lat.shape)

# over-write the data with an analytic fields so we can check the interpolation at arbitrary points
data2d.values = (np.sin(3*np.radians(lon.values)-np.pi/4.))**2 + (np.sin(4*np.radians(lat.values)))**3
print('true', (np.sin(3*np.radians(xc)-np.pi/4.))**2 + (np.sin(4*np.radians(yc)))**3)

#data_regridded = griddata((lon, lat), data2d, (xi[None,:], yi[:,None]), method='linear')


```

```python
%%timeit
data_regridded = interp_ap(xc, yc, data2d.values,lat.values,lon.values)
print('vals', data_regridded)
```

```python
%%timeit
data_regridded2 = interp_ap(xc, yc, data2d.values,lat.values,lon.values,method='cubic')
print('vals', data_regridded2)
```

```python
xim,yim = np.meshgrid(xi, yi)
print('xim',xim.shape,xim.min(),xim.max())
truth = (np.sin(3*np.radians(xim)-np.pi/4.))**2 + (np.sin(4*np.radians(yim)))**3
```

```python
%%timeit
data_regridded = interp_ap(xim, yim, data2d.values,lat.values,lon.values)
err = (data_regridded - truth)[((yim > -45) & (yim < 45))]
print('crude RMS and 1 norm', np.sqrt((err**2).mean()),err.min(),err.max())
```

```python
%%timeit
data_regridded = interp_ap(xim, yim, data2d.values,lat.values,lon.values,method='cubic')
err = (data_regridded - truth)[((yim > -45) & (yim < 45))]
print('crude RMS and 1 norm', np.sqrt((err**2).mean()),err.min(),err.max())
```

```python

```
