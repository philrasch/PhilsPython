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

print(indir)
#print('exists',os.path.exists(indir))
DS = xr.open_mfdataset(indir).chunk({'time': 20})
print('DS',DS)
weights = DS.area
weights.name = 'weights'
Var = DS.FSNT - DS.FLNT  # net radiative flux at Top of Model
print(Var.shape)
Varm2 = Var.weighted(weights).mean() # calculate global average
print('area weighted mean', Varm2.values)
```

```python
# interpolate xarray datavariable on eta hybrid surfaces to 20 pressure levels
Tin = DS.T
Tin = Tin.squeeze()
#print("Tin",Tin)

# calculate the pressure levels for the model from PS and the hybrid coordinates
Pin = (DS.hyam*DS.P0 + DS.hybm*DS.PS)/100.
Pin.attrs["units"] = 'hPa'
Pin.attrs["long_name"] = 'Pressure'
Pin = Pin.squeeze()

# specify the pressures to interpolate to
nzout = 20
pout = np.linspace(1.,1000.,nzout)


Tout = hy2plev(Tin, Pin, pout)
print('Tout',Tout)
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

data2d = DS['PS'].squeeze()
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
plotproj=ccrs.Orthographic(central_latitude=0,central_longitude=55)   # any projections should work 
ax = plt.axes(projection=plotproj)
ax.set_global()
ax.coastlines(linewidth=0.2)
clevs = findNiceContours(np.array([zmin,zmax]),nlevs=10)
pl = ax.contourf(xi, yi, data_regridded, clevs, vmin=zmin, vmax=zmax,  extend='both', transform=ccrs.PlateCarree())
# Add colorbar to plot
cb = plt.colorbar(
    pl, orientation='horizontal',ticks=clevs,
    label='%s (%s)'%(data2d.long_name, data2d.units), pad=0.05
)
plt.plot(xc,yc,marker='*',color='green',transform=ccrs.PlateCarree())

plt.show
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
