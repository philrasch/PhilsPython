---
jupyter:
  jupytext:
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.14.5
  kernelspec:
    display_name: Python [conda env:pjrpy3] *
    language: python
    name: conda-env-pjrpy3-py
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
    fig, axes = plt.subplots(ncols=1,
                             gridspec_kw={'width_ratios': [1]},
                             subplot_kw={'projection': ccrs.Mollweide()},
                             figsize=(6,3),
                            )
    fig.set_dpi(300.0)
    return fig, axes;
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
# plot a column specified by xlook and ylook for fields specified on a cubed sphere grid
xlook = 0.
ylook = -90.
lat = DS['lat'].values#[0,:]
lon = DS['lon'].values#[0,:]
#dist = np.abs(ylook-lat) + np.abs(xlook-lon)
dist = gcd(lon,lat, xlook, ylook)

print('dist',lon.min(),lat.min(),dist.min())
print('xxx',-88.939-ylook)
# find the index closest to xlook ylook
ind = np.where(dist == np.min(dist))
ind = int(ind[0])
print('ind',ind,lat[ind],lon[ind])
Toutm = Tout[:,ind]
pl = Toutm.plot(y="plev",yincrease=False)
plt.savefig('test.pdf',format='pdf',dpi=300)
```

```python
infile = "test.pdf"
add_prov(infile)
```

```python
# a simple method to plot a horizontal field using a map projection.
# this method is fine for plots where you dont care too much about the exact shape of the cells
# the plot uses triangulation from cell center to cell center. Other methods are better for knowing exactly the shape of each cell

#indir = os.path.expanduser('/lustre/choi040/20210920.F2010.1Nudg.ne30pg2_r05_oECv3/run/20210920.F2010.1Nudg.ne30pg2_r05_oECv3.eam.h0.2015-01.nc')
DS = xr.open_mfdataset(indir).chunk({'time': 20})
Tout = DS.PS.isel(time=0)
#indir = os.path.expanduser('/lustre/choi040/20210920.F2010.1Nudg.ne30pg2_r05_oECv3/run/20210920.F2010.1Nudg.ne30pg2_r05_oECv3.eam.h2.2015-01-01-00000.nc')
#DS = xr.open_mfdataset(indir).chunk({'time': 20})
#Tout = DS.T850.isel(time=0)

#gridfile = '/lustre/d3x345/maps/ne30pg2_scrip_c20191218.nc'

Tin = DS.T
Tin = Tin.squeeze()

Pin = (DS.hyam*DS.P0 + DS.hybm*DS.PS)/100.
Pin.attrs["units"] = 'hPa'
Pin.attrs["long_name"] = 'Pressure'
Pin = Pin.squeeze()

pout = 850.

Tout2 = hy2plev(Tin, Pin, pout).squeeze()

# Read data
data = Tout2

lon = DS['lon']
lat = DS['lat']

dataproj=ccrs.PlateCarree()    # data is always assumed to be lat/lon
plotproj=ccrs.Orthographic(central_latitude=-50)   # any projections should work 
ax = plt.axes(projection=plotproj)
ax.set_global()
ax.coastlines(linewidth=0.2)

tcoords = plotproj.transform_points(dataproj,np.array(lon[:]),np.array(lat[:]))
data2d = data
#xi=tcoords[:,0]!=np.inf
xi = np.where(~(np.isnan(tcoords[:,0])|np.isinf(tcoords[:,0])))[0] # this works for either
tc=tcoords[xi,:]
datai=data2d[:][xi]  # convert to numpy array, then subset
dmin = datai.min().values
dmax = datai.max().values
print('dmin,dmax',dmin,dmax)
pl = ax.tripcolor(tc[:,0],tc[:,1], datai,shading='gouraud',vmin=dmin,vmax=dmax) # looks good
#pl = ax.tripcolor(tc[:,0],tc[:,1], datai,shading='flat') # looks bad
# Add colorbar to plot
cb = plt.colorbar(
    pl, orientation='horizontal',
    label='%s (%s)'%(data.long_name, data.units), pad=0.05
)
plt.show
```

```python
# create a regular lat/lon grid
# OK for quick and dirty plotting, but leaves NaNs near domain edges
from scipy.interpolate import griddata
import numpy

xi = numpy.linspace(0, 360, 361)  # to regrid to 1/2 degree
yi = numpy.linspace(-90, 90, 181)  # to regrid to 1/2 degree

data_regridded = griddata((lon, lat), data2d, (xi[None,:], yi[:,None]), method='linear')
df = data_regridded.flatten()
dsub = df[np.isfinite(df)] # ignore NaN
zmax = dsub.max()
zmin = dsub.min()

dataproj=ccrs.PlateCarree()    # data is always assumed to be lat/lon
plotproj=ccrs.Orthographic(central_latitude=-50)   # any projections should work 
ax = plt.axes(projection=plotproj)
ax.set_global()
ax.coastlines(linewidth=0.2)
clevs = findNiceContours(np.array([dmin,dmax]),nlevs=10)
pl = ax.contourf(xi, yi, data_regridded, clevs, vmin=dmin, vmax=dmax,  extend='both', transform=ccrs.PlateCarree())
# Add colorbar to plot
cb = plt.colorbar(
    pl, orientation='horizontal',ticks=clevs,
    label='%s (%s)'%(data.long_name, data.units), pad=0.05
)
plt.show
```

```python

```
