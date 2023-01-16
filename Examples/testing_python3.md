---
jupyter:
  jupytext:
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.14.0
  kernelspec:
    display_name: Python [conda env:.conda-pjrpy3]
    language: python
    name: conda-env-.conda-pjrpy3-py
---

```python
import sys
print(sys.version)
%matplotlib inline
%run -i ~/Python/pjr3
from jupytext.config import find_jupytext_configuration_file
print('jupytext config file is ',find_jupytext_configuration_file('.'))
```

```python
print ("two vals of pi", np.pi, pi)
```

$$
f(t) = \Sigma_i (a_i cos({2\pi\over i} t) + b_i sin({2\pi\over i} t)))
$$


```python
hy2plev?
```

```python
plotZMf?
```

```python
from nco import Nco
```

```python
# identify a model case directory, and a directory to store remapped climo files
import os
#host = os.environ.get('HOST')
#print(host)
#print(os.environ)
import platform
host = platform.node()
print(host)
filename = os.path.expanduser('~/my_folder/output.txt')
if ('cori' in host):
    indir = '/global/cscratch1/sd/ogaruba/acme_scratch/cori-haswell/archive/E1850C5CLM45CNMC.ne30_oECv3_3/atm/hist/E1850C5CLM45CNMC.ne30_oECv3_3.cam.h0.*-06.nc'
    indir = '/global/cscratch1/sd/ogaruba/acme_scratch/cori-haswell/archive/E1850C5CLM45CNMC.ne30_oECv3_3/atm/hist/E1850C5CLM45CNMC.ne30_oECv3_3.cam.h0.*.nc'
    indir = '/global/cscratch1/sd/ogaruba/acme_scratch/cori-haswell/archive/E1850C5CLM45CNMC.ne30_oECv3_3/atm/hist/E1850C5CLM45CNMC.ne30_oECv3_3.cam.h0.0049-06.nc'
    indir = os.path.expanduser('~/NetCDF_Files/vd05_ANN_climo.nc')
else:
    indir = os.path.expanduser('~/NetCDF_Files/vd05_JJA_climo.nc')
    indir = os.path.expanduser('~/NetCDF_Files/F2010_PJR1.eam.h0.0001-01.nc')
    indir = os.path.expanduser('~/NetCDF_Files/F2010_PJR1.eam.h0.0001-01_fv192x288.nc')
print(indir)
#print('exists',os.path.exists(indir))
DS = xr.open_mfdataset(indir).chunk({'time': 20})
#print(DS)
weights = DS.area
weights.name = 'weights'
print(weights)
print('weights.sum',weights.sum().values,4.*pi)
#Var = DS.FSNT.isel(time=0)
Var = DS.FSNT - DS.FLNT
print(Var.shape)
#Varwt = Var.weighted(weights)
#print(Varwt)
#Varmean = Varwt.mean('ncol')
#print(Varmean)
#Varm2 = Var.weighted(weights).mean('ncol')
Varm2 = Var.weighted(weights).mean()
print('area weighted mean', Varm2.values)
```

```python
print(Varm2.values)
#Varm2.plot()
```

```python
FSNT = DS.FSNT
#FSNTg = FSNT.weighted(weights).mean('ncol')
FSNTg = FSNT.weighted(weights).mean()
FSNTg.plot()
```

```python

inCmd='ncdump -v time '+indir+' | grep "FRAC.*="'
outCmd = os.popen(inCmd).read()
print(inCmd)
print(outCmd)
```

```python
T = DS.T
print(T)
```

```python
#DS = xr.open_dataset('~/NetCDF_Files/vd05_ANN_climo.nc')
#print (DS.T) 
T = DS.T.isel(time=0)
print('T',T)
#T?
TZ = T.mean(dim='lon')
#TZ?
lev = TZ['lev']
lat = TZ['lat']
vals = TZ.values
print('plotting on eta levels')
plotZMf(vals, lat, lev)

```

```python
# demonstrate shape specification to force correct broadcasting
# add a vector of numbers to a particular axis (in this case axis=1)
pout = np.arange(20)
x = np.zeros([2,20,100])
newshape = [1,20,1]
z = x + pout.reshape(newshape)
print('z shape', z.shape)
print(z[0,:,3])
print(z[1,:,-1])
```

```python
# interpolate xarray datavariable on eta hybrid surfaces to pressure

indir = os.path.expanduser('~/NetCDF_Files/F2010*-01.nc')
#indir = os.path.expanduser('/lustre/choi040/20210920.F2010.1Nudg.ne30pg2_r05_oECv3/run/20210920.F2010.1Nudg.ne30pg2_r05_oECv3.eam.h2.2015-01-01-00000.nc')
#indir = os.path.expanduser('~/NetCDF_Files/*F2010*01.nc')
#indir = os.path.expanduser('~/NetCDF_Files/vd05_ANN_climo.nc')

print(indir)
#print('exists',os.path.exists(indir))
DS = xr.open_mfdataset(indir).chunk({'time': 20})


    
Tin = DS.T
Tin = Tin.squeeze()
#print("Tin",Tin)
Pin = (DS.hyam*DS.P0 + DS.hybm*DS.PS)/100.
Pin.attrs["units"] = 'hPa'
Pin.attrs["long_name"] = 'Pressure'
Pin = Pin.squeeze()
nzout = 20
pout = np.linspace(1.,1000.,nzout)
#pout = [850.]
#pout = 850.
#print ("pout", pout.shape, pout)

Tout = hy2plev(Tin, Pin, pout)

```

```python
%run -i ~/Python/pjr3
from cartopy import crs

indir = os.path.expanduser('/lustre/choi040/20210920.F2010.1Nudg.ne30pg2_r05_oECv3/run/20210920.F2010.1Nudg.ne30pg2_r05_oECv3.eam.h0.2015-01.nc')
indir = os.path.expanduser('~/NetCDF_Files/F2010_PJR1.eam.h0.0001-01.nc')

DS = xr.open_mfdataset(indir).chunk({'time': 20})
Tout = DS.PS.isel(time=0)
#indir = os.path.expanduser('/lustre/choi040/20210920.F2010.1Nudg.ne30pg2_r05_oECv3/run/20210920.F2010.1Nudg.ne30pg2_r05_oECv3.eam.h2.2015-01-01-00000.nc')
#DS = xr.open_mfdataset(indir).chunk({'time': 20})
#Tout = DS.T850.isel(time=0)

Tin = DS.T
Tin = Tin.squeeze()
#print("Tin",Tin)
Pin = (DS.hyam*DS.P0 + DS.hybm*DS.PS)/100.
Pin.attrs["units"] = 'hPa'
Pin.attrs["long_name"] = 'Pressure'
Pin = Pin.squeeze()
pout = 850.

Tout2 = hy2plev(Tin, Pin, pout).squeeze()

# Read data
data = Tout2
print('plotting ', data)

lon = DS['lon']
lat = DS['lat']

dataproj=crs.PlateCarree()    # data is always assumed to be lat/lon
plotproj=crs.Orthographic(central_latitude=-50)   # any projections should work 
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
# plot a column specified by xlook and ylook for fields specified on a cubed sphere grid
xlook = 0.
ylook = -90.
lat = DS['lat'].values#[0,:]
#print('latshape',lat.shape)
#print('lat4',lat[0:3])
lon = DS['lon'].values#[0,:]
#print('lon4',lon[0:3])
#print('lat',lat.min(),lat.max())
#print('lon',lon.min(),lon.max())
dist = np.abs(ylook-lat) + np.abs(xlook-lon)
#print('dist',dist)
ind = np.where(dist == np.min(dist))
ind = int(ind[0])
#print('ind',ind, ind.shape, ind[0])
print('ind',ind,lat[ind],lon[ind])
print('Tout', Tout.shape, Tout.squeeze().shape)
#Toutm = Tout.mean(dim=['ncol','time'])
#Tinm = Tin[0,:,ind]
#print('Tinm',Tinm.values)
#Pinm = Pin[0,:,ind]
#print('Pinm',Pinm.values)
#tind = 0
#for index, item in enumerate(Pinm):
#    print('Pinx ',index, Pinm[index].values, Tinm[index].values)
#Toutm = Tout[:,ind]
Toutm = Tout[ind]
#print('Toutm',Toutm)
#for index, item in enumerate(Toutm):
#    print('Poutx ',index, pout[index], Toutm[index].values)

Toutm.plot()
```

```python
# fast 1d interpolation of fields along an axis
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d as scipy1d

# toy coordinates and data
nx, ny, nz = 25, 30, 10
x = np.arange(nx)
y = np.arange(ny)
z = np.tile(np.arange(nz), (nx,ny,1)) + np.random.randn(nx, ny, nz)*.1
testdata = np.random.randn(nx,ny,nz) # x,y,z

# Desired z-coordinates (must be between bounds of z)
znew = np.tile(np.linspace(2,nz-2,50), (nx,ny,1)) + np.random.randn(nx, ny, 50)*0.01

# Inverse the coordinates for testing
z = z[..., ::-1]
znew = znew[..., ::-1]

# Now use own routine 
ynew = interp_along_axis(testdata, z, znew, axis=2, inverse=True, method='cubic')

# Check some random profiles
for i in range(5):
    randx = np.random.randint(nx)
    randy = np.random.randint(ny)

    checkfunc = scipy1d(z[randx, randy], testdata[randx,randy], kind='cubic')
    checkdata = checkfunc(znew)

    fig, ax = plt.subplots()
    ax.plot(testdata[randx, randy], z[randx, randy], 'x', label='original data')
    ax.plot(checkdata[randx, randy], znew[randx, randy], label='scipy')
    ax.plot(ynew[randx, randy], znew[randx, randy], '--', label='Peter')
    ax.legend()
    plt.show()

```

```python

```
