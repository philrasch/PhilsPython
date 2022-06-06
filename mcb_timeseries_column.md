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
import sys
print(sys.version)
%matplotlib inline
%run -i ~/Python/pjr3
```

```python
xc,yc = great_circle(190.,0,230.,35,n=40)
print('xc,yc',xc,yc)
```

```python
# process file holding data only over a region at every timestep of a model run
indir = '/global/cscratch1/sd/pjr/E3SMv2/v2.LR.histAMIP_e1/tests/M_1x2_ndays/run/v2.LR.histAMIP_e1.eam.h1.2015-07-*.nc'
regtag = '_190e_to_250e_0n_to_35n'
xr.set_options(keep_attrs=True)
DS0 = xr.open_mfdataset(indir)
DS = DS0
# reorder coords so ncol is alway first dim 
# to allow lat/lon interpolation across multiple dimensions
DS = DS0.transpose('ncol'+regtag,...) 

weights = np.cos(DS['lat'+regtag]/57.296) # approximate area weight
weights.attrs['name'] = 'weights'
weights.attrs['units'] = '1'

lon = DS['lon'+regtag].isel(time=0)#.squeeze()
print('lon',lon.shape,lon.min().values,lon.max().values)
lat = DS['lat'+regtag].isel(time=0)
print('lat',lat.shape,lat.min().values,lat.max().values)

Varname = 'ZCONVT'
Var = xr_getvar(Varname,DS,regtag)

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

lat = DS['lat'+regtag].isel(time=0).values
lon = DS['lon'+regtag].isel(time=0).values
print('shape of lat and lon arrays',lat.shape)
# calculate the GCD between tlon,tlat, and all the lons and lats from the model grid
tlon = 220
tlat = 25.
tlon = 222.
tlat = 25.
dd = gcd(lon,lat, tlon, tlat) # great circle dist from tlon, tlat

# identify the index for the nearest lat and lon and print it
ind = np.where(dd == np.min(dd))
ind = int(ind[0])
print('ind',ind,lat[ind],lon[ind])
```

```python


Tin = xr_getvar('CLDLIQ',DS,regtag=regtag)
PS = xr_getvar('PS',DS,regtag=regtag)
Pin = xr_getvar('P3',DS,regtag=regtag)

PCONVT = xr_getvar('PCONVT',DS,regtag=regtag)
PRECT = xr_getvar('PRECT',DS,regtag=regtag)
PRECC = xr_getvar('PRECC',DS,regtag=regtag)

Z3 = DS['Z3'+regtag]

Tinc = Tin.isel(ncol_190e_to_250e_0n_to_35n=ind)
Pinc = Pin.isel(ncol_190e_to_250e_0n_to_35n=ind)
Z3c = Z3.isel(ncol_190e_to_250e_0n_to_35n=ind)
PCONVTc = PCONVT.isel(ncol_190e_to_250e_0n_to_35n=ind)
PSc = PS.isel(ncol_190e_to_250e_0n_to_35n=ind)

print('pconvt, psc ', PCONVTc, PSc)
print('pconvt.min', PCONVTc.min().values)

tind = np.where(PCONVTc.values == np.min(PCONVTc.values))
tind = int(tind[0])
print('tind of pconvt.min is ', tind)
print('confirming', PCONVTc.isel(time=tind).values)
Pscx = Pinc.isel(time=tind)
PCONVTcx = PCONVTc.isel(time=tind)
print('Pscx',Pscx)
pd = np.argmin(np.abs(Pscx.values-PCONVTcx.values))
print('pd',pd)
#print('Pinc', Pinc.isel(time=tind).values)
Z3cx = Z3c.isel(time=tind)
print('Z3cx', Z3cx.values)
print('Pscxpd,Z3cxd', Pscx[pd].values, Z3cx[pd].values)
print('aaa',xr_getvar('ZCONVT',DS,regtag=regtag).isel(ncol_190e_to_250e_0n_to_35n=ind,time=tind).values)
print('PRECc, PREct', PRECC.isel(ncol_190e_to_250e_0n_to_35n=ind,time=tind).values,
      PRECT.isel(ncol_190e_to_250e_0n_to_35n=ind,time=tind).values)
PRECT.isel(ncol_190e_to_250e_0n_to_35n=ind).plot()
PRECC.isel(ncol_190e_to_250e_0n_to_35n=ind).plot()
plt.show()
PCONVTc.plot()
plt.show()

Varlist = ["ICWNC","ICWMR","CLOUDFRAC_CLUBB","PS","CAPE","T","CLOUD","CONCLD","ZCONVT"]
#Varlist = ["ZMDLF"]
for Vname in Varlist:
    Var = xr_getvar(Vname,DS,regtag=regtag)
    Varc = Var.isel(ncol_190e_to_250e_0n_to_35n=ind)
    print(Vname+' range', Varc.min().values, Varc.max().values)
    Varc.plot()
    plt.show()
```

```python
# next cells not used
help(xr_getvar)
1./0.
```

```python
Varlist = ["ICWNC","CLOUD","CLOUDFRAC_CLUBB"]
#Varlist = ["ZMDLF"]
for Vname in Varlist:
    Var = xr_getvar(Vname,DS,regtag=regtag)
    Varc = Var.isel(ncol_190e_to_250e_0n_to_35n=ind)
    print(Vname+' range', Varc.min().values, Varc.max().values)
#    Varc.plot()
#    plt.show()
    if Vname == "ICWNC":
        ICWNCc = Varc
    if Vname == "CLOUD":
        CLOUDc = Varc
    if Vname == "CLOUDFRAC_CLUBB":
        CFCc = Varc


af = ICWNCc.values.flatten()
afd = af[np.isfinite(af)]
print("ICWNCc.shape",ICWNCc.shape)
print('afd max',afd.max())

i,j = np.where(ICWNCc.values == afd.max())
print('i,j',i,j)
print('CLOUD, ICWNC max, CLUBB_CF', CLOUDc.values[i,j],  ICWNCc.values[i,j], CFCc.values[i,j])
print(CFCc)
(CLOUDc-CFCc).sel(lev=slice(500.,1000.)).plot()
plt.show()
```

```python

# specify the pressures to interpolate to
nzout = 40
pout = np.linspace(1.,1000.,nzout)


Tout = hy2plev(Tin, Pin, pout)
#Toutc = Tout.isel(time=0,ncol_190e_to_250e_0n_to_35n=ind).values

print("Tout",Tout)

```
