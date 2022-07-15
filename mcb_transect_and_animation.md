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
indir = os.path.expanduser('/global/cscratch1/sd/pjr/E3SMv2/v2.LR.histAMIP_x3/tests/M_1x2_ndays/run/v2.LR.histAMIP_x3.eam.h1.2015-07-*.nc')
print(indir)
#print('exists',os.path.exists(indir))
DS = xr.open_mfdataset(indir)#.chunk({'time': 20})
xr.set_options(keep_attrs=True)
#print('DS',DS)
#weights = DS.area
regtag = '_190e_to_250e_0n_to_35n'
lat = xr_getvar('lat',DS,regtag=regtag).isel(time=0)
lon = xr_getvar('lon',DS,regtag=regtag).isel(time=0)
print('lat',lat.min().values,lat.max().values)
print('lon',lon.min().values,lon.max().values)
weights = np.cos(lat/57.296) # approximate area weight
weights.attrs['name'] = 'weights'
weights.attrs['units'] = '1'

Var = xr_getvar('ZCONVT',DS,regtag=regtag)


Var = Var.mean(dim='time')

print('Var',Var)
Varm2 = Var.weighted(weights).mean() # calculate global average
#zavg = -7.1e3*np.log(Varm2.values/1000.)
#print('area weighted mean', Varm2.values, zavg)
```

```python

dinc = 1.  # increment of mesh in degrees
lon_h=np.arange(np.floor(lon.values.min()),np.ceil(lon.values.max()+dinc), dinc)
lat_h=np.arange(np.floor(lat.values.min()),np.ceil(lat.values.max()+dinc), dinc)
xv,yv=np.meshgrid(lon_h,lat_h)

```

```python
#%%timeit
data_regridded = interp_ap(xv, yv, Var.values,lat.values,lon.values)
print(Var)

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
clat = (lat.values.min()+lat.values.max())/2.
clon = (lon.values.min()+lon.values.max())/2.
#plotproj=ccrs.NearsidePerspective(central_longitude=clon, central_latitude=clat)
plotproj=ccrs.Mercator()
ax = plt.axes(projection=plotproj)
ax.set_extent([lon.values.min(), 260., lat.values.min(), lat.values.max()])
#ax.set_global()
clevs = findNiceContours(np.array([zmin,zmax]),nlevs=10)
pl = ax.contourf(xv, yv, data_regridded, clevs, vmin=zmin, vmax=zmax,  extend='both', transform=ccrs.PlateCarree())
# Add colorbar to plot
cb = plt.colorbar(
    pl, orientation='horizontal',ticks=clevs,
    label='%s (%s)'%(Var.long_name, Var.units), pad=0.1
)
cb.ax.tick_params(labelsize=8)

xc,yc = great_circle(190.,0,230.,35,n=40)
print('xc,yc',xc,yc)
plt.plot(xc,yc,marker='*',color='red',transform=ccrs.PlateCarree())
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=2, color='gray', alpha=0.5)
ax.coastlines(linewidth=1,color='blue')

plt.savefig("transect.pdf")
plt.show
```

```python
# interpolate xarray datavariable on eta hybrid surfaces to 20 pressure levels
Tin = xr_getvar("T",DS,regtag=regtag)
Tin = xr_getvar("CLDLIQ",DS,regtag=regtag)
Pin = xr_getvar('P3',DS,regtag=regtag)
# make sure the two quantities have the same coordinate order
ldims = list(Tin.dims)
Pin = Pin.transpose(*ldims)
#print('newPin.dims', Pin.dims)
#print('newPin.shape', Pin.shape)

Var = xr_getvar('PCONVT',DS,regtag=regtag)

# specify the pressures to interpolate to
nzout = 40
pout = np.linspace(1.,1000.,nzout)


Tout = hy2plev(Tin, Pin, pout)
#print('Tout',Tout)
Tinc = Tin.isel(time=0,ncol_190e_to_250e_0n_to_35n=100).values
Pinc = Pin.isel(time=0,ncol_190e_to_250e_0n_to_35n=100).values
Toutc = Tout.isel(time=0,ncol_190e_to_250e_0n_to_35n=100).values

Varc = Var.isel(time=0,ncol_190e_to_250e_0n_to_35n=100).values
print('pconvt ', Varc)

plt.plot(Tinc,Pinc)
plt.plot(Toutc,pout)
plt.show()
```

```python
Tin = Tin.transpose('ncol'+regtag,...)
Tin.shape
Tout = Tout.transpose('ncol'+regtag,...)
```

```python
xc,yc = great_circle(190.,0,230.,35,n=44)
data_regridded = interp_ap(xc, yc, Tout,lat.values,lon.values)
#transect coord, time, pressure)
print('data ordered (cols, time, level)',data_regridded.shape)
```

```python
xc,yc = great_circle(190.,0,230.,35,n=40)
#print('xc',xc)
data2 = interp_ap(xc, yc, Tin.values,lat.values,lon.values)
#transect coord, time, pressure)
print('data2 ordered (cols, time, level)',data2.shape)
```

```python
# example of locating max of a multidimensional array
a = np.array([[1,2,3],[4,3,1]])
a = np.array([[1,4,3,np.nan],[np.nan,4,3,1]])
af = a.flatten()
afd = af[np.isfinite(af)]
print('a',a)
print('afd max',afd.max())
i,j = np.where(a==a.max())
print('i,j',i,j)
i,j = np.where(a==afd.max())
print('i,j',i,j)
a[i,j]
```

```python
import matplotlib.pyplot as plt
import matplotlib.animation
import numpy as np
import time
import datetime
import matplotlib.animation as animation

print('Tin.dims',Tin.dims)
lev = DS['lev']
#print('lev.values',lev.values)
timev = DS['time'].values
#print('time.values',time)
dtime = (timev - timev[0])/datetime.timedelta(hours=1) # convert time variable to hours as a float
levs = lev.values

l = len(timev)
print('l=',l)

data2_hztf = data2.flatten()
data2_hztfd = data2_hztf[np.isfinite(data2_hztf)]
zmin = data2_hztfd.min()
zmax = data2_hztfd.max()
print('zmin,zmax', zmin, zmax)
clevs = findNiceContours(np.array([zmin,zmax]),nlevs=10)
print('clevs',clevs)
i,j,k = np.where(data2 == zmax) # coord where data2 takes on max val
print('max at ', i, j, k)
print('data2 max', data2[i[0],j[0],k[0]])

xv,yv=np.meshgrid(xc,levs)
fig, axes = plt.subplots(nrows=2
                        ,gridspec_kw={'height_ratios': [1, 3]}
                        ,figsize=(7,7)
                        )

def getvals(tind):   
    #print('tind',tind,dtime[tind])
    #ax = plt.axes()
    
    data2_hz = data2[:,tind,:].transpose()

    PRECT = xr_getvar('PRECT',DS,regtag=regtag).isel(time=tind)
    ptvals = interp_ap(xc,yc,PRECT.values,lat.values,lon.values)

    title = 'hr='+"%05.2f" % dtime[tind]

    PCONVT = xr_getvar('PCONVT',DS,regtag=regtag).isel(time=tind)
    pctvals = interp_ap(xc, yc, PCONVT.values,lat.values,lon.values)

    PRECC = xr_getvar('PRECC',DS,regtag=regtag).isel(time=tind)
    pcvals = interp_ap(xc,yc,PRECC.values,lat.values,lon.values)
    
    return title, ptvals, pcvals, pctvals, data2_hz,

axes[1].set_ylim([1000., 100.])
axes[0].set_xlim([xc.min(), xc.max()])
axes[1].set_xlim([xc.min(), xc.max()])
axes[0].legend(['total','convective'])
axes[0].set_ylim([0.01, 20.])
axes[0].set_yscale('log')

props = dict(boxstyle='round', facecolor='wheat')
timelabel = axes[1].text(0.5,0.9, "", transform=ax.transAxes, ha="right", bbox=props)

# grab data  used for initial figure
title, v1,v2,v3,v4 = getvals(0) #time,prect, precv, pconvt
# create objects used to animate the data display
p = [axes[1].contourf(xv, yv, v4, clevs, vmin=zmin, vmax=zmax,  extend='both')]
lct, = axes[1].plot(xc,v3,linewidth=2,color="white")
rct, = axes[0].plot(xc,v1,linewidth=2,color="red",label='total')
bct, = axes[0].plot(xc,v2,linewidth=2,color="blue",label='convective')
timelabel.set_text(title)
axes[0].legend(['total','convective'],loc='upper right')


cb = fig.colorbar(
        p[0], orientation='horizontal',ticks=clevs,
        label='%s (%s)'%(Tin.long_name, Tin.units), pad=0.20, ax=axes[1]
    )


def update(i):
    for tp in p[0].collections:
        tp.remove()
    title, v1,v2,v3,v4 = getvals(i)
    p[0] = axes[1].contourf(xv, yv, v4, clevs, vmin=zmin, vmax=zmax,  extend='both')
    #replace line data with new values
    lct.set_data(xc, v3)
    rct.set_data(xc, v1)
    #print('v1 lim',i,v1[np.isfinite(v1)].min(),v1[np.isfinite(v1)].max())
    bct.set_data(xc, v2)
    timelabel.set_text(title)
    return p[0].collections+[timelabel]

#plt.show()
ani = matplotlib.animation.FuncAnimation(fig, update, frames=len(dtime), 
                                         interval=10, blit=True, repeat=True)
ani.save('t1.gif', writer = 'imagemagick')
```

```python
#debug cells below
1./0.
```

```python
import matplotlib.pyplot as plt
import matplotlib.animation
import numpy as np
import time

x= np.linspace(0,3*np.pi)
X,Y = np.meshgrid(x,x)
f = lambda x,y, alpha, beta :(np.sin(X+alpha)+np.sin(Y*(1+np.sin(beta)*.4)+alpha))**2
alpha=np.linspace(0, 2*np.pi, num=34)
levels= 10
cmap=plt.cm.magma


fig, ax=plt.subplots()
props = dict(boxstyle='round', facecolor='wheat')
timelabel = ax.text(0.9,0.9, "", transform=ax.transAxes, ha="right", bbox=props)
t = np.ones(10)*time.time()
p = [ax.contourf(X,Y,f(X,Y,0,0), levels, cmap=cmap ) ]


def update(i):
    for tp in p[0].collections:
        tp.remove()
    p[0] = ax.contourf(X,Y,f(X,Y,alpha[i],alpha[i]), levels, cmap= cmap) 
    t[1:] = t[0:-1]
    t[0] = time.time()
    timelabel.set_text("{:.3f} fps".format(-1./np.diff(t).mean())+', hr='+"%.2f" % float(i))
    x1 = p[0].collections+[timelabel]
    return p[0].collections+[timelabel]

ani = matplotlib.animation.FuncAnimation(fig, update, frames=len(alpha), 
                                         interval=10, blit=True, repeat=True)
#plt.show()
ani.save('celluloid_minimal.gif', writer = 'imagemagick')
print('pp2',p)
```

```python
# use the camera package to create an animation. 
# I haven't been able to write an animated text stream to the animation yet so I am not using this technique
#
import datetime
from celluloid import Camera

print('Tin.dims',Tin.dims)
lev = DS['lev']
#print('lev.values',lev.values)
time = DS['time'].values
#print('time.values',time)
dtime = (time - time[0])/datetime.timedelta(hours=1) # convert time variable to hours as a float
levs = lev.values

l = len(time)
print('l=',l)

data2_hztf = data2.flatten()
data2_hztfd = data2_hztf[np.isfinite(data2_hztf)]
zmin = data2_hztfd.min()
zmax = data2_hztfd.max()
print('zmin,zmax', zmin, zmax)
clevs = findNiceContours(np.array([zmin,zmax]),nlevs=10)
print('clevs',clevs)
i,j,k = np.where(data2 == zmax) # coord where data2 takes on max val
print('max at ', i, j, k)
print('data2 max', data2[i[0],j[0],k[0]])

xv,yv=np.meshgrid(xc,levs)
fig, axes = plt.subplots(nrows=2)
camera = Camera(fig)

props = dict(boxstyle='round', facecolor='wheat')
timelabel = axes[0].text(0.4,0.4, "", transform=axes[0].transAxes, ha="right", bbox=props)

for tind in range(l):

    data2_hz = data2[:,tind,:].transpose()
    pl = axes[1].contourf(xv, yv, data2_hz, clevs, vmin=zmin, vmax=zmax,  extend='both')
    # Add colorbar to plot
    PRECT = DS['PRECT'+regtag].isel(time=tind)
    PRECT = PRECT*8.64e7
    PRECT.attrs['units'] = 'mm/day'
    vals = interp_ap(xc,yc,PRECT.values,lat.values,lon.values)
    axes[0].plot(xc,vals,linewidth=2,color="red",)

    title = 'hr='+"%.2f" % dtime[tind]
    #axes[0].set_title(title)
    timelabel.set_text(title)  
    #print('title',title)


    PCONVT = DS['PCONVT'+regtag]
    PCONVT = PCONVT/100.
    PCONVT.attrs['units'] = 'hPa'
    pconvtvals = interp_ap(xc, yc, PCONVT.values,lat.values,lon.values)
    pconvtval0 = pconvtvals[:,tind]
    #print('pconvtvals', pconvtval0)
    axes[1].plot(xc,pconvtval0,linewidth=2,color="white")

    PRECC = DS['PRECC'+regtag].isel(time=tind)
    PRECC = PRECC*8.64e7
    PRECC.attrs['units'] = 'mm/day'
    vals = interp_ap(xc,yc,PRECC.values,lat.values,lon.values)
    #print('PRECC',vals)
    axes[0].plot(xc,vals,linewidth=2,color="green")
    camera.snap()
 
cb = plt.colorbar(
        pl, orientation='horizontal',ticks=clevs,
        label='%s (%s)'%(Tin.long_name, Tin.units), pad=0.25, ax=axes[1]
    )
axes[1].set_ylim([1000., 100.])
axes[0].set_xlim([xc.min(), xc.max()])
axes[1].set_xlim([xc.min(), xc.max()])
axes[0].legend(['total','convective'])

animation = camera.animate()
animation.save('celluloid_minimal.gif', writer = 'imagemagick')
```

```python
# partial implementation, not satisfactory
import datetime
import matplotlib.animation as animation
#from celluloid import Camera

print('Tin.dims',Tin.dims)
lev = DS['lev']
#print('lev.values',lev.values)
time = DS['time'].values
#print('time.values',time)
dtime = (time - time[0])/datetime.timedelta(hours=1) # convert time variable to hours as a float
levs = lev.values

l = len(time)
print('l=',l)

data2_hztf = data2.flatten()
data2_hztfd = data2_hztf[np.isfinite(data2_hztf)]
zmin = data2_hztfd.min()
zmax = data2_hztfd.max()
print('zmin,zmax', zmin, zmax)
clevs = findNiceContours(np.array([zmin,zmax]),nlevs=10)
print('clevs',clevs)
i,j,k = np.where(data2 == zmax) # coord where data2 takes on max val
print('max at ', i, j, k)
print('data2 max', data2[i[0],j[0],k[0]])

xv,yv=np.meshgrid(xc,levs)
fig, axes = plt.subplots(nrows=2)

def draw(tind):   
    #print('tind',tind,dtime[tind])
    #ax = plt.axes()
    
    data2_hz = data2[:,tind,:].transpose()
    pl = axes[1].contourf(xv, yv, data2_hz, clevs, vmin=zmin, vmax=zmax,  extend='both')
    # Add colorbar to plot
    PRECT = DS['PRECT'+regtag].isel(time=tind)
    PRECT = PRECT*8.64e7
    PRECT.attrs['units'] = 'mm/day'
    vals = interp_ap(xc,yc,PRECT.values,lat.values,lon.values)
    axes[0].plot(xc,vals,linewidth=2,color="red",)

    title = 'hr='+"%.2f" % dtime[tind]
    axes[0].set_title(title)
    print('title',title)


    PCONVT = DS['PCONVT'+regtag]
    PCONVT = PCONVT/100.
    PCONVT.attrs['units'] = 'hPa'
    pconvtvals = interp_ap(xc, yc, PCONVT.values,lat.values,lon.values)
    pconvtval0 = pconvtvals[:,tind]
    #print('pconvtvals', pconvtval0)
    axes[1].plot(xc,pconvtval0,linewidth=2,color="white")

    PRECC = DS['PRECC'+regtag].isel(time=tind)
    PRECC = PRECC*8.64e7
    PRECC.attrs['units'] = 'mm/day'
    vals = interp_ap(xc,yc,PRECC.values,lat.values,lon.values)
    #print('PRECC',vals)
    axes[0].plot(xc,vals,linewidth=2,color="green")

 
cb = fig.colorbar(
        pl, orientation='horizontal',ticks=clevs,
        label='%s (%s)'%(Tin.long_name, Tin.units), pad=0.25, ax=axes[1]
    )
axes[1].set_ylim([1000., 100.])
axes[0].set_xlim([xc.min(), xc.max()])
axes[1].set_xlim([xc.min(), xc.max()])
axes[0].legend(['total','convective'])


def animate(i):
    draw(i)
    return axes,


draw(0)
print('np',np.arange(1,4))
#ani = animation.FuncAnimation(fig, animate, np.arange(1, 6), interval=10, blit=True)
#ani.save('celluloid_minimal.gif', writer = 'imagemagick')
```
