---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.14.4
  kernelspec:
    display_name: Python [conda env:.conda-pjrpy3] *
    language: python
    name: conda-env-.conda-pjrpy3-py
---

```python
import sys
print(sys.version)
%matplotlib inline
#from xhistogram.xarray import histogram
%run -i ~/Python/pjr3
```

```python
Varlist = np.array(['T','Q','CLOUD','CLDLIQ','ICWMR','CLDICE','RELHUM','NUMICE','NUMLIQ','Mass_bc'])
#                    RESTOM','FLNT','FSNT','TS','TMQ','PRECT','AEROD_v','CLDLOW','CLDTOT','LWCF','SWCF','TGCLDIWP','TGCLDLWP','SHFLX','LHFLX','PBLH','PCONVT','PRECC','PRECS'])
#Varlist = np.array(['TS','TMQ','PRECT'])
#Varlist = np.array(['RESTOM','LWCF','SWCF','FLNT','FSNT'])

Varname = 'TS' # Varlist[0]

"""
case_start1 = "/home/jupyter-haruki/work/CESM_MCB/Fixed_SST/"
case_start1 = "/scratch2/PJR/haruki_workdir/CESM_MCB/Fixed_SST/"
case_end1 = "Fixed_SST.cam.h0."+Varname+".y1-19.nc"
ind1 = case_start1+case_end1
DS1 = xr.open_mfdataset(ind1)
DS1 = center_time(DS1)

case_start2 = "/home/jupyter-haruki/work/CESM_MCB/MCB_R1R2R3_CN375cm/" 
case_start2 = "/scratch2/PJR/haruki_workdir/CESM_MCB/MCB_R1R2R3_CN600cm/"
case_end2 = "MCB_R1R2R3_CN600cm.cam.h0."+Varname+".y1-10.nc"
ind2 = case_start2+case_end2
DS2 = xr.open_mfdataset(ind2)
DS2 = center_time(DS2)
"""

caseid1 = "F2010_maint2.0_MCB-SSLT-EM_CNTL_test01"
case_start1 = "/compyfs/d3x345/climo_files/"+caseid1+"/"
case_end1 = caseid1+"_01_200401_200401_climo_fv192x288.nc"
#case_end1 = caseid1+"_01_200401_200401_climo.nc"
ind1 = case_start1+case_end1
#ind1 = '/compyfs/d3x345/vd05_ANN_climo.nc'
#print('ind1',ind1)
DS1 = xr.open_mfdataset(ind1)
#print('DS1',DS1)
#DS1 = center_time(DS1)

caseid2 = "F2010_maint2.0_MCB-SSLT-EM_CNTL_test01"
case_start2 = "/compyfs/d3x345/climo_files/"+caseid1+"/"
case_end2 = caseid2+"_01_200401_200401_climo_fv192x288.nc"
ind2 = case_start2+case_end2
DS2 = xr.open_mfdataset(ind2)
#DS2 = center_time(DS2)

Var1 = DS1[Varname]
#print('xxx',Var1)
Var2 = DS2[Varname]
#print('yyy',Var2)

```

```python tags=[]
def xr_llhplot(xrVar, plotproj=None, ax=None, cax=None,ylabels=None,clevs=None, cmap=None, title=None):
    """xr_llhplot xarray lat lon horizontal plot
    """
    #print(' entering xr_llhplot', xrVar)
    
    lon=xrVar['lon'].values
    lat=xrVar['lat'].values
    xv,yv=np.meshgrid(lon,lat)
    data_regridded = xrVar.values
    #print('aaa',data_regridded.shape, xv.shape, yv.shape)
    df = data_regridded.flatten()
    dsub = df[np.isfinite(df)] # ignore NaN
    zmax = dsub.max()
    zmin = dsub.min()
    #print('masked interpolated range',zmin,zmax)
    dataproj=ccrs.PlateCarree()    # data is always assumed to be lat/lon
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
    clat = (lat.min()+lat.max())/2.
    clon = (lon.min()+lon.max())/2.
    if plotproj is None:
        plotproj = ccrs.PlateCarree()
    #ax.set_extent([lon.values.min(), 260., lat.values.min(), lat.values.max()])
    #ax.set_global()
    #print('plotproj is ',plotproj)
    #rint('ax',ax)
 
    # if no ax argument, could get current axis, or create it
    if ax is None:
        #print('grab current axis')
        #ax = plt.gca()
        ax = plt.axes(projection=plotproj)

    if cax is None: cax = ax
    pl = ax.contourf(xv, yv, data_regridded, levels=clevs, # vmin=zmin, vmax=zmax,
                     norm=norm, cmap=cmap,
                     extend=extend, transform=ccrs.PlateCarree())

    # Add colorbar to plot
    cb = plt.colorbar(
        pl, orientation='horizontal',ticks=clevs,ax=cax,
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

Var = Var1.mean(dim='time',keep_attrs=True)

plotproj = ccrs.PlateCarree()
#plotproj=ccrs.Mercator()
#plotproj=ccrs.Orthographic(central_latitude=0,central_longitude=55)   # any projections should work


if True:
    fig, axes = plt.subplots(ncols=3
                         ,gridspec_kw={'width_ratios': [1, 1, 1]}
                         ,subplot_kw={'projection': plotproj}
                         ,figsize=(16,5)
                        )
 
xr_llhplot(Var,ax=axes[1])
#xr_llhplot(Var)

```

```python
regtag = ""
weights = None

Varlist = np.array(['RESTOM','FLNT','FSNT','TS','TMQ','PRECT','AEROD_v','CLDLOW','CLDTOT','LWCF','SWCF','TGCLDIWP','TGCLDLWP',
                    'SHFLX','LHFLX','PBLH','PCONVT','PRECC','PRECS'])
Varlist = np.array(['FLNT','FSNT','TS','PRECC','PRECL','AODVIS','CLDLOW','CLDTOT','LWCF','SWCF','TGCLDIWP','TGCLDLWP',
                    'SHFLX','LHFLX','PBLH','PRECSC','PRECSL'])
#Varlist = np.array(['TS','TMQ','PRECT'])
#Varlist = np.array(['RESTOM','LWCF','SWCF','FLNT','FSNT'])
#Varlist = np.array(['AEROD_v'])
#Varlist = np.array(['FSNT','TGCLDLWP'])
#Varlist = np.array(['AODVIS'])

case_start1 = "/home/jupyter-haruki/work/CESM_MCB/MCB_R1R2R3_CN375cm/MCB_R1R2R3_CN375cm.cam.h0." 
case_start1 = "/home/jupyter-haruki/work/CESM_MCB/MCB_R1R2R3_CN600cm/MCB_R1R2R3_CN600cm.cam.h0."
case_start1 = "/scratch2/PJR/haruki_workdir/CESM_MCB/Fixed_SST/"
case_end1 = ".y1-10.nc"
pref1='expCN600'

case_start2 = "/home/jupyter-haruki/work/CESM_MCB/Fixed_SST/Fixed_SST.cam.h0."
case_start2 = "/scratch2/PJR/haruki_workdir/CESM_MCB/MCB_R1R2R3_CN600cm/"
case_end2 = ".y1-19.nc"
pref2='control'

case_start1 = "/home/jupyter-haruki/work/CESM_MCB/MCB_R1R2R3_CN375cm/MCB_R1R2R3_CN375cm.cam.h0." 
case_start1 = "/scratch2/PJR/haruki_workdir/E3SM_MCB/20221018.v2.LR.F2010.E1_R1-3_CDNC600.eam.h0.y1-5.FORCING.nc"
case_end1 = ""
pref1='expCN600'

case_start2 = "/scratch2/PJR/haruki_workdir/E3SM_MCB/20220930.v2.LR.F2010.E1_CNTL.eam.h0.y1-14.FORCING.nc"
case_end2 = ""
pref2='control'

caseid1 = "F2010_maint2.0_MCB-SSLT-EM_CNTL_test01"
case_start1 = "/compyfs/d3x345/climo_files/"+caseid1+"/"+caseid1+"_01_200401_200401_climo_fv192x288.nc"
case_end1 = "" 

caseid2 = "F2010_maint2.0_MCB-SSLT-EM_CNTL_test01"
case_start2 = "/compyfs/d3x345/climo_files/"+caseid2+"/"+caseid2+"_01_200401_200401_climo_fv192x288.nc"
case_end2 = ""

for Varname in Varlist:
    print()
    print('-------------------------------'+Varname)
    
    ind1 = case_start1+Varname+case_end1
    ind1 = case_start1
    print('opening',ind1)
    DS1 = xr.open_mfdataset(ind1)
    #DS1 = center_time(DS1)
    Var1 = xr_getvar(Varname,DS1)
    V1 = Var1.mean(dim='time',keep_attrs=True)


    ind2 = case_start2+Varname+case_end2
    ind2 = case_start2
    print('opening',ind2)
    DS2 = xr.open_mfdataset(ind2)
    #DS2 = center_time(DS2)
    Var2 = xr_getvar(Varname,DS2)
    #print('yyy',Var2)
    V2 = Var2.mean(dim='time',keep_attrs=True)

    DV = V1-V2
    
    if weights is None:
        lat = Var1['lat'].values
        lon = Var1['lon'].values
        area = make_fvarea(lon,lat)
        weights = V1.copy()
        weights.data =area
        weights.attrs['units']='steradians'
        #area = xr_getvar('area',DS1).where(pmask)
        #print('weights',weights)
        #print('weights shape',weights.shape)

    
    print(Varname, V1.attrs['long_name'],'Range V1 and V2 ',V1.min().values, V1.max().values, V2.min().values, V2.max().values)
    V1A = V1.weighted(weights).mean()
    V2A = V2.weighted(weights).mean()
    DVA = V1A-V2A
    print('V1A %5.3f' % (V1A.values),' V2A %5.3f' % (V2A.values),' DVA %5.3f' % (DVA.values))

    if V1.min().values == V1.max().values:
        print('constant field skipping plot ')
    else:
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

        xr_llhplot(V1, ax=axes[0],clevs=clevs,title=pref1)
        xr_llhplot(V2, ax=axes[1],clevs=clevs,ylabels=False,title=pref2)
        xr_llhplot(DV, ax=axes[2],clevs=dlevs,cmap=dmap,title=pref1+'-'+pref2)
        plt.savefig('test_'+Varname+'.jpg',format='jpg')
        plt.show()
        
    print('field processing complete')

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
