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

**compare UKESM SYN (Synthetic) estimate of fields to E3SM and CESM CON (Concurrent) syms**


```python
import sys
print(sys.version)
%matplotlib inline
%run -i ~/Python/pjr3
```

```python
def setfig ():
    """
    return fig and axes for a single panel figure
    """
    plotproj = ccrs.Mollweide()
    plotproj._threshold /= 100.
    fig, axes = plt.subplots(ncols=3,nrows=2,
                             #gridspec_kw={'width_ratios': [1]},
                             subplot_kw={'projection': plotproj},
                             figsize=(8.,4.5),
                            )
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.02, hspace=0.05)

    fig.set_dpi(300.0)
    return fig, axes;

def pltllbox2(xri, yri,ax=None):
    if xri[1] < xri[0]:
        xri[1] += 360.
    regcx = [xri[0],xri[1],xri[1],xri[0],xri[0]]
    regcy = [yri[0],yri[0],yri[1],yri[1],yri[0]]
    if ax is None:
        ax = gca()
    ax.plot(regcx,regcy,color='red',transform=ccrs.PlateCarree())
```

```python
def make_AA_tser(Var):
    ''' make Annual Average from tseries files
        Var: Xarray object
        returns object with annual averages
    '''
    month_length = Var.time.dt.days_in_month
    twgts = month_length.groupby("time.year") / month_length.groupby("time.year").sum()
    V1AY = (Var*month_length).groupby("time.year").sum()/month_length.groupby("time.year").sum()
    #print('V1AY.values', V1AY.values)
    return V1AY
```

```python
def xr_llhplot2 (xrVar, cbar='default', plotproj=None, ax=None, cax=None,
                 ylabels=None, clevs=None, cmap=None, title=None, cbartitle=None,regmark=False):
    """xr_llhplot xarray lat lon horizontal plot
       returns a "mappable" containing the artist info that is needed to plot a colorbar
       that is
       mpbl = xr_llhplot2()
       plt.colorbar(mpbl, orientation='horizontal',cax=cax,...)
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
        #cmap = mpl.cm.get_cmap()
        cmap = plt.get_cmap()
        #print('bbb',cmap.N)
    #print('cmap',cmap)
    extend = 'both'
    norm = mpl.colors.BoundaryNorm(clevs,cmap.N,extend=extend)
    #print('norm',norm(clevs))
    clat = (lat.min()+lat.max())/2.
    clon = (lon.min()+lon.max())/2.
    if plotproj is None:
        plotproj = ccrs.PlateCarree()
        plotproj = ccrs.Mollweide()
 
    # if no ax argument, could get current axis, or create it
    if ax is None:
        #print('grab current axis')
        #ax = plt.gca()
        ax = plt.axes(projection=plotproj)
        
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,
                      linewidth=2, color='gray', alpha=0.5)
    pl = ax.contourf(xv, yv, data_regridded, levels=clevs, # vmin=zmin, vmax=zmax,
                     norm=norm, cmap=cmap,
                     extend=extend, transform=ccrs.PlateCarree())
    
    gl.left_labels=ylabels
    gl.right_labels=ylabels
    ax.coastlines(linewidth=1,color='blue')
 
    ## Find the location of the main plot axes
    ## has to be done after some plotting is done in projection space
    posn = ax.get_position()
    
    ax2 = fig.add_axes([0,0,0.1,0.1])
    ax2.set_position([posn.x0-0.005, posn.y0-0.005, posn.width+0.01, posn.height+0.01])
    ax2.patch.set_alpha(0.0)
    ax2.set_axis_off()
    ax2.set_xlim([0,1])
    ax2.set_ylim([0,1])
    if regmark:
        # print some registration marks to help in lining up figures
        ax2.scatter([0,0,1,1], [0,1,0,1], c="r", s=100)

    if not title is None:
        #ax2.set_title(title)
        ax2.text(0.01,0.93,title,fontsize=6)

        
    if cbar == 'default':
        # Add colorbar to plot
        if cbartitle is None:
            cbartitle = xrVar.long_name
            
        if cax is not None:
            cax = ax
        else:
            # create an colorbar axis
            cax = fig.add_axes([0,0,0.1,0.1])
            ## Adjust the positioning and orientation of the colorbar
            xof = 0.15
            cax.set_position([posn.x0-xof, posn.y0-0.05, posn.width+2*xof, 0.035])
            #cax.set_position([posn.x0, posn.y0-0.02, posn.width, 0.015])
    
        cb = plt.colorbar(
             pl, orientation='horizontal',ticks=clevs,cax=cax, 
             #label='%s (%s)'%(cbartitle, xrVar.units)
             )
        cb.ax.tick_params(labelsize=7)
        cax.text(0.,-1.6,'%s (%s)'%(cbartitle, xrVar.units),va='center',ha='center')
    
        
    return pl
```

```python
#Var = DS0.sel(time=slice("2020-01-01", "2030-01-01")).FSNT
#print(Var)
#VarAA = make_AA_tser(Var)
#print(VarAA)
```

```python
# begin cells for CESM fv grid
```

```python
casename0 = "b.e21.BSSP245smbb.f09_g17.001"  # reference run
casename1 = "b.e21.BSSP245smbb_MCBss7TgYr_R1R2R3.f09_g17.LE2-1011.001"
pref1='R1'
prefmod='CESM_Coupled_MCBSSE123'

Varname='FSNT'
#Varname='SWCF'
Varname='TS'
#Varname='PRECT'


def bld_fname2(casename, Varname):
    fname = "/e3sm_prod/phil/timeseries/cesm2-mcb-reshaped/"+casename+"/"+casename+".cam.h0.2015-*."+Varname+".nc"
    return fname

def bld_fname3(casename, Varname):
    fname = "/e3sm_prod/phil/timeseries/cesm2-mcb-reshaped/CESM2_SSP245/ens001/"+casename+".cam.h0."+Varname+".201501-206412.nc"
    return fname

ind0 = bld_fname3(casename0, Varname)
print('ind0',ind0)
DS0 = center_time(xr.open_mfdataset(ind0))
ind1 = bld_fname2(casename1, Varname)
print('ind1',ind1)
DS1 = center_time(xr.open_mfdataset(ind1))

#print('DS0',DS0)
#print('DS1',DS1)
DS0, DS1 = reconcile_xr_coords(DS0, DS1)
print('DS0.time',len(DS0.time.values))
print('DS1.time',len(DS1.time.values))
# grab part of the data
#DS0 = DS0.isel(time=slice(0,120))
yb = '2020-01-01'
ye = '2030-01-01'
DS0 = DS0.sel(time=slice(yb,ye))
DS1 = DS1.sel(time=slice(yb,ye))


long_name = None
if Varname == 'PRECT': long_name = 'Precip'
if Varname == 'TS': long_name = "Surface Temperature"
    
C0 = xr_getvar(Varname, DS0,long_name=long_name).mean(dim="time")
C1 = xr_getvar(Varname, DS1,long_name=long_name).mean(dim="time")

if 'area' in DS1:
    area = DS1['area']
else:
    print('calculating weights')
    lat = DS1['lat'].values
    lon = DS1['lon'].values
    aread = make_fvarea(lon,lat)
    area = xr.DataArray(aread, dims=['lat','lon'], coords={'lon':lon,'lat':lat})
    area.attrs['units']='steradians'
    #print('area',area)

wdims = area.dims
#print('wdims',wdims)
weights = area/(area.sum(dim=wdims))
#print('weights sum',weights.shape,weights.sum(dim=wdims).values)
    

#dlevs=None
#clevs=None

C0A = C0.weighted(weights).mean()
print("C0A",C0A.values)

C1A = C1.weighted(weights).mean()
C1AD = C1A - C0A
#print('C1A,V2A',C1A.values,V2A.values)
#print(C1AD)
sC1AD = ' (%5.2f)' % C1AD
dlevs = None

fig, axes = setfig()
print('axes', axes.shape, axes)
axf = axes.flatten()

D1 = C1-C0

print('min',D1.min().values)
drmin = np.min([D1.min()])
drmax = np.max([D1.max()])
print('drange', drmin, drmax)
factor = 0.9
dlevs = findNiceContours(np.array([drmin,drmax])*factor,nlevs = 15,rmClev=0.,sym=True)
print("dlevs",dlevs)
dmap = diverge_map()

if True:
    xr_llhplot2(D1, ax=axf[1],clevs=dlevs,cmap=dmap,title=pref1+sC1AD, ylabels=False)#,cbar=None)
    pltllbox2([-150.,-110.],[0.,30.],ax=axf[1])
    pltllbox2([-110.,-70.],[-30.,0.],ax=axf[1])
    pltllbox2([-25.,15.],[-30.,0.],ax=axf[1])


#plt.savefig(prefmod+'_'+Varname+'-D.pdf',format='pdf',dpi=300,transparent=True)#,facecolor='xkcd:mint green')
plt.show()
```

```python
# begin cells for UKESM
```

```python
def xr_getvar_sl(VN, DS1, method='surface', verbose=False):
    """ get a field from a netCDF file.
    If it is a multi-level field, do something useful to provide single level information
    This function assumes that the x-xoordinate increases monotonically

    """
    Var1 = xr_getvar(VN,DS1)
    #V1 = Var1.mean(dim='time',keep_attrs=True)
    dimlist = Var1.dims
    #print('dimlist',dimlist)
    if 'model_level_number' in dimlist:
        level_height = xr_getvar('level_height',DS1)
        sigma = xr_getvar('sigma',DS1)
        #print('sigma',sigma.values)
        surface_altitude = xr_getvar('surface_altitude',DS1)
        altitude = level_height + (sigma * surface_altitude)
        altitude.attrs['long_name'] = 'altitude above mean sea-level'
        #print('altitude',altitude)
        if method == 'surface':
            print('method:surface')
            V1 = Var1.isel(model_level_number=0)
            V1.attrs['long_name'] = V1.attrs['long_name'] + ' (surface level)'
        elif method == 'maxb850':
            i = 0
            j = 0
            #print('method:max below 850')
            dz1 = altitude - surface_altitude   # height above surface
            #print('altitude',altitude[:,i,j].values)
            #print('surface_altitude',surface_altitude[i,j].values,surface_altitude)
            #print('dz1(i,j)',dz1[:,i,j].values)
            dz2 = (sigma - 1)*surface_altitude + level_height  # height above sea level
            #print('dz2(i,j)',dz2[:,i,j].values)
            pmb = 850.
            psmb = 1000.
            scaleheight = 8.4e3
            altmb = -np.log(pmb/psmb)*scaleheight
            #print('altmb is ', altmb)
            V1 = Var1.copy()
            #V1.values = dz1.values ### for the moment, overwrite field with height above surface height
            #print('V1',V1[:,i,j].values)
            #V2 = V1.where(dz1 <= altmb+50.)
            V2 = V1.where(altitude <= altmb+50.)
            V3 = V2.max(dim='model_level_number')
            #V3.attrs['long_name'] = 'max value below 850hPa of '+V2.attrs['long_name']
            V3.attrs['long_name'] = 'max value below 850hPa of '+V2.name
            V1 = V3
            #print('V3',V3[i,j].values)
            #1./0.
    else:
        V1 = Var1
    return V1
```

```python

def fix_UKMO_ds(filename, dsin: xr.Dataset) -> xr.Dataset:
    """
    Rename some variables and rescale 
    :param ds: some xarray dataset
    :return: a normalized xarray dataset, or the original one
    """

    #print('filename is ', filename)
    name_dict = dict()
    ds = dsin.copy()
    ds = _normalize_lat_lon(ds) # harmonize the lat, lon fields
    ds.coords['lon'] = (ds.coords['lon'] + 180) % 360 - 180
    ds = ds.sortby(ds.lon)
    #print('ds',ds)
    for Vname in ds:
        #print('Vname',Vname)
        if (('PBL' in filename) & ('height_after' in Vname)):
            ds[Vname].attrs['UKESM name'] = Vname
            ds[Vname].attrs['units'] = 'kg/m2'
            ds[Vname].attrs['long_name'] = 'PBL height'
            ds = ds.rename({Vname:'PBLH'})
            #print('fixed PBLH')
            
        if (('LWP' in filename) & (Vname == 'm01s02i391')):
            ds[Vname].attrs['UKESM name'] = Vname
            ds[Vname].attrs['units'] = 'kg/m2'
            ds[Vname].attrs['long_name'] = 'Liq Water Path'
            ds = ds.rename({Vname:'TGCLDLWP'})
            #print('fixed LWP')
            
        if (('AOD' in filename) & (Vname == 'unknown')):
            ds[Vname].attrs['UKESM name'] = Vname
            ds[Vname].attrs['units'] = '1'
            ds[Vname].attrs['long_name'] = 'AOD 550nm'
            ds = ds.rename({Vname:'AOD'})
    
        if (('_r_e' in filename) & (Vname == 'unknown')):
            ds[Vname].attrs['UKESM name'] = Vname
            ds[Vname].attrs['long_name'] = 'Effective Radius at Cloud-top'
            ds[Vname].attrs['units'] = '$\mu m$'
            ds = ds.rename({Vname:'REPJR'})
            #print('fixed RE')

        if (('precip_rate' in filename) & (Vname == 'precipitation_flux')):
            ds[Vname].attrs['UKESM name'] = Vname
            ds[Vname].attrs['units'] = 'm/s'
            ds[Vname] = ds[Vname]/1000.
            ds[Vname].attrs['long_name'] = 'Total Precipitation'
            ds = ds.rename({Vname:'PRECT'})
            #print('fixed PRECT')
           
        if (('p_surface_Pa' in filename) & (Vname == 'surface_air_pressure')):
            ds[Vname].attrs['UKESM name'] = Vname
            ds[Vname] = ds[Vname]/100.
            ds[Vname].attrs['units'] = 'hPa'
            ds[Vname].attrs['long_name'] = ds[Vname].attrs['standard_name']
            ds = ds.rename({Vname:'PS'})
            #print('fixed PS')
        
        if (('cloud_fraction' in filename) & (Vname == 'cloud_area_fraction_in_atmosphere_layer')):
            ds[Vname].attrs['UKESM name'] = Vname
            #ds[Vname] = ds[Vname]/100.
            #ds[Vname].attrs['units'] = 'hPa'
            #ds[Vname].attrs['long_name'] = ds[Vname].attrs['standard_name']
            ds = ds.rename({Vname:'CLOUD'})
        
        if (('T_surface_K' in filename) & (Vname == 'surface_temperature')):
            ds[Vname].attrs['UKESM name'] = Vname
            ds[Vname].attrs['long_name'] = ds[Vname].attrs['standard_name']
            ds = ds.rename({Vname:'TS'})
            #print('fixed TS')
            
        if (('net_ToA_SW_W' in filename) & (Vname == 'unknown')):
            ds[Vname].attrs['UKESM name'] = Vname
            ds[Vname].attrs['units'] = 'W/m2'
            ds[Vname].attrs['long_name'] = 'Net TOA Shortwave'
            ds = ds.rename({Vname:'FSNT'})
            #print('fixed FSNT')
            
        if (('net_ToA_SW_C' in filename) & (Vname == 'unknown')):
            ds[Vname].attrs['UKESM name'] = Vname
            ds[Vname].attrs['units'] = 'W/m2'
            ds[Vname].attrs['long_name'] = 'Net TOA Shortwave Clear'
            ds = ds.rename({Vname:'FSNTC'})
            #print('fixed FSNT'
            
        if (('_SW_Clear' in filename) & (Vname == 'toa_outgoing_shortwave_flux_assuming_clear_sky')):
            ds[Vname] = -ds[Vname]
            ds[Vname].attrs['UKESM name'] = Vname
            ds[Vname].attrs['units'] = 'W/m2'
            ds[Vname].attrs['long_name'] = 'outgoing SW assuming clearsky (upward -ive)'
            ds = ds.rename({Vname:'MFSUTC'})
              

    return ds


```

```python
Vdict = {'net_ToA_LW_W_m2':'toa_outgoing_longwave_flux'
        ,'net_ToA_SW_W_m2':'FSNT'
        ,'net_ToA_SW_Clear_W_m2':'FSNTC'
        ,'AOD_550nm':'AOD'
        ,'LWP_kg_m2':'TGCLDLWP'
        ,'p_surface_Pa':'PS'
        ,'T_surface_K':'TS'
        ,'precip_rate_kg_m2_sec':'PRECT'
        ,'PBL_depth_metres':'PBLH'
        ,'cloudtop_r_e_microns':'REPJR'
        ,'cloud_fraction':'CLOUD'
        ,'Outgoing_SW_Clear_W_m2':'MFSUTC'
        }


```

```python
pref_fn = 'UKESM'

def make_ind1 (REG_ID, Varname, filetype=None):
    case_start1 = '~/NetCDF_Files/UKESM1_data/'+REG_ID+'_20450101_20490101_mean_'
    case_start1 = '~/NetCDF_Files/UKESM1_data_v2/Coupled_50Tg/'+REG_ID+'_coupled_50Tgy_20410101_20500101_mean_'
    case_end1 = ".nc"
    fstring1 ='%s%s%s' 
    pref1=REG_ID+'_UKESM1_50Tgpy_Cpld'
    
    if filetype == 'Fixed_SST':
        # fixed SST simulations
        case_start1 = '~/NetCDF_Files/UKESM1_data/'+REG_ID+'_AtmosOnly_19840101_19880101_mean_'
        pref1=REG_ID+'_50Tgpy_FixedSST'
        case_start1 = '~/NetCDF_Files/UKESM1_data_v2/AtmosOnly_25Tg_1979-1989/'+REG_ID+'_AtmosOnly_25Tgy_19790101_19890101_mean_'
        pref1=REG_ID+'_25Tgpy_FixedSST'
        #case_start1 = '~/NetCDF_Files/UKESM1_data_v2/AtmosOnly_50Tg_1979-1989/'+REG_ID+'_AtmosOnly_50Tgy_19790101_19890101_mean_'
        #pref1=REG_ID+'_50Tgpy_FixedSST'
        case_end1 = ".nc"
        fstring1 ='%s%s%s' 
        pref1=REG_ID+'_UKESM1_25Tgpy_Cpld'

    ind1 = fstring1 % (case_start1,Varname,case_end1)
    return ind1

def make_ind2(REG_ID, Varname, filetype='Fixed_SST'):    
    case_start2 = '~/NetCDF_Files/UKESM1_data/CTL_20450101_20490101_mean_'
    case_start2 = '~/NetCDF_Files/UKESM1_data_v2/Coupled_Control/CTL_coupled_20410101_20500101_mean_'
    case_end2 = ".nc"
    pref2='UKESM1_control'
    fstring2 ='%s%s%s' 

    if filetype == 'Fixed_SST':
        case_start2 = '~/NetCDF_Files/UKESM1_data/CTL_AtmosOnly_19840101_19880101_mean_'
        case_start2 = '~/NetCDF_Files/UKESM1_data_v2/AtmosOnly_Control_1979-1989/'+'CTL_AtmosOnly_19790101_19890101_mean_'
        case_end2 = ".nc"
        pref2='Control'
        fstring2 ='%s%s%s' 
        
    ind2 = fstring2 % (case_start2,Varname,case_end2)
    return ind2


Varname='LWP_kg_m2'
REG_ID = 'R1_NEP'
filetype = 'Fixed_SST'
#filetype = 'Coupled'
ind1 = make_ind1(REG_ID,Varname,filetype)

print('example string used for file 1 open',ind1)
ind2 = make_ind2(REG_ID,Varname,filetype)

difftitle='Syn R1+R2+R3(each@25Tgpyr)-Ctl'
difftitle=''
print('example string used for file 2 open',ind2)


```

```python
# make a "synthetic" estimate of the change in field by accumulating differences 
# set for calculation over 3 areas
# if both FSNT and FSNTC are requested we also calculate SWCRE

#Varlist = np.array(['p_surface_Pa'])
#Varlist = np.array(['Outgoing_SW_Clear_W_m2','p_surface_Pa','T_surface_K','precip_rate_kg_m2_sec','PBL_depth_metres','cloudtop_r_e_microns','AOD_550nm','LWP_kg_m2','net_ToA_LW_W_m2','net_ToA_SW_W_m2'])
Varlist = np.array(['Outgoing_SW_Clear_W_m2','precip_rate_kg_m2_sec','PBL_depth_metres','cloudtop_r_e_microns','AOD_550nm','LWP_kg_m2','net_ToA_LW_W_m2','net_ToA_SW_W_m2',
                    'net_ToA_SW_Clear_W_m2','cloud_fraction'])
Varlist = np.array(['net_ToA_SW_Clear_W_m2','net_ToA_SW_W_m2'])
#Varlist = np.array(['net_ToA_SW_Clear_W_m2','net_ToA_SW_W_m2'])
#Varlist = np.array(['AOD_550nm','LWP_kg_m2','net_ToA_LW_W_m2','net_ToA_SW_W_m2'])
Varlist = np.array(['T_surface_K'])
#Varlist = np.array(['LWP_kg_m2','net_ToA_SW_Clear_W_m2','net_ToA_SW_W_m2','cloud_fraction'])
#Varlist = np.array(['cloud_fraction'])
#Varlist = np.array(['net_ToA_SW_W_m2'])
Varlist = np.array(['precip_rate_kg_m2_sec'])

FSNT1 = None
FSNT2 = None
FSNTC1 = None
FSNTC2 = None

# specify regions (assume lon always specified as west, then east limit)
xreg = np.array([[-150.,-110.],[-110,-70],[-25.,15.],[170.,-120.],[-170.,-90.]])%360.
yreg = np.array([[0.,30.],     [-30.,0.], [-30.,0.], [30.,50.],   [-50.,-30.] ])
namereg = ['NEP','SEP','SEA','NP','SP']
#xreg = [[0.,360.]]
#yreg = [[-90.,91.]]


reglist = np.array(['R1_NEP','R2_SEP','R3_SEA'])
#reglist = np.array(['R2_SEP'])
#reglist = np.array(['R3_SEA'])
#reglist = np.array(['R1_NEP'])

filetype = None
filetype = 'Fixed_SST'
filetype = 'Coupled'

for Varname in Varlist:
    print()
    print('-------------------------------'+Varname)
    nreg = 0 # is it the start of the region summation
    for REG_ID in reglist:
        #ind1 = fstring1 % (case_start1,Varname,case_end1)
        ind1 = make_ind1(REG_ID,Varname,filetype)
        print('ind1 opening',ind1)
        DS1 = xr.open_mfdataset(ind1)
        #print('xxx',DS1.time_bnds.values)
        DS1 = fix_UKMO_ds(ind1, DS1)
        #print('DS1.lon',DS1.lon.values)
        #DS1 = center_time(DS1)
        VN = Vdict[Varname]
        print('VN is ',VN)
        V1 = xr_getvar_sl(VN,DS1,method='maxb850')
        #print('V1',V1)
        ind2 = make_ind2(REG_ID,Varname,filetype)
        print('opening ind2',ind2)
        #DS2 = xr.open_mfdataset(ind2)
        DS2 = xr.open_mfdataset(ind2)
        DS2 = fix_UKMO_ds(ind2, DS2)
        V2 = xr_getvar_sl(VN,DS2,method='maxb850')

        DV = V1-V2
        print('DV range', DV.min().values, DV.max().values)
        DVB = DV.sel(lon=-110.,lat=-15.,method='nearest')
        print('DVB',DVB.values)
        weights = None
        if 'area' in DS1:
            area = DS1['area']
        elif 'area' in DS2:
            area = DS2['area']
        else:
            print('calculating areas')
            lat = V1['lat'].values
            lon = V1['lon'].values
            area = make_fvarea(lon,lat)
        weights = V1.copy()
        weights.data =area
        weights.attrs['units']='steradians'

        print(Varname, V1.attrs['long_name'],'Range V1 and V2 ',V1.min().values, V1.max().values, V2.min().values, V2.max().values)
        V1A = V1.weighted(weights).mean()
        sV1A = ' (%5.2f)' % V1A
        V2A = V2.weighted(weights).mean()
        sV2A = ' (%5.2f)' % V2A
        DVA = V1A-V2A
        sDVA = ' (%5.2f)' % DVA
        print('nreg',nreg,'area avgs (V1A, V2A, DVA) %5.2f' % (V1A.values), '%5.2f' % (V2A.values),' Delta %5.2f' % (DVA.values))
        if nreg == 0:
            V1S = V1
            V2S = V2
            DVS = DV
            nreg = 1
        else:
            V1S = V1S + V1
            V2S = V2S +V2
            DVS = DVS + DV
            nreg = nreg+1
    print(' all regions accumulated')
    #V1S = V1S/nreg
    #V2S = V2S/nreg
    print('DVSB',DVS.sel(lon=-110.,lat=-15.,method='nearest').values)
    long_name = None
    if VN == 'PRECT': long_name = 'Precip'
    if VN == 'TS': long_name = "Surface Temperature"
    if not long_name is None:
        DVS.attrs['long_name'] = long_name
        print('aaa')
    
    DVA = DVS.weighted(weights).mean()
    sDVA = ' ({:5.2f}{:s})'.format(DVA.values, DVA.units)
    print('yyy', sDVA)
    print('accumulated area Delta %5.2f' % (DVA.values))
    fname = pref_fn+'_'+Varname+'_'+DV.name+'-D.pdf'

# now do some plotting
    fig, axes = setfig()
    axf = axes.flatten()

    prefsum = ''

    print('min',DVS.min().values)
    drmin = np.min([DVS.min()])
    drmax = np.max([DVS.max()])
    print('drange', drmin, drmax)
    factor = 0.8
    dlevs = findNiceContours(np.array([drmin,drmax])*factor,nlevs = 15,rmClev=0.,sym=True)
    print("dlevs",dlevs)
    dmap = diverge_map()

    if True:
        xr_llhplot2(DVS, ax=axf[1],clevs=dlevs,cmap=dmap,title=prefsum+sDVA, ylabels=False)#,cbar=None)
        pltllbox2([-150.,-110.],[0.,30.],ax=axf[1])
        pltllbox2([-110.,-70.],[-30.,0.],ax=axf[1])
        pltllbox2([-25.,15.],[-30.,0.],ax=axf[1])


    #plt.savefig(prefmod+'_'+Varname+'-D.pdf',format='pdf',dpi=300,transparent=True)#,facecolor='xkcd:mint green')
    plt.show()

    if VN == 'FSNT':
        print('save FSNT')
        FSNT1=V1S
        FSNT2=V2S
    elif VN == 'FSNTC':
        print('save FSNTC')
        FSNTC1=V1S
        FSNTC2=V2S
        
print('field processing complete')

if ((FSNT1 is None) or (FSNTC1 is None)):
    print ('fields for SWCRE not requested')
else:
    print (' calculating SWCRE')
    SWCRE1 = FSNT1-FSNTC1
    SWCRE1.attrs['long_name'] = 'TOA shortwave CRE'
    SWCRE1 = SWCRE1.rename('SWCRE')
    SWCRE2 = FSNT2-FSNTC2
    SWCRE2.attrs['long_name'] = 'TOA shortwave CRE'
    SWCRE2 = SWCRE2.rename('SWCRE')
    DSWCRE = SWCRE1-SWCRE2
    DSWCREA = DSWCRE.weighted(weights).mean()
    #sDSWCREA = ' (%5.2f %s)' % DSWCREA, DSWCREA.units
    sDSWCREA = ' ({:5.2f}{:s})'.format(DSWCREA.values, DSWCREA.units)
    print('accumulated SWCRE change',sDSWCREA)
    fname = pref_fn+'_'+Varname+'_'+'SWCRE'+'-D.pdf'
    #pltfld(DSWCRE, difftitle+sDSWCREA,fname)
    
    fig, axes = setfig()
    axf = axes.flatten()

    prefsum = 'Syn UKESM'

    print('min',DSWCRE.min().values)
    drmin = np.min([DSWCRE.min()])
    drmax = np.max([DSWCRE.max()])
    print('drange', drmin, drmax)
    factor = 0.8
    dlevs = findNiceContours(np.array([drmin,drmax])*factor,nlevs = 15,rmClev=0.,sym=True)
    print("dlevs",dlevs)
    dmap = diverge_map()
    
    xr_llhplot2(DSWCRE, ax=axf[4],clevs=dlevs,cmap=dmap,title=difftitle+sDSWCREA, ylabels=False)#,cbar=None)
    pltllbox2([-150.,-110.],[0.,30.],ax=axf[4])
    pltllbox2([-110.,-70.],[-30.,0.],ax=axf[4])
    pltllbox2([-25.,15.],[-30.,0.],ax=axf[4])
    plt.show()
```

```python
1./0.
```