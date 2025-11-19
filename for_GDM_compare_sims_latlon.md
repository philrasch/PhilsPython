---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.15.2
  kernelspec:
    display_name: Python [conda env:pjrpy3] *
    language: python
    name: conda-env-pjrpy3-py
---

**compare two cases over the globe assuming they are on lat/lon grid at same resolution**

*configured for the GMD paper*

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
    plotproj = ccrs.Mollweide(central_longitude=-80)
    plotproj._threshold /= 100.
    fig, axes = plt.subplots(ncols=1,
                             gridspec_kw={'width_ratios': [1]},
                             subplot_kw={'projection': plotproj},
                             figsize=(6,4.1),
                            )
    fig.set_dpi(300.0)
    return fig, axes;

def pltllbox(xri, yri,ax=None):
    if ax is None:
        ax = plt.gca()
    if xri[1] < xri[0]:
        xri[1] += 360.
    regcx = [xri[0],xri[1],xri[1],xri[0],xri[0]]
    regcy = [yri[0],yri[0],yri[1],yri[1],yri[0]]
    ax.plot(regcx,regcy,color='red',transform=ccrs.PlateCarree())
```

```python
def xr_llhplot2 (xrVar, cbar='default', plotproj=None, ax=None, cax=None, fig=None,
                 ylabels=False, clevs=None, cmap=None, title=None, cbartitle=None):
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
    
    # print some registration marks to help in lining up figures
    ax2 = fig.add_axes([0,0,0.1,0.1])
    ax2.set_position([posn.x0-0.005, posn.y0-0.005, posn.width+0.01, posn.height+0.01])
    ax2.patch.set_alpha(0.0)
    ax2.scatter([0,0,1,1], [0,1,0,1], c="r", s=100)
    ax2.set_axis_off()
    ax2.set_xlim([0,1])
    ax2.set_ylim([0,1])

    if not title is None:
        #ax.set_title(title)
        ax2.text(0.01,0.9,title)
    
    # Add colorbar to plot
    if cbartitle is None:
        cbartitle = xrVar.long_name
        
    if cbar == 'default':
        if cax is not None:
            cax = ax
        else:
            # create an colorbar axis
            cax = fig.add_axes([0,0,0.1,0.1])
            ## Adjust the positioning and orientation of the colorbar
            cax.set_position([posn.x0, posn.y0-0.07, posn.width, 0.05])

        cb = plt.colorbar(
             pl, orientation='horizontal',ticks=clevs,cax=cax,
             label='%s (%s)'%(cbartitle, xrVar.units),
             )
        cb.ax.tick_params(labelsize=11)
        #cb.ax.set_yticklabels(['{:.0f}'.format(x) for x in clevs])#, fontsize=16, weight='bold')
        if len(clevs) > 15:
            clevs2 = findNiceContours(clevs,nlevs = 10, rmClev=0.,sym=True)
            cb.set_ticks(clevs2)
            cb.set_ticklabels(clevs2)
            #cb.ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
            cb.ax.set_xticklabels(["{:.0f}".format(i) for i in clevs2]) # set ticks of your format


        
    return
```

```python
def pltfld(DV, titled, fname=None):
    
    cbartitle = DV.long_name
    
    if DV.min().values == DV.max().values:
        print('constant field skipping plot ')
    else:
        dlev_rng = {'CDNUMC':np.array([0.,3.e11])/2.,'FSNT':np.array([-45.,45.]),
                   'TGCLDLWP':np.array([-50.,50.]),'PRECL':np.array([-1.,1.]),
                    #'TGCLDLWP':np.array([-80.,80.]),'PRECL':np.array([-1.,1.]),
                    'PRECT':np.array([-5.,5.]),'SWCF':np.array([-45.,45.]),
                    'CLDLOW':np.array([-10.,10.]),'XXX':np.array([-45.,45.]),
                    'CLDTOT':np.array([-10.,10.]),'FSNTC':np.array([-20.,20.]),

                   }
        if DV.name in dlev_rng:
            dlevs = findNiceContours(dlev_rng[DV.name],nlevs = 15,rmClev=0.,sym=True)
        else:
            dlevs = findNiceContours(np.array([DV.min().values,DV.max().values]),nlevs = 15, rmClev=0.,sym=True)
        #dlevs = [-5.,-2.,-1.,-0.5,-0.2,-0.1,0.1,0.2,0.5,1.,2.,5.]
        #print('xxx',dlevs)
        dmap = diverge_map()

        plconf = '3-1x1'
        #plconf = '1x3'
        # good setup for 1 row of 3 columns
        # good setup for 3 rows of 1 columns
        if plconf == '3-1x1':
            fig, axes = setfig3b1x1()
            xr_llhplot2(DV, fig=fig, ax=axes,clevs=dlevs,cmap=dmap,title=titled, cbartitle=cbartitle)
            pltllbox([-150.,-110.],[0.,30.],ax=axes)
            pltllbox([-110.,-70.],[-30.,0.],ax=axes)
            pltllbox([-25.,15.],[-30.,0.],ax=axes)
            if fname is not None:
                print('fname',fname)
                plt.savefig(fname,dpi=300,transparent=True)
            plt.show()


```

```python


regtag = ""
weights = None

Varlist = np.array(['RESTOM','FLNT','FSNT','TS','TMQ','PRECT','AEROD_v','CLDLOW','CLDTOT','LWCF','SWCF','TGCLDIWP','TGCLDLWP',
                    'SHFLX','LHFLX','PBLH','PCONVT','PRECC','PRECS'])
Varlist = np.array(['FLNT','FSNT','TS','PRECC','PRECL','CLDLOW','CLDTOT','LWCF','SWCF','TGCLDIWP','TGCLDLWP',
                    'SHFLX','LHFLX','PBLH','PRECSC','PRECSL'])
#Varlist = np.array(['TS','TMQ','PRECT'])
#Varlist = np.array(['RESTOM','LWCF','SWCF','FLNT','FSNT'])
#Varlist = np.array(['AEROD_v'])
Varlist = np.array(['FSNT','FSNTC','FLNT','FLNTC'])
Varlist = np.array(['SWCF','TGCLDLWP','CLDLOW'])
#Varlist = np.array(['FSUTOA','FSNT','FSNTC'])
Varlist = np.array(['PRECC'])
Varlist = np.array(['SWCF','CLDTOT','TGCLDLWP','FSNT','FSNTC','AODVIS'])
Varlist = np.array(['TGCLDLWP'])
Varlist = np.array(['CLDTOT','TGCLDLWP','FSNT','FSNTC'])
#Varlist = np.array(['CLDTOT'])



# specify regions (assume lon always specified as west, then east limit)
xreg = np.array([[-150.,-110.],[-110,-70],[-25.,15.],[170.,-120.],[-170.,-90.]])%360.
yreg = np.array([[0.,30.],     [-30.,0.], [-30.,0.], [30.,50.],   [-50.,-30.] ])
namereg = ['NEP','SEP','SEA','NP','SP']
#xreg = [[0.,360.]]
#yreg = [[-90.,91.]]

case_start1 = "/home/jupyter-haruki/work/CESM_MCB/MCB_R1R2R3_CN375cm/MCB_R1R2R3_CN375cm.cam.h0." 
case_start1 = "/home/jupyter-haruki/work/CESM_MCB/MCB_R1R2R3_CN600cm/MCB_R1R2R3_CN600cm.cam.h0."
case_start1 = "/e3sm_prod/phil/climo/cesm/MCB_R1R2R3_CN600cm/fv192x288/MCB_R1R2R3_CN600cm.cam.h0.1-10."
case_end1 = ".nc"
pref1='CESM_CN600'
fstring1 ='%s%s%s'

case_start2 = "/home/jupyter-haruki/work/CESM_MCB/Fixed_SST/Fixed_SST.cam.h0."
case_start2 = "/e3sm_prod/phil/climo/cesm/Fixed_SST/fv192x288/Fixed_SST.cam.h0.1-20."
case_end2 = ".nc"
pref2='CESMcontrol'
fstring2 ='%s%s%s' 

if False:

    case_start1 = "/e3sm_prod/PJR/haruki_workdir/E3SM_MCB/F2010.E1_R1-3_C600_remapped/20221018.v2.LR.F2010.E1_R1-3_CDNC600.eam.h0.y1-5.FORCING.nc"
    case_end1 = ""
    pref1='E3SM_CN600'
    fstring1 ='%s%.0s%.0s' 

    case_start2 = "/e3sm_prod/PJR/haruki_workdir/E3SM_MCB/F2010.E1_R1-3_CNTL_remapped/20220930.v2.LR.F2010.E1_CNTL.eam.h0.y1-14.FORCING.nc"
    case_end2 = ""
    fstring2 ='%s%.0s%.0s' 
    pref2='E3SMcontrol'

if False: # E3SM CDNC2000 -1.8W/m2 CI
    #case_start1 = "/scratch2/PJR/haruki_workdir/E3SM_MCB/F2010.E1_R1-3_C600_remapped/20221018.v2.LR.F2010.E1_R1-3_CDNC600.eam.h0.1-11"
    case_start1 = "/e3sm_prod/phil/climo/e3sm/20221123.v2.LR.F2010.E1_R1-3_CDNC2000/fv192x288/20221123.v2.LR.F2010.E1_R1-3_CDNC2000.eam.h0.1-11."
    case_end1 = ".nc"
    pref1='E3SM_CN2000'
    fstring1 ='%s%.0s%.0s' 
    fstring1 ='%s%s%s' 

    #case_start2 = "/scratch2/PJR/haruki_workdir/E3SM_MCB/F2010.E1_R1-3_CNTL_remapped/20220930.v2.LR.F2010.E1_CNTL.eam.h0.1-14"
    case_start2 = "/e3sm_prod/phil/climo/e3sm/20220930.v2.LR.F2010.E1_CNTL/fv192x288/20220930.v2.LR.F2010.E1_CNTL.eam.h0.1-14."
    case_end2 = ".nc"
    fstring2 ='%s%.0s%.0s' 
    fstring2 ='%s%s%s' 
    pref2='E3SMcontrol'
    
if True: # E3SM 50Tg/yr
    #case_start1 = "/scratch2/PJR/haruki_workdir/E3SM_MCB/F2010.E1_R1-3_C600_remapped/20221018.v2.LR.F2010.E1_R1-3_CDNC600.eam.h0.1-11"
    case_start1 = "/e3sm_prod/phil/climo/e3sm/20230426.v2.LR.F2010.MCB-SSLT-EM.R1-3.test01/20230426.v2.LR.F2010.MCB-SSLT-EM.R1-3.test01_ANN_000101_001112_climo_fv192x288.nc"
    case_end1 = ""
    pref1='E3SM_50Tgpy'
    fstring1 ='%s%.0s%.0s' 
    #fstring1 ='%s%s%s' 

    #case_start2 = "/scratch2/PJR/haruki_workdir/E3SM_MCB/F2010.E1_R1-3_CNTL_remapped/20220930.v2.LR.F2010.E1_CNTL.eam.h0.1-14"
    case_start2 = "/e3sm_prod/phil/climo/e3sm/20220930.v2.LR.F2010.E1_CNTL/fv192x288/20220930.v2.LR.F2010.E1_CNTL.eam.h0.1-14."
    case_end2 = ".nc"
    fstring2 ='%s%.0s%.0s' 
    fstring2 ='%s%s%s' 
    pref2='E3SMcontrol'

if True: # CESM 7.5 Tg/yr
    case_start1 = "/e3sm_prod/phil/climo/cesm/F2010climo.ss_NEP_SEP_SEA.1.5Tg/fv192x288/F2010climo.ss_NEP_SEP_SEA.1.5Tg.cam.h0."
    case_end1 = ".1-25.nc"
    pref1='CESM_7.5Tgpyr'
    fstring1 ='%s%.0s%.0s' 
    fstring1 ='%s%s%s' 

    case_start2 = "/e3sm_prod/phil/climo/cesm/Fixed_SST/fv192x288/Fixed_SST.cam.h0.1-20."
    case_end2 = ".nc"
    fstring2 ='%s%.0s%.0s' 
    fstring2 ='%s%s%s' 
    pref2='CESMcontrol'



Varname='<Varname>'
ind1 = fstring1 % (case_start1,Varname,case_end1)
print('example string used for file open',ind1)

if False:
    indlf = fstring1 % (case_start1,'LANDFRAC',case_end1)
    print('indlf',indlf)
    DSLF = xr.open_mfdataset(indlf)
    lf = xr_getvar('LANDFRAC',DSLF).squeeze()
    #indof = fstring1 % (case_start1,'OCNFRAC',case_end1)
    #DSOF = xr.open_mfdataset(indof)
    #of = xr_getvar('OCNFRAC',DSOF).squeeze()
    of = 1.-lf # make ocean fraction the complement of lf (ignore ice)
    indif = fstring1 % (case_start1,'ICEFRAC',case_end1)
    DSIF = xr.open_mfdataset(indif)
    ifr = xr_getvar('ICEFRAC',DSIF).squeeze()


for Varname in Varlist:
    print()
    print('-------------------------------'+Varname)    
    ind1 = fstring1 % (case_start1,Varname,case_end1)
    print('opening',ind1)
    DS1 = xr.open_mfdataset(ind1)
    #print('xxx',DS1.time)
    DS1 = center_time(DS1)
    Var1 = xr_getvar(Varname,DS1)
    V1 = Var1.mean(dim='time',keep_attrs=True)

    ind2 = fstring2 % (case_start2,Varname,case_end2)
    print('opening ind2',ind2)
    #DS2 = xr.open_mfdataset(ind2)
    DS2 = xr.open_mfdataset(ind2)

    #DS2 = center_time(DS2)
    Var2 = xr_getvar(Varname,DS2)
    #print('yyy',Var2)
    #print('yy2',Var2.lat)
    #print('yy3',Var2.time)
    V2 = Var2.mean(dim='time',keep_attrs=True)

    DV = V1-V2
    print('DV range', DV.min().values, DV.max().values)
    if 'area' in DS1:
        area = DS1['area']
    elif 'area' in DS2:
        area = DS2['area']
    else:
        print('calculating weights')
        lat = Var1['lat'].values
        lon = Var1['lon'].values
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
    #sDVA = ' (%5.2f)' % DVA
    sDVA = ' ({:5.2f}{:s})'.format(DVA.values, DVA.units)

    print('area avgs '+pref1+' %5.2f' % (V1A.values),' '+pref2+' %5.2f' % (V2A.values),' Delta %5.2f' % (DVA.values))

    if V1.min().values == V1.max().values:
        print('constant field skipping plot ')
    else:
        clev_rng = {'CDNUMC':np.array([0.,3.e11]),'FSNT':np.array([40.,360]),
                    'TGCLDLWP':np.array([0.,280.]),'PRECL':np.array([0.,10.]),
                    'PRECC':np.array([0.,16.]),'SWCF':np.array([-140.,0.]),
                    'CLDLOW':np.array([0.,90.]),'AODVIS':np.array([0,0.2]),
                    'XXX':np.array([-45.,45.]),
                   }
        dlev_rng = {'CDNUMC':np.array([0.,3.e11])/2.,'FSNT':np.array([-45.,45.]),
                   'TGCLDLWP':np.array([-50.,50.]),'PRECL':np.array([-1.,1.]),
                    'PRECC':np.array([-1.,1.]),'SWCF':np.array([-45.,45.]),
                    'CLDLOW':np.array([-10.,10.]),'AODVIS':np.array([-.9,.9]),
                    'CLDTOT':np.array([-10.,10.]),'FSNTC':np.array([-20.,20.])
                   }
        if Varname in clev_rng:
            clevs = findNiceContours(clev_rng[Varname],nlevs = 12)
        else:
            clevs = findNiceContours(np.array([V1.values,V2.values]),nlevs = 12)
        if Varname in dlev_rng:
            dlevs = findNiceContours(dlev_rng[Varname],nlevs = 15,rmClev=0.,sym=True)
        else:
            dlevs = findNiceContours(np.array([DV.min().values,DV.max().values]),nlevs = 15, rmClev=0.,sym=True)
        #dlevs = [-5.,-2.,-1.,-0.5,-0.2,-0.1,0.1,0.2,0.5,1.,2.,5.]
        #print('xxx',dlevs)
        dmap = diverge_map()
        
        if Varname == 'CLDTOT':
            DV.attrs['long_name'] = "Cloud Cover"

        plconf = '3-1x1'
        #plconf = '1x3'

        # good setup for 1 row of 3 columns
        if plconf == '1x3':
            fig, axes = plt.subplots(ncols=3
                                     ,gridspec_kw={'width_ratios': [1, 1, 1]}
                                     ,subplot_kw={'projection': ccrs.Mollweide()}
                                     ,figsize=(16,5)
                                    )

            xr_llhplot(V1, ax=axes[0],clevs=clevs,title=pref1+sV1A)
            xr_llhplot(V2, ax=axes[1],clevs=clevs,ylabels=False,title=pref2+sV2A)
            xr_llhplot(DV, ax=axes[2],clevs=dlevs,cmap=dmap,title=pref1+'-'+pref2+sDVA)
            #plt.savefig(pref1+'_'+Varname+'.pdf',format='pdf',dpi=300)
            plt.show()
            
        # good setup for 3 rows of 1 columns
        if plconf == '3-1x1':
            
            if False:
                fig, axes = setfig3b1x1()
                #print('V1XXX',V1)
                xr_llhplot(V1, ax=axes,clevs=clevs,title=pref1+sV1A)
                pltllbox([-150.,-110.],[0.,30.])
                pltllbox([-110.,-70.],[-30.,0.])
                pltllbox([-25.,15.],[-30.,0.])
                plt.savefig(pref1+'_'+Varname+'.pdf',format='pdf',dpi=300)
                plt.show()

                fig, axes = setfig3b1x1()
                xr_llhplot(V2, ax=axes,clevs=clevs,ylabels=False,title=pref2+sV2A)
                pltllbox([-150.,-110.],[0.,30.])
                pltllbox([-110.,-70.],[-30.,0.])
                pltllbox([-25.,15.],[-30.,0.])
                plt.savefig(pref2+'_'+Varname+'.pdf',format='pdf',dpi=300)
                plt.show()

            if False:
                fig, axes = setfig3b1x1()
                xr_llhplot(DV, ax=axes,clevs=dlevs,cmap=dmap,title=pref1+'-'+pref2+sDVA)
                pltllbox([-150.,-110.],[0.,30.])
                pltllbox([-110.,-70.],[-30.,0.])
                pltllbox([-25.,15.],[-30.,0.])
                plt.savefig(pref1+'_'+Varname+'-D.pdf',format='pdf',dpi=300)
                plt.show()

            pref_fn = pref1
            difftitle = ''
            fname = pref_fn+'_'+DV.name+'-D.pdf'
            pltfld(DV, difftitle+sDVA,fname)

        
    print('field processing complete')

```

```python
print('rest of script not needed')
1./0.
```

```python


def getvarDSM(Varname,fstring,case_start,case_end):
    """getvar DSM
       get variable from file specifying the formatting
    """
    ind = fstring % (case_start,Varname,case_end)
    print('getvarDSM opening',ind)
    DS = xr.open_mfdataset(ind)
    DS.coords['lon'] = (DS.coords['lon'] + 180) % 360 - 180
    DS = DS.sortby(DS.lon)
    Var = xr_getvar(Varname,DS)
    VM = Var.mean(dim='time',keep_attrs=True)
    return VM;

def derfld(VN, fstring1, case_start1, case_end1):
    """ calculate some derived fields
    ICNUMLIQ850 is in-cloud number conc at 850 lev hybrid surface
    ICNUMLIQPBLT is in-cloud number conc near PBL top
    PRECT is total precip
    """
    if VN == 'ICNUMLIQ850':
        
        mylev=850.
        FREQL = getvarDSM('FREQL', fstring1, case_start1, case_end1).sel(lev=mylev,method='nearest')
        CLOUD = getvarDSM('CLOUD', fstring1, case_start1, case_end1).sel(lev=mylev,method='nearest')/100.

        NUMLIQ = getvarDSM('NUMLIQ', fstring1, case_start1, case_end1).sel(lev=mylev,method='nearest')*1.e-6
        ICNUMLIQ = NUMLIQ/(CLOUD+1.e-2)
        ICNUMLIQ = ICNUMLIQ.rename('ICNUMLIQ')
        ICNUMLIQ.attrs['long_name'] = 'approx in-cl number @850hPa hybsfc'
        ICNUMLIQ.attrs['units'] = '#/cc'
        return ICNUMLIQ
    elif VN == 'PRECT':
        print("PRECT")
        PRECL = getvarDSM('PRECL', fstring1, case_start1, case_end1)
        PRECC = getvarDSM('PRECC', fstring1, case_start1, case_end1)
        PRECT = PRECL+PRECC
        PRECT = PRECT.rename('PRECT')
        PRECT.attrs['long_name'] = 'total precipitation (liq + ice)'
        #PRECT.attrs['units'] = '#/cc'
        return PRECT
    elif VN == 'ICNUMLIQPBLT':
        PBLH = getvarDSM('PBLH', fstring1, case_start1, case_end1)
        Z3 = getvarDSM('Z3', fstring1, case_start1, case_end1)
        #P3 = getvarDSM('P3', fstring1, case_start1, case_end1)
        PHIS = getvarDSM('PHIS', fstring1, case_start1, case_end1)
        # find the 1D array of indices closest to the PBLH
        # I think Z3 is height above sea-level, and PBLH is height above surface
        Z3D = Z3 - PHIS/9.8 - PBLH
        llmin1 =  np.abs(Z3D.values).argmin(axis=0)
        lind = llmin1.flatten()
        indi = np.arange(0,len(lind))
        ni,nj,nk = np.shape(Z3.values)
        NUMLIQ = getvarDSM('NUMLIQ', fstring1, case_start1, case_end1)*1.e-6
        CLOUD = getvarDSM('CLOUD', fstring1, case_start1, case_end1)/100.
        ICNUMLIQ = NUMLIQ/(CLOUD+1.e-2)
        ICNUMLIQ = ICNUMLIQ.rename('ICNUMLIQ')
        ICNUMLIQ.attrs['long_name'] = 'approx in-cloud number conc'
        ICNUMLIQ.attrs['units'] = '#/cc'
        VAR = ICNUMLIQ.copy()
        # now extract data for each column at level lind
        data = VAR.values
        datas = data.reshape((ni,nj*nk))
        dataz = datas[lind,indi]
        datazr = dataz.reshape((nj,nk))
        VARs = VAR.sel(lev=1000.,method='nearest')
        VARs.attrs['long_name'] = VARs.attrs['long_name']+' near PBLH'
        VARs.data = datazr
        return VARs
    else:
        1./0.
        
# lines to support debugging the above functions

# sample multilevel field on hybrid levels
#fig, axes = setfig3b1x1()

#PRECT = derfld('PRECT',fstring1, case_start1, case_end1)
#PRECTA = PRECT.weighted(weights).mean()
#print('PRECTA',PRECTA.values)
#xr_llhplot(PRECT, ax=axes)#,clevs=clevs,title=pref1+sV1A)
if False:
    ICNUMLIQ1 = derfld('ICNUMLIQPBLT',fstring1, case_start1, case_end1)
    xr_llhplot(ICNUMLIQ1, ax=axes)#,clevs=clevs,title=pref1+sV1A)
    pltllbox([-150.,-110.],[0.,30.])
    pltllbox([-110.,-70.],[-30.,0.])
    pltllbox([-25.,15.],[-30.,0.])
    plt.show()
    
if False:
    ICNUMLIQ2 = derfld('ICNUMLIQPBLT',fstring2, case_start2, case_end2)
    xr_llhplot(ICNUMLIQ2, ax=axes)#,clevs=clevs,title=pref1+sV1A)
    #plt.savefig(pref1+'_'+Varname+'.pdf',format='pdf',dpi=300)
    pltllbox([-150.,-110.],[0.,30.])
    pltllbox([-110.,-70.],[-30.,0.])
    pltllbox([-25.,15.],[-30.,0.])
    plt.show()

#PRECT = derfld('PRECT',fstring1, case_start1, case_end1)

```

```python
# sample multilevel field at level nearest PBLH
fig, axes = setfig3b1x1()

PBLH = getvarDSM('PBLH', fstring1, case_start1, case_end1)
Z3 = getvarDSM('Z3', fstring1, case_start1, case_end1)
#P3 = getvarDSM('P3', fstring1, case_start1, case_end1)
PHIS = getvarDSM('PHIS', fstring1, case_start1, case_end1)

# find the 1D array of indices closest to the PBLH
# I think Z3 is height above sea-level, and PBLH is height above surface
Z3D = Z3 - PHIS/9.8 - PBLH
llmin1 =  np.abs(Z3D.values).argmin(axis=0)
#print('llmin',np.shape(llmin1))
lind = llmin1.flatten()
#print('lind',np.shape(lind),lind[0:3])
indi = np.arange(0,len(lind))
#print('indi',np.shape(indi))
ni,nj,nk = np.shape(Z3.values)
#print('ni,nj,nk',ni,nj,nk)

NUMLIQ = getvarDSM('NUMLIQ', fstring1, case_start1, case_end1)*1.e-6
CLOUD = getvarDSM('CLOUD', fstring1, case_start1, case_end1)/100.
ICNUMLIQ = NUMLIQ/(CLOUD+1.e-2)
ICNUMLIQ = ICNUMLIQ.rename('ICNUMLIQ')
ICNUMLIQ.attrs['long_name'] = 'approx in-cloud number conc'
ICNUMLIQ.attrs['units'] = '#/cc'

VAR = ICNUMLIQ.copy()
# now extract data for each column
data = VAR.values
datas = data.reshape((ni,nj*nk))
#print('datas',np.shape(datas))
dataz = datas[lind,indi]
#print('dataz',np.shape(dataz))
datazr = dataz.reshape((nj,nk))
#print('datazr shape ij',datazr.shape, datazr[i,j])
#datam = data[llmin1,:,:]
#print('datam',datam.shape)
VARs = VAR.sel(lev=1000.,method='nearest')
VARs.attrs['long_name'] = VARs.attrs['long_name']+' near PBLH'
VARs.data = datazr

xr_llhplot(VARs, ax=axes)#,clevs=clevs,title=pref1+sV1A)
pltllbox([-150.,-110.],[0.,30.])
pltllbox([-110.,-70.],[-30.,0.])
pltllbox([-25.,15.],[-30.,0.])

#plt.savefig(pref1+'_'+Varname+'.pdf',format='pdf',dpi=300)
plt.show()

```

```python

lon = xr_getvar('lon',DS1)
lat = xr_getvar('lat',DS1)
indof = fstring1 % (case_start1,'OCNFRAC',case_end1)
DSOF = xr.open_mfdataset(indof)
of = xr_getvar('OCNFRAC',DSOF).squeeze()
print('of',of)
DVcum = 0.
#DV = V1
print('analyzing variable DV.name=',DV.name)
for i in range(len(xreg)):
    xri = xreg[i]
    yri = yreg[i]
    print('namereg, xreg',namereg[i],xri)
    #print('lon',lon.values)
    #R2 = DV.sel(lon=slice(220,250),lat=slice(15,35))
    #%print('R2',R2)
    # allow longitudes to cross greenwich
    print('DV range', DV.min().values, DV.max().values)

    if xri[1]>xri[0]: 
        DVI = DV.sel(lon=(lon < xri[1]) & (lon >= xri[0]),
                     lat=(lat >=yri[0]) & (lat < yri[1])
                    ).load()
        areai = area.sel(lon=(lon < xri[1]) & (lon >= xri[0]),
                     lat=(lat >=yri[0]) & (lat < yri[1])
                    ).load()
    else:
        DVI = DV.sel(lon=(lon >= xri[0]) | (lon < xri[1]),
                     lat=(lat >=yri[0]) & (lat < yri[1])
                    ).load()
        areai = area.sel(lon=(lon >= xri[0]) | (lon < xri[1]),
                     lat=(lat >=yri[0]) & (lat < yri[1])
                    ).load()

    #print('DVI',DVI)
    #print('DVI range', DVI.min().values, DVI.max().values)

    DVIG = (DVI*areai).sum()/area.sum()
    DVcum = DVcum + DVIG.values
    print('i, areai, area, fglob, focean ', i, areai.sum().values, area.sum().values,
         (areai*of).sum().values/area.sum().values,
         (areai*of).sum().values/(area*of).sum().values)
    #print('i, DVI.sum, DVI*areai,', DVI.sum().values, (DVI*areai).sum().values)
    print('i, DVIG and cum',i, DVIG.values, DVcum)
    print(' ')
    #break
```

```python
xmesh, ymesh = np.meshgrid(lon, lat)
indlf = fstring1 % (case_start1,'LANDFRAC',case_end1)
DSLF = xr.open_mfdataset(indlf)
lf = xr_getvar('LANDFRAC',DSLF).squeeze()
indof = fstring1 % (case_start1,'OCNFRAC',case_end1)
DSOF = xr.open_mfdataset(indof)
of = xr_getvar('OCNFRAC',DSOF).squeeze()
indif = fstring1 % (case_start1,'ICEFRAC',case_end1)
DSIF = xr.open_mfdataset(indif)
ifr = xr_getvar('ICEFRAC',DSIF).squeeze()

DVcum = 0.
i = 2
if True:
#for i in range(len(xreg)):
    xri = xreg[i].copy()
    yri = yreg[i]
    # find indices in region of interest
    if xri[1]>=xri[0]:
        inds = ((xmesh>xri[0]) & (xmesh <= xri[1]) & (ymesh > yri[0]) & (ymesh <= yri[1]))
    else:
        inds = (((xmesh >= xri[0]) | (xmesh < xri[1])) &
                (ymesh >=yri[0]) & (ymesh < yri[1]))
    # remove points with land from the mask
    inds = inds & (lf < 0.02)
    #print('DISABLED land mask')
    #DVx = DV.where((xmesh>xri[0]) & (xmesh <= xri[1]) & (ymesh > yri[0]) & (ymesh <= yri[1]))
    DVx = DV.where(inds)
    areax = area.where(inds)
    DVIG = (DVx*areax).sum()/area.sum()
    DVIL = (DVx*areax).sum()/areax.sum()
    DVcum = DVcum + DVIG.values
    fglob = (areax.sum()/area.sum()).values
    focn = (areax.sum()/(area*of).sum()).values
    print('areax',areax.sum().values)
    print('area',area.sum().values)
    print('areao',(area*of).sum().values)
    print('areai',(area*ifr).sum().values)
    print('areal',(area*lf).sum().values)
    print(i," DVIG = {:.2f}".format(DVIG.values), 
          "DVcum = {:.2f}".format(DVcum),
          "DVIL = {:.2f}".format(DVIL.values), 
          "fglob = {:.2%}".format(fglob),
          "focn = {:.2%}".format(focn)
         )
    if xri[1] < xri[0]:
        xri[1] += 360.
    regcx = [xri[0],xri[1],xri[1],xri[0],xri[0]]
    regcy = [yri[0],yri[0],yri[1],yri[1],yri[0]]
    #fig, axes = setfig3b1x1()
    map_proj = ccrs.Mollweide()
    map_proj._threshold /= 100.
    fig, axes = plt.subplots(ncols=1,
                             gridspec_kw={'width_ratios': [1]},
                             subplot_kw={'projection':map_proj},
                             figsize=(6,3),
                            )
    fig.set_dpi(300.0)
    xr_llhplot(DVx, ax=axes)
    plt.plot(regcx,regcy,color='green',transform=ccrs.PlateCarree())
    plt.show()
    #break
```

```python
xmesh, ymesh = np.meshgrid(lon, lat)
indlf = fstring1 % (case_start1,'LANDFRAC',case_end1)
DSLF = xr.open_mfdataset(indlf)

# add new lon coordinate that runs from -180 to 180. for use when region crosses the gm
lf = xr_getvar('LANDFRAC',DSLF).squeeze()
lf = lf.assign_coords(longm=(((lf.lon +180.) % 360) - 180))

area = xr_getvar('area',DSLF).squeeze()
area = area.assign_coords(longm=(((area.lon + 180) % 360) - 180))

DV = DV.assign_coords(longm=(((DV.lon + 180) % 360) - 180))
#print('DV',DV)

DVcum = 0.
i = 2
if True:
#for i in range(len(xreg)):
    print(' ')
    print('xreg i ', xreg[i])
    xri = xreg[i].copy()
    xrigm = ((xreg[i]+180)%360)-180.
    print('xrigm',xrigm)
    yri = yreg[i]
    # find indices in region of interest
    indlat = np.where((DV.lat.values>yri[0]) & (DV.lat.values <= yri[1]))[0]
    print('indlat',indlat)
    if xri[1]>=xri[0]:
        indlon = np.where((DV.lon.values>xri[0]) & (DV.lon.values<= xri[1]))[0]
        DVxx = DV.isel(lon=indlon,lat=indlat).copy()
        print('DVxx',DVxx)
        areaxx = area.isel(lon=indlon,lat=indlat)
        lfxx = lf.isel(lon=indlon,lat=indlat)
    else:
        indlon = np.where((DV.longm.values>xrigm[0]) & (DV.longm.values<= xrigm[1]))[0]
        DVxx = DV.isel(lon=indlon,lat=indlat).copy()
        print('DVxx',DVxx)
        areaxx = area.isel(lon=indlon,lat=indlat)
        lfxx = lf.isel(lon=indlon,lat=indlat)
    print('indlon',indlon)
    print('DVxx.longm',DVxx.longm.values)
    print('DVnew', DVxx.sortby("longm"))

    """
    # remove points with land from the mask
    inds = inds & (lf < 0.02)
    #print('DISABLED land mask')
    #DVx = DV.where((xmesh>xri[0]) & (xmesh <= xri[1]) & (ymesh > yri[0]) & (ymesh <= yri[1]))
    DVx = DV.where(inds)
    areax = area.where(inds)
    DVIG = (DVx*areax).sum()/area.sum()
    DVIL = (DVx*areax).sum()/areax.sum()
    DVcum = DVcum + DVIG.values
    fglob = (areax.sum()/area.sum()).values
    focn = (areax.sum()/(area*of).sum()).values
    print('areax',areax.sum().values)
    print('area',area.sum().values)
    print('areao',(area*of).sum().values)
    print('areai',(area*ifr).sum().values)
    print('areal',(area*lf).sum().values)
    print(i," DVIG = {:.2f}".format(DVIG.values), 
          "DVcum = {:.2f}".format(DVcum),
          "DVIL = {:.2f}".format(DVIL.values), 
          "fglob = {:.2%}".format(fglob),
          "focn = {:.2%}".format(focn)
         )
    if xri[1] < xri[0]:
        xri[1] += 360.
    regcx = [xri[0],xri[1],xri[1],xri[0],xri[0]]
    regcy = [yri[0],yri[0],yri[1],yri[1],yri[0]]
    #fig, axes = setfig3b1x1()
    map_proj = ccrs.Mollweide()
    map_proj._threshold /= 100.
    print('xri',xri)
    print('2nd xreg i', xreg[i])
    fig, axes = plt.subplots(ncols=1,
                             gridspec_kw={'width_ratios': [1]},
                             subplot_kw={'projection':map_proj},
                             figsize=(6,3),
                            )
    fig.set_dpi(300.0)
    xr_llhplot(DVx, ax=axes)
    plt.plot(regcx,regcy,color='green',transform=ccrs.PlateCarree())
    plt.show()
    #break
    """
```

```python

if True:
    map_proj = ccrs.Mollweide()
    #map_proj = ccrs.PlateCarree()

    map_proj._threshold /= 100.
    fig, axes = plt.subplots(ncols=1,
                             gridspec_kw={'width_ratios': [1]},
                             subplot_kw={'projection':map_proj},
                             figsize=(6,3),
                            )
    fig.set_dpi(300.0)
    #xr_llhplot(DVxx, ax=axes)
    DVnew = DVxx.sortby("longm")
    #DVnew = DVxx.sortby("lon")
    xmesh, ymesh = np.meshgrid(DVnew.longm.values, DVnew.lat.values)
    axes.coastlines()
    axes.contourf(xmesh, ymesh, DVnew.values, transform=ccrs.PlateCarree())
    ylabels = True 
    gl = axes.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,
                      linewidth=2, color='gray', alpha=0.5)
    gl.left_labels=ylabels
    gl.right_labels=ylabels
    gl.bottom_labels= True
    #axes.set_extent([-180., 180., -90., 90.], crs=ccrs.PlateCarree())
    #axes.set_extent([-180., 180., -90., 90.], crs=map_proj)
    axes.set_global()
    plt.show()
```

```python
lonsub = DS1['lon'+regtag]#.where(pmask)#.isel(time=0)
latsub = DS1['lat'+regtag]#.where(pmask)#.isel(time=0)
lat = latsub.values
print('subreg',lonsub.min().values,lonsub.max().values, latsub.min().values,latsub.max().values)
#print('subreg size',lonsub.shape)
print('shape and size of variables',lonsub.shape, lonsub.size,' number of unmasked cells ',np.count_nonzero(lonsub.notnull().values))

if 'area' in DS1:
    print('area present in DS1')
    area = xr_getvar('area',DS1,regtag=regtag)#.where(pmask)
    wtsh = area.fillna(0).values
else:
    print('area not found')
    wtsh = make_fvarea(lonsub.values,latsub.values)

#print('weights shape',wtsh.shape,wtsh.sum())

Varlist = np.array(['T','AREL','CLOUD'])
#Varlist = np.array(['AREL'])

#print('wtsh',wtsh)

for Vname in Varlist:
    print()
    print('-------------------------------')
    ind1 = fstring1 % (case_start1,Vname,case_end1)
    #print('ind1 is ', ind1)
    DS1 = xr.open_mfdataset(ind1).transpose(...,'lev') # make sure lev is last index
    DS1 = center_time(DS1)
    ind2 = fstring2 % (case_start2,Vname,case_end2)  # make sure lev is last index
    DS2 = xr.open_mfdataset(ind2).transpose(...,'lev') 
    DS2 = center_time(DS2)
    
    V1 = xr_getvar(Vname,DS1,regtag=regtag)#.where(pmask)
    if V1.min().values == V1.max().values:
        print('constant field skipping plot ')
    else:
        V1 = xr_getvar(Vname, DS1, regtag).mean(dim='time',keep_attrs=True)#.where(pmask).squeeze()
        #print('V1xxx ', V1)
        #print('V1',V1.min().values,V1.max().values)
        #print('V1 shape, size, realsize', V1.shape, np.size(V1.values), np.count_nonzero(V1.notnull().values) )
        V2 = xr_getvar(Vname, DS2, regtag).mean(dim='time',keep_attrs=True)#.where(pmask).squeeze()
        DV = V1-V2
        print(Vname, V1.attrs['long_name'],'Range V1 and V2 ',V1.min().values, V1.max().values, V2.min().values, V2.max().values)

        # """

        DPOG1 = xr_getvar('DPOG',DS1,regtag=regtag).mean(dim='time',keep_attrs=True)
        weights1 = wtsh[...,np.newaxis]*DPOG1
        DPOG2 = xr_getvar('DPOG',DS2,regtag=regtag).mean(dim='time',keep_attrs=True)
        weights2 = wtsh[...,np.newaxis]*DPOG2
        V1A = V1.weighted(weights1).mean()
        V2A = V2.weighted(weights2).mean()
        print('Area/Mass weighted avg '+pref1+' %5.3f' % (V1A.values),' '+pref2+' %5.3f' % (V2A.values))
        
        # """
        
        data1 = V1.mean(dim='lon',keep_attrs=True).transpose().values
        #print('data1 shape',data1.shape)
        data2 = V2.mean(dim='lon',keep_attrs=True).transpose().values
        datad = data1-data2
        

        lev = V1['lev'].values
        #print('lev shape',lev.shape)
    
#        plotZMf(data1, lat_h, lev)
        fig, axes = plt.subplots(ncols=3
                                 ,gridspec_kw={'width_ratios': [1, 1, 1]}
#                                 ,subplot_kw={'projection': ccrs.PlateCarree()}
                                 ,figsize=(16,5)
                                )
        ytop = 20.
        plotZMf(data1, lat, lev,axesa=axes[0],plotOpt={'colorbar':"botnd",'units':V1.units,'ltitle':pref1,'ytop':ytop})
        plotZMf(data2, lat, lev,axesa=axes[1],plotOpt={'colorbar':"bot",'units':V2.units,'ltitle':pref2,'rtitle':V2.long_name,'ytop':ytop})
        dlevs = findNiceContours(np.array([datad.min(),datad.max()]),nlevs = 10, rmClev=0.,sym=True)
        dmap = diverge_map()
        plotZMf(datad, lat, lev,axesa=axes[2],plotOpt={'clevs':dlevs,'cmap':dmap,'colorbar':"bot",'units':V2.units,'ytop':ytop,'ltitle':pref1+'-'+pref2})

        #print('attribute check on xarray',hasattr(V1,'units'))
        #plt.savefig('test_'+Vname+'.pdf',format='pdf')
        #plt.tight_layout()
        plt.show()
      

```

```python
# next cells not used
help(xr_getvar)
1./0.
```

```python
Varlist = np.array(['T','Q','CLOUD','CLDLIQ','ICWMR','CLDICE','RELHUM','NUMICE','NUMLIQ','Mass_bc'])
#                    RESTOM','FLNT','FSNT','TS','TMQ','PRECT','AEROD_v','CLDLOW','CLDTOT','LWCF','SWCF','TGCLDIWP','TGCLDLWP','SHFLX','LHFLX','PBLH','PCONVT','PRECC','PRECS'])
#Varlist = np.array(['TS','TMQ','PRECT'])
#Varlist = np.array(['RESTOM','LWCF','SWCF','FLNT','FSNT'])

Varname = 'TS' # Varlist[0]

case_start1 = "/home/jupyter-haruki/work/CESM_MCB/Fixed_SST/"
case_end1 = "Fixed_SST.cam.h0."+Varname+".y1-19.nc"
ind1 = case_start1+case_end1
DS1 = xr.open_mfdataset(ind1)


Var1a = DS1[Varname]
Var1am = Var1a.mean(dim='time')
DS1 = center_time(DS1)
Var1 = DS1[Varname]
#print('xxx',Var1)

case_start2 = "/home/jupyter-haruki/work/CESM_MCB/MCB_R1R2R3_CN375cm/" 
case_end2 = "MCB_R1R2R3_CN375cm.cam.h0."+Varname+".y1-10.nc"
ind2 = case_start2+case_end2
DS2 = xr.open_mfdataset(ind2)
DS2 = center_time(DS2)
Var2 = DS2[Varname]
#print('yyy',Var2)

```
