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

```python
import sys
print(sys.version)
%matplotlib inline
from xhistogram.xarray import histogram
%run -i ~/Python/pjr3
```

```python
# process file holding data only over a region at every timestep of a model run
indir = '/global/cscratch1/sd/pjr/E3SMv2/v2.LR.histAMIP_e1/tests/M_1x2_ndays/run/v2.LR.histAMIP_e1.eam.h1.2015-07-*0.nc'
indir = '/global/cscratch1/sd/pjr/E3SMv2/v2.LR.histAMIP_x3/tests/M_1x2_ndays/run/v2.LR.histAMIP_x3.eam.h1.2015-07-*0.nc'
prefix = 'exp_'
indir = '/global/cscratch1/sd/pjr/E3SMv2/v2.LR.histAMIP_x4/tests/M_1x10_ndays/run/v2.LR.histAMIP_x4.eam.h1.2015-07-*0.nc'
prefix = 'con_'
indir = '/global/cscratch1/sd/pjr/E3SMv2/v2.LR.histAMIP_e4/tests/M_1x10_ndays/run/v2.LR.histAMIP_e4.eam.h1.2015-07-*0.nc'
regtag = '_190e_to_250e_0n_to_35n'
xr.set_options(keep_attrs=True)
DS0 = xr.open_mfdataset(indir)
DS = DS0
# reorder coords so ncol is alway first dim 
# to allow lat/lon interpolation across multiple dimensions
DS = DS0.transpose('ncol'+regtag,...) 
tlast = DS.time[-1]
print('tlast',tlast)

lon = DS['lon'+regtag].isel(time=0)#.squeeze()
print('lon',lon.shape,lon.min().values,lon.max().values)
lat = DS['lat'+regtag].isel(time=0)
print('lat',lat.shape,lat.min().values,lat.max().values)

# create a mask to isolate a region of interest
pmask = ((lon > 220) & (lon < 250) & (lat > 15))#[0] # select a subregion
#pmask = (lon > -999) # select all points
print('pmask',pmask)
lonsub = DS['lon'+regtag].isel(time=0).where(pmask)
latsub = DS['lat'+regtag].isel(time=0).where(pmask)
print('subreg',lonsub.min().values,lonsub.max().values, latsub.min().values,latsub.max().values)
print('subreg size',lonsub.shape)

PRECC = xr_getvar("PRECC", DS, regtag).where(pmask)
#PRECC[:,:] = 1.
#print('PRECC',PRECC)
print('PRECC',PRECC.min().values,PRECC.max().values)
print('PRECC shape, size, realsize', PRECC.shape, np.size(PRECC.values), np.count_nonzero(PRECC.notnull().values) )
xr_cshplot(PRECC.isel(time=0), lonsub, latsub)
plt.show()

PRECC.plot()
plt.show()
PRECS = (xr_getvar("PRECT", DS, regtag).where(pmask) - PRECC)
PRECS = PRECS.rename("PRECS")
PRECS.attrs['long_name']='stratiform precipitation'
PRECS.attrs['basename']='PRECS'
print('PRECS',PRECS.min().values,PRECS.max().values)
PRECS.plot()
pmax = np.max([PRECS.max().values, PRECC.max().values])
print('pmax',pmax)
plt.show()      
bins = np.linspace(0, 30, 31)
bins = np.array([0.,0.01,0.02,0.05,0.1,0.2,0.5,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,20.,30.])
bins = np.array([0.01,0.02,0.05,0.1,0.2,0.5,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,20.,30.])
bins = np.exp(findNiceContours(np.log(np.array([0.01,30.])),nlevs=30))
bins = np.insert(np.append(bins,pmax),0,0)
#bins = findNiceContours(np.array([-pmax*1.01,pmax*1.01]),nlevs=100)


print("bins",bins)
fig, axes = plt.subplots(figsize=(7,7))
#print('axes',axes)                      )

hc = histogram(PRECC, bins=[bins])
print('hc',hc)
print('hc.values',hc.values)
print('hc 1-val',hc.sel(PRECC_190e_to_250e_0n_to_35n_bin=[1.],method='nearest').values)
print('hc.sum.values',hc.sum().values)

# count the drizzling samples (precip less than 1 mm/day)
ndriz = hc.sel(PRECC_190e_to_250e_0n_to_35n_bin=slice(0.01,1.)).sum().values 
# total samples in the bins (there can be samples outside of the bins that aren't counted)
tcnt = hc.sum().values
print('ndriz',ndriz,tcnt,ndriz/tcnt)


#display(h)
hc[1:].plot(xscale='log')
hs = histogram(PRECS, bins=[bins])
print('hs',hs)
print('hs 1-val',hs.sel(PRECS_bin=[1.],method='nearest').values)
print('hs.sum.values',hs.sum().values)

#display(h)
hs[1:].plot(xscale='log',add_legend=True)
axes.legend(['convective','stratiform'])
plt.plot([1.,1.,1.],[0.,hc[1:].max().values,hs[1:].max().values])
plt.savefig(prefix+'1dprechisto'+'.jpg',format='jpg')
plt.show()
      
```

```python
# compare sums using original and histogram representations
print('PRECC.sum.values',PRECC.sum().values)
print('hc*bin.sum.values',(hc*hc.PRECC_190e_to_250e_0n_to_35n_bin).sum().values)
bl = bins[0:-1]
br = bins[1:]
print('lower bound on sum',(hc*bl).sum().values)
print('upper bound on sum',(hc*br).sum().values)

```

```python

histogram(PRECS, PRECC, bins=[bins[1:], bins[1:]]).plot(xscale='log',yscale='log')
plt.plot([0.01,30.],[0.01,30.],color='white')
plt.plot(np.array([0.01,30.]),np.array([0.01,30./10]),color='red')
```

```python
Varlist = np.array(['num_a3SFSBS','num_a3SFSBC','num_a3SFSIS','num_a3SFSIC',
                    'num_c3SFSBS','num_c3SFSBC','num_c3SFSIS','num_c3SFSIC'])
#Varlist = np.array(['num_a3SFSBS','num_a3SFSIS',
#                    'num_c3SFSBS','num_c3SFSIS'])
#Varlist = [Varlist[0]]
for Vname in Varlist:
    print()
    print('-------------------------------')
    V1 = -xr_getvar(Vname,DS,regtag=regtag).where(pmask)
    print(Vname, V1.attrs['long_name'],'Range ',V1.min().values, V1.max().values)
    print('size of Variable',V1.shape, V1.size,' number of unmasked cells ',np.count_nonzero(V1.notnull().values))
    print('Sum %5.3f' % (V1.sum().values*1.e-6), 'million, and avg =', V1.mean().values)
    if V1.min().values == V1.max().values:
        print('constant field skipping histogram ')
    else:
        fig, axes = plt.subplots(ncols=2
                        ,gridspec_kw={'width_ratios': [1, 1]}
                        ,figsize=(16,5)
                        )
        v1max = V1.max().values
        v1d4 = v1max*1.e-4

        binv1 = np.exp(findNiceContours(np.log(np.array([v1d4,v1max])),nlevs=40))
        binv1 = np.insert(np.append(binv1,v1max),0,0.)
        #print('variable bins ',binv1.shape,binv1)
        dbinv1 = np.diff(binv1)
        #print('dbinv1',dbinv1)

        hs2d = histogram(PRECS, V1, bins=[bins, binv1])
        #print('hs2d',hs2d)
        #print('hs prec bin ',hs2d.PRECS_bin, hs2d.PRECS_bin[0:3])
        #print('bin midpoints', (bins[1:4]+bins[0:3])*0.5)
        #print('number of hits in histogram hs.sum',hs2d.sum().values, 'number of values outside of histogram range',dsub.size-hs2d.sum().values)
        print('number of hits in lowest stratiform precip bin ',bins[0:2],hs2d[0,:].sum().values)
        #print('hs_sub.PRECS_bin',hs_sub.PRECS_bin.values)
        hs2d[1:,1:].plot(xscale='log',yscale='log',ax=axes[0])
        #print('hs2d shape',hs2d.shape)
        prshs = hs2d*hs2d.PRECS_bin[:]
        #print('prshs shape', prshs.shape)
        print('preshs sum (tot, then excluding bin0)',prshs.sum().values, prshs[1:,:].sum().values)
        #hsprecc = hs*
        hc2d = histogram(PRECC, V1, bins=[bins, binv1])
        hc2d[1:,1:].plot(xscale='log',yscale='log',ax=axes[1])
        #plt.plot([0.01,30.],[0.01,30.],color='white')
        #plt.plot(np.array([0.01,30.]),np.array([0.01,30./10]),color='red')
        plt.savefig(prefix+'2d_'+Vname+'_histo'+'.jpg',format='jpg')
        plt.tight_layout()
        plt.show()
        #fig, axes = plt.subplots(figsize=(7,7))
        #hs2div = hs2d.sum(dim='num_a3SFSBS_190e_to_250e_0n_to_35n_bin')
        #print('hs2div',hs2div)
        #hs2div[1:].plot(xscale='log',add_legend=True)
        #hs[1:].plot(xscale='log',add_legend=True)
        #axes.legend(['2div','1d'])
        #plt.show()
        #hs2dwf = hs2d*(binv1[1:]+binv1[0:-1])[np.newaxis,:]*0.5
        #hs2dwf[1:,1:].plot(xscale='log',yscale='log')
        fig, axes = plt.subplots(ncols=2
                                 ,gridspec_kw={'width_ratios': [1, 1]}
                                 ,figsize=(16,5)
                                )

        hs2dwf = hs2d*(binv1[1:]+binv1[0:-1])[np.newaxis,:]*0.5
        hs2dwf.attrs = V1.attrs
        hs2dwf.attrs['long_name'] = 'Count-weighted '+hs2dwf.attrs['long_name']
        hc2dwf = hc2d*(binv1[1:]+binv1[0:-1])[np.newaxis,:]*0.5
        hc2dwf.attrs = V1.attrs
        hc2dwf.attrs['long_name'] = 'Count-weighted '+hc2dwf.attrs['long_name']
        
        hs2dwf[1:,1:].plot(xscale='log',yscale='log',ax=axes[0])
        hc2dwf[1:,1:].plot(xscale='log',yscale='log',ax=axes[1])
        
        plt.tight_layout()
        plt.savefig(prefix+'2dwt_'+Vname+'_histo'+'.jpg',format='jpg')
        plt.show()
        # next calculations not currently used
        # weighted  histogram
        #whs = hs*dbins[:,np.newaxis]*dbinv1[np.newaxis,:]
        #print('whs',whs)
        # sub histogram
        #hs_sub = hs.sel(PRECS_bin=slice(0.,1.))
        #print('sub hist of PRECCS btw 0., 1. mm/d (limits sampled bin midpoint val)',hs_sub)

```

```python
# next cells not used
help(xr_getvar)
1./0.
```

```python
Varlist = np.array(['num_a3SFSBS','num_a3SFSBC','num_a3SFSIS','num_a3SFSIC',
                    'num_c3SFSBS','num_c3SFSBC','num_c3SFSIS','num_c3SFSIC'])
#Varlist = np.array(['num_a3SFSBS','num_a3SFSIS',
#                    'num_c3SFSBS','num_c3SFSIS'])
Varlist = [Varlist[0]]
for Vname in Varlist:
    print()
    print('-------------------------------')
    V1 = -xr_getvar(Vname,DS,regtag=regtag).where(pmask)
    print(Vname, V1.attrs['long_name'],'Range ',V1.min().values, V1.max().values)
    print('size of Variable',V1.shape, V1.size,' number of unmasked cells ',np.count_nonzero(V1.notnull().values))
    print('Sum %5.3f' % (V1.sum().values*1.e-6), 'million, and avg =', V1.mean().values)
    if V1.min().values == V1.max().values:
        print('constant field skipping histogram ')
    else:
        fig, axes = plt.subplots(ncols=2
                        ,gridspec_kw={'width_ratios': [1, 1]}
                        ,figsize=(16,5)
                        )
        v1max = V1.max().values
        v1d4 = v1max*1.e-4

        binv1 = np.exp(findNiceContours(np.log(np.array([v1d4,v1max])),nlevs=40))
        binv1 = np.insert(np.append(binv1,v1max),0,0.)
        #print('variable bins ',binv1.shape,binv1)
        dbinv1 = np.diff(binv1)
        #print('dbinv1',dbinv1)

        hs2d = histogram(PRECS, V1, bins=[bins, binv1])
        print('hs2d',hs2d)
        #print('hs prec bin ',hs2d.PRECS_bin, hs2d.PRECS_bin[0:3])
        #print('bin midpoints', (bins[1:4]+bins[0:3])*0.5)
        #print('number of hits in histogram hs.sum',hs2d.sum().values, 'number of values outside of histogram range',dsub.size-hs2d.sum().values)
        print('number of hits in lowest stratiform precip bin ',bins[0:2],hs2d[0,:].sum().values)
        #print('hs_sub.PRECS_bin',hs_sub.PRECS_bin.values)
        hs2d[1:,1:].plot(xscale='log',yscale='log',ax=axes[0])
        print('hs2d shape',hs2d.shape)
        prshs = hs2d*hs2d.PRECS_bin[:]
        print('prshs shape', prshs.shape)
        print('preshs sum (tot, then excluding bin0)',prshs.sum().values, prshs[1:,:].sum().values)
        #hsprecc = hs*
        hc2d = histogram(PRECC, V1, bins=[bins, binv1])
        hc2d[1:,1:].plot(xscale='log',yscale='log',ax=axes[1])
        #plt.plot([0.01,30.],[0.01,30.],color='white')
        #plt.plot(np.array([0.01,30.]),np.array([0.01,30./10]),color='red')
        plt.savefig(prefix+'2d_'+Vname+'_histo'+'.jpg',format='jpg')
        plt.show()
        #fig, axes = plt.subplots(figsize=(7,7))
        #hs2div = hs2d.sum(dim='num_a3SFSBS_190e_to_250e_0n_to_35n_bin')
        #print('hs2div',hs2div)
        #hs2div[1:].plot(xscale='log',add_legend=True)
        #hs[1:].plot(xscale='log',add_legend=True)
        #axes.legend(['2div','1d'])
        #plt.show()
        #hs2dwf = hs2d*(binv1[1:]+binv1[0:-1])[np.newaxis,:]*0.5
        #hs2dwf[1:,1:].plot(xscale='log',yscale='log')
        fig, axes = plt.subplots(ncols=2
                                 ,gridspec_kw={'width_ratios': [1, 1]}
                                 ,figsize=(16,5)
                                )

        hs2dwf = hs2d*(binv1[1:]+binv1[0:-1])[np.newaxis,:]*0.5
        #print('V1.units',V1)
        hs2dwf.attrs = V1.attrs
        print('hs2dwf',hs2dwf)
        hc2dwf = hc2d*(binv1[1:]+binv1[0:-1])[np.newaxis,:]*0.5
        hc2dwf.attrs = V1.attrs
        hs2dwf[1:,1:].plot(xscale='log',yscale='log',ax=axes[0])
        hc2dwf[1:,1:].plot(xscale='log',yscale='log',ax=axes[1])
        plt.tight_layout()
        plt.savefig(prefix+'2dwt_'+Vname+'_histo'+'.jpg',format='jpg')
        plt.show()
        # next calculations not currently used
        # weighted  histogram
        #whs = hs*dbins[:,np.newaxis]*dbinv1[np.newaxis,:]
        #print('whs',whs)
        # sub histogram
        #hs_sub = hs.sel(PRECS_bin=slice(0.,1.))
        #print('sub hist of PRECCS btw 0., 1. mm/d (limits sampled bin midpoint val)',hs_sub)
```

```python
Varlist = np.array(['num_a3SFSBS','num_a3SFSBC','num_a3SFSIS','num_a3SFSIC',
                    'num_c3SFSBS','num_c3SFSBC','num_c3SFSIS','num_c3SFSIC'])
#Varlist = [Varlist[0]]
for Vname in Varlist:
    V1 = -xr_getvar(Vname,DS,regtag=regtag).where(pmask)
    print(Vname,'Range', V1.attrs['long_name'],V1.min().values, V1.max().values)
    print('sum %5.3f %s ' % (V1.sum().values*1.e-6, 'million'))
    if V1.min().values == V1.max().values:
        print('skipping histogram ')
    else:
        fig, axes = plt.subplots(ncols=2
                        ,gridspec_kw={'width_ratios': [1, 1]}
                        ,figsize=(14,5)
                        )
        v1max = V1.max().values
        v1d4 = v1max*1.e-4
        bins = np.exp(findNiceContours(np.log(np.array([0.01,30.])),nlevs=30))
        binv1 = np.exp(findNiceContours(np.log(np.array([v1d4,v1max])),nlevs=30))
        #print('binv1',binv1)
        hs = histogram(PRECS, V1, bins=[bins, binv1])
        hs.plot(xscale='log',yscale='log',ax=axes[0])
        hc = histogram(PRECC, V1, bins=[bins, binv1])
        hc.plot(xscale='log',yscale='log',ax=axes[1])
        #plt.plot([0.01,30.],[0.01,30.],color='white')
        #plt.plot(np.array([0.01,30.]),np.array([0.01,30./10]),color='red')
        plt.show()
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
