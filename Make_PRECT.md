---
jupyter:
  jupytext:
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.15.2
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
# mess around with the file names until it matches the correct convention, then set the last logical to
# true to write out the dataset
if True:
    # create a PRECT dataset from PRECL and PRECC
    
    def bld_fname_e1(casename, Varname):
        fname = '/e3sm_prod/phil/timeseries/e3sm-reshaped/'+casename+"/"+casename+".eam.h0.2015-*."+Varname+".nc"
        return fname

    def bld_fname_e2(casename, Varname):
        fname = "/e3sm_prod/phil/timeseries/e3sm-reshaped/"+casename+"/"+Varname+"_201501_*.nc"
        return fname

    def bld_fname_out(casename, Varname):
        #fname = '/e3sm_prod/phil/timeseries/e3sm-reshaped/'+casename+"/"+casename+".eam.h0.2015-2029."+Varname+".nc"
        #fname = '/e3sm_prod/phil/timeseries/e3sm-reshaped/'+casename+"/"+casename+".eam.h0.2015-2046."+Varname+".nc"
        fname = "/e3sm_prod/phil/timeseries/e3sm-reshaped/"+casename+"/"+Varname+"_201501_204412.nc"
        return fname
    
    casename_ctl = '20221014.v2.LR.WCYCLSSP245.E2_CNTL_01'
    casename1 = casename_ctl
    #casename_ptb='20230724.v2.LR.WCYCLSSP245.MCB-SSLT-EM.R1-3.test01'
    casename_ptb='20231122.v2.LR.WCYCLSSP245.MCB-SSLT-EM.R1245.test01'
    casename_ptb='20230714.v2.LR.WCYCLSSP245.MCB-SSLT-EM.R1.test01'
    casename_ptb='20230724.v2.LR.WCYCLSSP245.MCB-SSLT-EM.R5.test01'
    casename1 = casename_ptb

    ind1 = bld_fname_e2(casename1, 'PRECL')
    print('ind1',ind1)
    DS1 = xr.open_mfdataset(ind1)
    ind2 = bld_fname_e2(casename1, 'PRECC')
    DS2 = xr.open_mfdataset(ind2)
    PRECT = DS1['PRECL']+DS2['PRECC']
    PRECT = PRECT.rename('PRECT')
    PRECT.attrs['long_name'] = 'Total Precipitation'

    time_bnds = DS1['time_bnds']
    print('time_bnds last year month',time_bnds.values[-1,0])
    ind3 = bld_fname_out(casename1, 'PRECT')
    #ind3 = ind1.replace('PRECL','PRECT')
    print('ind3 output',ind3)
    # Save one DataArray as dataset
    DSOUT = PRECT.to_dataset(name = 'PRECT')
    # Add second DataArray to existing dataset (ds)
    DSOUT['time_bnds'] = time_bnds
if True:
    DSOUT.to_netcdf(ind3)
```

```python

```
