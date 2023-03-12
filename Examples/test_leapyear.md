---
jupyter:
  jupytext:
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
import xarray as xr
import numpy as np
 
indir = '/global/cfs/cdirs/m1199/wilbert/JRA55-do/ts/ts_input4MIPs_atmosphericState_OMIP_MRI-JRA55-do-1-5-0_gr_201601010000-201612312100.nc'
print(indir)
#print('exists',os.path.exists(indir))
DS = xr.open_mfdataset(indir).chunk({'time': 20})
#print(DS)

Var = DS.ts
#print(Var.shape)
print('Var',Var)
# select times that are Feb 29
f29 = ((Var.time.dt.month == 2) & (Var.time.dt.day == 29))
if len(np.where(f29)[0]) > 0:
    myinds = np.where(~f29)[0]
    print('mydays',myinds,np.shape(myinds))
    Var_nof29 = Var.isel(time=myinds)
    print('Var_nof29',Var_nof29)
    Var_nof29.to_netcdf('/tmp/test.nc')
else:
    print('nofeb29')

```
