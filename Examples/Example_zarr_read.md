---
jupyter:
  jupytext:
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

```python
import xarray as xr
import fsspec

def s3load_zarr(s3path):
    return xr.open_dataset(fsspec.get_mapper(s3path),engine = 'zarr')

experiment = 'b.e21.BSSP245smbb_MCBsse2.5Tg_R1.LE2-1031.002'
var = 'TREFHT'
s3_path = f's3://mcb-zarr/CESM2/{experiment}/atm/proc/tseries/month_1/{experiment}.cam.h0.{var}.zarr'

## data is a lazy-loaded xarray dataset
data = s3load_zarr(s3_path)
```
