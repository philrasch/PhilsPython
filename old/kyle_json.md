---
jupyter:
  jupytext:
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.15.2
  kernelspec:
    display_name: Python [conda env:.conda-pjrpy3] *
    language: python
    name: conda-env-.conda-pjrpy3-py
---

```python
import json
import xarray as xr
import glob
d_varlist = json.load(open('/home/jupyter-adminphil/misc/e3sm_variablelist.json'))
#print('d_varlist',d_varlist)
simtype = 'coupled'
path = '/home/ec2-user/e3sm_cases/20221014.v2.LR.WCYCLSSP245.E2_CNTL_01/run/'
path = '/scratch2/geostrat/E3SMv2_AWS/20221014.v2.LR.WCYCLSSP245.E2_CNTL_01/run/'
for component in d_varlist[simtype]:
    for freq in d_varlist[simtype][component]:
        print(component,freq)
        if True: continue
        if freq == 'module': continue;
        print(component,freq)
        fname = path + '20221014.v2.LR.WCYCLSSP245.E2_CNTL_01.{0}.{1}.2016-01*'.format(d_varlist[simtype][component]['module'],freq)
        #print('fname',fname)
        l_fname = glob.glob(fname)
        #print('l_fname',l_fname)
        fname = l_fname[0]
        print('fname',fname)
        data = xr.open_dataset(fname)
        vlist = list(data.variables)
        diff = [var for var in vlist if not var in d_varlist[simtype][component][freq]['varlist']]
        print('diff',diff)                                                                                                                                                                                                                                                                  


```

```python

```
