---
jupyter:
  jupytext:
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.15.2
  kernelspec:
    display_name: Python 3 (ipykernel)
    language: python
    name: python3
---

```python
import numpy as np
import regionmask
import cartopy.crs as ccrs
import matplotlib.pyplot as plt

regionmask.__version__

regList = []
nameList = []
abbrevList = []

name, abbrev, reg = "North_Pacific", "NP", np.array([[160.,30.],[-130.%360,30.],[-130.%360,50.],[160.,50.]])
name, abbrev, reg = "North_Pacific", "NP", np.array([[170.,30.],[-120.%360,30.],[-120.%360,50.],[170.,50.]])
regList, nameList, abbrevList = regList + [reg], nameList + [name], abbrevList + [abbrev]

name, abbrev, reg = "SE_Pacific", "SEP", np.array([[250.,-30.],[290.,-30.],[290.,0.],[250.,0.]])
regList, nameList, abbrevList = regList + [reg], nameList + [name], abbrevList + [abbrev]

name, abbrev, reg = "S_Pacific", "SP", np.array([[-170%360.,-50.],[-90.%360,-50.],[-90.%360,-30.],[-170%360.,-30.]])
regList, nameList, abbrevList = regList + [reg], nameList + [name], abbrevList + [abbrev]

name, abbrev, reg = "NE_Pacific", "NEP", np.array([[-160%360.,0.],[-120.%360,0.],[-120.%360,30.],[-160%360.,30.]])
name, abbrev, reg = "NE_Pacific", "NEP", np.array([[-150%360.,0.],[-110.%360,0.],[-110.%360,30.],[-150%360.,30.]])
regList, nameList, abbrevList = regList + [reg], nameList + [name], abbrevList + [abbrev]

name, abbrev, reg = "SE_Atlantic", "SEA", np.array([[-25.,-30.],[15,-30.],[15,0.],[-25.,0.]])
regList, nameList, abbrevList = regList + [reg], nameList + [name], abbrevList + [abbrev]

o = 360.
name, abbrev, reg = "Northern_Oceans", "NO", np.array([[-179.9+o,60.],[179.9+o,60.],[179.+o,90.],[-179.9+o,90.]])
regList, nameList, abbrevList = regList + [reg], nameList + [name], abbrevList + [abbrev]

SDregions = regionmask.Regions(regList, names=nameList, abbrevs=abbrevList, name="SEED",overlap=True)
SDregions

projection = ccrs.Mollweide(central_longitude=-80)
ax = SDregions.plot(label="abbrev",projection=projection)
ax.set_global()
ax.coastlines(linewidth=0.2,color='blue')
varname = 'regions'
plotfile = '%s.pdf'%varname
plt.savefig(plotfile, bbox_inches='tight')
print('plotfile',plotfile)
```
