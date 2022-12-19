---
jupyter:
  jupytext:
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.14.4
  kernelspec:
    display_name: Python [conda env:pjrpy]
    language: python
    name: conda-env-pjrpy-py
---

```python
import sys
print(sys.version)
%matplotlib inline
%run -i ~/Python/pjr3.py
from jupytext.config import find_jupytext_configuration_file
print('jupytext config file is ',find_jupytext_configuration_file('.'))
```

```python
da = xr.DataArray(np.sin(0.3 * np.arange(12).reshape(4, 3)),
                  [('time', np.arange(4)), 
                   ('space', [0.1, 0.2, 0.3])])
```

```python
n = 1
da = xr.DataArray(np.sin(np.linspace(0, 2 * np.pi, n)), dims='x',
                  coords={'x': np.linspace(0, 1, n)})
print('da',da)
da.plot.line('o', label='original')
```

```python
da.plot()
```

```python
da.plot.hist()
```

```python
import xarray as xr
import matplotlib.pyplot as plt

plt.close('all')
data = xr.DataArray(np.random.randn(2, 3), dims=('x', 'y'),coords={'x': [10, 20],'y' : [0,50,100]})
data.plot() # Will produce a square plot

```

```python
data.plot.contourf()
```
