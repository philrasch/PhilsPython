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
import sys
print(sys.version)
%matplotlib inline
%run -i ~/Python/pjr3
from jupytext.config import find_jupytext_configuration_file
print('jupytext config file is ',find_jupytext_configuration_file('.'))
```

```python
def setfig3b1x1 ():
    """
    return fig and axes for a single panel figure
    """
    fig, axes = plt.subplots(ncols=1,
                             gridspec_kw={'width_ratios': [1]},
                             subplot_kw={'projection': ccrs.Mollweide()},
                             figsize=(6,3),
                            )
    fig.set_dpi(300.0)
    return fig, axes;
```

```python
from matplotlib.backends.backend_pdf import PdfPages
with PdfPages('test.pdf') as pdf:
    for page in np.arange(0,3):
        print('page ',page)

        fig, axes = setfig3b1x1()
        ax = axes
        ax.coastlines(linewidth=1,color='blue')
        title =' Page (%5.0f)' % page
        ax.set_title(title)
        #xr_llhplot(V1, ax=axes,clevs=clevs,title=pref1+sV1A, cbartitle=cbartitle)
        pdf.savefig(dpi=300)
        plt.show()
plt.close()
```

```python
infile = "test.pdf"
add_prov(infile)
```

```python

```
