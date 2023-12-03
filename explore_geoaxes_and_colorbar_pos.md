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

import matplotlib as mpl

plotproj = ccrs.Mollweide()
plotproj._threshold /= 100.
ylabels=False
cbar='default'
cbartitle='cbartitle'
units='units'
fig, axes = plt.subplots(ncols=1,
                         gridspec_kw={'width_ratios': [1]},
                         subplot_kw={'projection': plotproj},
                         figsize=(6,6),
                        )
fig.set_dpi(300.0)
axes.patch.set_alpha(0.0)



dataproj=ccrs.PlateCarree()    # data is always assumed to be lat/lon

cmap = plt.get_cmap()

extend = 'both'
clevs = [0.,10.,20.,30.]
units = 'K'
norm = mpl.colors.BoundaryNorm(clevs,cmap.N,extend=extend)
title = 'A title'

axes.set_title(title)
x = np.arange(-180., 190, 10)
y = np.arange(-90., 100, 10)
xv, yv = np.meshgrid(x, y)
data = 30*np.exp(-(xv/40.)**2 - (yv/20.)**2)
gl = axes.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,
                  linewidth=2, color='gray', alpha=0.5)
gl.left_labels=ylabels
gl.right_labels=ylabels
if True:
    pl = axes.contourf(xv, yv, data, levels=clevs, # vmin=zmin, vmax=zmax,
                        norm=norm, cmap=cmap,
                        extend=extend, transform=ccrs.PlateCarree())

x = axes.coastlines(linewidth=1,color='blue')

## Find the location of the main plot axes
## has to be done after some plotting is done in projection space
posn = axes.get_position()
print('posn',posn)


# print some registration marks to help in lining up figures
ax2 = fig.add_axes([0,0,0.1,0.1])
ax2.set_position([posn.x0-0.005, posn.y0-0.005, posn.width+0.01, posn.height+0.01])
ax2.patch.set_alpha(0.0)
ax2.scatter([0,0,1,1], [0,1,0,1], c="r", s=100)
ax2.set_axis_off()
ax2.set_xlim([0,1])
ax2.set_ylim([0,1])

# Add colorbar to plot
# separate colorbar axis following this example
# https://stackoverflow.com/questions/30030328/correct-placement-of-colorbar-relative-to-geo-axes-cartopy

## Where do you want the colorbar? 'bottom' or 'right' as examples
cbar_loc = 'bottom'
#cbar_loc = 'right'
cbar_loc = 'off'
cbar_lab = 'Variable Name [units]'

if cbar_loc != 'off':
    cax = axes
    ## Create colorbar axes (temporarily) anywhere
    cax = fig.add_axes([0,0,0.1,0.1])

## Adjust the positioning and orientation of the colorbar, and draw it
if cbar_loc == 'bottom':
    cax.set_position([posn.x0, posn.y0-0.06, posn.width, 0.03])
    cb = plt.colorbar(
         pl, orientation='horizontal',ticks=clevs,cax=cax,
         label='%s (%s)'%(cbartitle, units)
         )
    cb.ax.tick_params(labelsize=8)
    # works plt.colorbar(pl,cax=cax, orientation='horizontal', label=cbar_lab)
    #plt.colorbar(pl,cax=cax, orientation='horizontal', label=cbar_lab)
elif cbar_loc == 'right':
    cax.set_position([posn.x0+posn.width+0.05, posn.y0, 0.04, posn.height])
    plt.colorbar(pl, cax=cax, orientation='vertical', label=cbar_lab)
    


#y = fig.patch.set_facecolor('xkcd:mint green') # This changes the background
plt.savefig('transparent.pdf',format='pdf',dpi=300,transparent=True)#,facecolor='xkcd:mint green')
#plt.savefig('transparent.png', format='PNG', transparent = True)
z = plt.show()
```

```python
import matplotlib as mpl
#mpl.use('pgf')

import matplotlib.pyplot as plt

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(range(10))
#fig.savefig(f"test.png", transparent=True)
fig.savefig(f"test.pdf", transparent=True)

```
