# module of phils useful stuff
"""Phils useful functions"""

import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as mcolors
import numpy as np
import string
import copy
import cartopy.crs as ccrs
import xarray as xr
import os
import sys
import glob
#import cdms2
#import cdutil
from scipy.interpolate import interp1d
from scipy.interpolate import griddata
from mpl_toolkits.axes_grid1 import make_axes_locatable
#from metrics.packages.amwg.derivations import *
#from metrics.packages.amwg.derivations import qflx_lhflx_conversions as flxconv

# from mpl_toolkits.basemap import Basemap
#from netCDF4 import Dataset, date2index
from datetime import datetime
import matplotlib.colors as mcolors


import warnings

from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()

# original code from 
# https://stackoverflow.com/questions/28934767/best-way-to-interpolate-a-numpy-ndarray-along-an-axis
def interp_along_axis(y, x, newx, axis, inverse=False, method='linear', verbose=False):
    """ Interpolate vertical profiles, e.g. of atmospheric variables
    using vectorized numpy operations

    This function assumes that the x-xoordinate increases monotonically

    ps:
    * Updated to work with irregularly spaced x-coordinate.
    * Updated to work with irregularly spaced newx-coordinate
    * Updated to easily inverse the direction of the x-coordinate
    * Updated to fill with nans outside extrapolation range
    * Updated to include a linear interpolation method as well
        (it was initially written for a cubic function)

    Peter Kalverla
    March 2018

    --------------------
    More info:
    Algorithm from: http://www.paulinternet.nl/?page=bicubic
    It approximates y = f(x) = ax^3 + bx^2 + cx + d
    where y may be an ndarray input vector
    Returns f(newx)

    The algorithm uses the derivative f'(x) = 3ax^2 + 2bx + c
    and uses the fact that:
    f(0) = d
    f(1) = a + b + c + d
    f'(0) = c
    f'(1) = 3a + 2b + c

    Rewriting this yields expressions for a, b, c, d:
    a = 2f(0) - 2f(1) + f'(0) + f'(1)
    b = -3f(0) + 3f(1) - 2f'(0) - f'(1)
    c = f'(0)
    d = f(0)

    These can be evaluated at two neighbouring points in x and
    as such constitute the piecewise cubic interpolator.
    """
    # View of x and y with axis as first dimension
    if inverse:
        _x = np.moveaxis(x, axis, 0)[::-1, ...]
        _y = np.moveaxis(y, axis, 0)[::-1, ...]
        _newx = np.moveaxis(newx, axis, 0)[::-1, ...]
    else:
        _y = np.moveaxis(y, axis, 0)
        _x = np.moveaxis(x, axis, 0)
        _newx = np.moveaxis(newx, axis, 0)
    # Sanity checks
    if np.any(_newx[0] < _x[0]) or np.any(_newx[-1] > _x[-1]):
        # raise ValueError('This function cannot extrapolate')
        if verbose:
            warnings.warn("Some values are outside the interpolation range. "
                          "These will be filled with NaN")
    if np.any(np.diff(_x, axis=0) < 0):
        raise ValueError('x should increase monotonically')
    if np.any(np.diff(_newx, axis=0) < 0):
        raise ValueError('newx should increase monotonically')
    # Cubic interpolation needs the gradient of y in addition to its values
    if method == 'cubic':
        # For now, simply use a numpy function to get the derivatives
        # This produces the largest memory overhead of the function and
        # could alternatively be done in passing.
        ydx = np.gradient(_y, axis=0, edge_order=2)

    # This will later be concatenated with a dynamic '0th' index
    ind = [i for i in np.indices(_y.shape[1:])]
    # Allocate the output array
    original_dims = _y.shape
    newdims = list(original_dims)
    newdims[0] = len(_newx)
    newy = np.zeros(newdims)

    # set initial bounds
    i_lower = np.zeros(_x.shape[1:], dtype=int)
    i_upper = np.ones(_x.shape[1:], dtype=int)
    x_lower = _x[0, ...]
    x_upper = _x[1, ...]
    for i, xi in enumerate(_newx):
        # Start at the 'bottom' of the array and work upwards
        # This only works if x and newx increase monotonically

        # Update bounds where necessary and possible
        needs_update = (xi > x_upper) & (i_upper+1<len(_x))
        # print x_upper.max(), np.any(needs_update)
        while np.any(needs_update):
            i_lower = np.where(needs_update, i_lower+1, i_lower)
            i_upper = i_lower + 1
#orig            x_lower = _x[[i_lower]+ind]
#orig            x_upper = _x[[i_upper]+ind]
            x_lower = _x[tuple([i_lower]+ind)]
            x_upper = _x[tuple([i_upper]+ind)]
            # Check again
            needs_update = (xi > x_upper) & (i_upper+1<len(_x))
        # Express the position of xi relative to its neighbours
        xj = (xi-x_lower)/(x_upper - x_lower)

        # Determine where there is a valid interpolation range
        within_bounds = (_x[0, ...] < xi) & (xi < _x[-1, ...])
        if method == 'linear':
#orig            f0, f1 = _y[[i_lower]+ind], _y[[i_upper]+ind]
            f0, f1 = _y[tuple([i_lower]+ind)], _y[tuple([i_upper]+ind)]
            a = f1 - f0
            b = f0

            newy[i, ...] = np.where(within_bounds, a*xj+b, np.nan)

        elif method=='cubic':
#orig            f0, f1 = _y[[i_lower]+ind], _y[[i_upper]+ind]
#orig            df0, df1 = ydx[[i_lower]+ind], ydx[[i_upper]+ind]
            f0, f1 = _y[tuple([i_lower]+ind)], _y[tuple([i_upper]+ind)]
            df0, df1 = ydx[tuple([i_lower]+ind)], ydx[tuple([i_upper]+ind)]
            a = 2*f0 - 2*f1 + df0 + df1
            b = -3*f0 + 3*f1 - 2*df0 - df1
            c = df0
            d = f0

            newy[i, ...] = np.where(within_bounds, a*xj**3 + b*xj**2 + c*xj + d, np.nan)

        else:
            raise ValueError("invalid interpolation method"
                             "(choose 'linear' or 'cubic')")
    if inverse:
        newy = newy[::-1, ...]

    return np.moveaxis(newy, 0, axis)

import numpy as np
import warnings

def interp_along_axis_V2(y, x, newx, axis, inverse=False, method='linear'):
    """ Interpolate vertical profiles, e.g. of atmospheric variables
    using vectorized numpy operations

    This function assumes that the x-xoordinate increases monotonically

    ps:
    * Updated to work with irregularly spaced x-coordinate.
    * Updated to work with irregularly spaced newx-coordinate
    * Updated to easily inverse the direction of the x-coordinate
    * Updated to fill with nans outside extrapolation range
    * Updated to include a linear interpolation method as well
        (it was initially written for a cubic function)

    Peter Kalverla
    March 2018

This version was downloaded from
https://stackoverflow.com/questions/28934767/best-way-to-interpolate-a-numpy-ndarray-along-an-axis


    --------------------
    More info:
    Algorithm from: http://www.paulinternet.nl/?page=bicubic
    It approximates y = f(x) = ax^3 + bx^2 + cx + d
    where y may be an ndarray input vector
    Returns f(newx)

    The algorithm uses the derivative f'(x) = 3ax^2 + 2bx + c
    and uses the fact that:
    f(0) = d
    f(1) = a + b + c + d
    f'(0) = c
    f'(1) = 3a + 2b + c

    Rewriting this yields expressions for a, b, c, d:
    a = 2f(0) - 2f(1) + f'(0) + f'(1)
    b = -3f(0) + 3f(1) - 2f'(0) - f'(1)
    c = f'(0)
    d = f(0)

    These can be evaluated at two neighbouring points in x and
    as such constitute the piecewise cubic interpolator.
    """

    # View of x and y with axis as first dimension
    if inverse:
        _x = np.moveaxis(x, axis, 0)[::-1, ...]
        _y = np.moveaxis(y, axis, 0)[::-1, ...]
        _newx = np.moveaxis(newx, axis, 0)[::-1, ...]
    else:
        _y = np.moveaxis(y, axis, 0)
        _x = np.moveaxis(x, axis, 0)
        _newx = np.moveaxis(newx, axis, 0)

    # Sanity checks
    if np.any(_newx[0] < _x[0]) or np.any(_newx[-1] > _x[-1]):
        # raise ValueError('This function cannot extrapolate')
        warnings.warn("Some values are outside the interpolation range. "
                      "These will be filled with NaN")
    if np.any(np.diff(_x, axis=0) < 0):
        raise ValueError('x should increase monotonically')
    if np.any(np.diff(_newx, axis=0) < 0):
        raise ValueError('newx should increase monotonically')

    # Cubic interpolation needs the gradient of y in addition to its values
    if method == 'cubic':
        # For now, simply use a numpy function to get the derivatives
        # This produces the largest memory overhead of the function and
        # could alternatively be done in passing.
        ydx = np.gradient(_y, axis=0, edge_order=2)

    # This will later be concatenated with a dynamic '0th' index
    ind = [i for i in np.indices(_y.shape[1:])]

    # Allocate the output array
    original_dims = _y.shape
    newdims = list(original_dims)
    newdims[0] = len(_newx)
    newy = np.zeros(newdims)

    # set initial bounds
    i_lower = np.zeros(_x.shape[1:], dtype=int)
    i_upper = np.ones(_x.shape[1:], dtype=int)
    x_lower = _x[0, ...]
    x_upper = _x[1, ...]

    for i, xi in enumerate(_newx):
        # Start at the 'bottom' of the array and work upwards
        # This only works if x and newx increase monotonically

        # Update bounds where necessary and possible
        needs_update = (xi > x_upper) & (i_upper+1<len(_x))
        # print x_upper.max(), np.any(needs_update)
        while np.any(needs_update):
            i_lower = np.where(needs_update, i_lower+1, i_lower)
            i_upper = i_lower + 1
            x_lower = _x[[i_lower]+ind]
            x_upper = _x[[i_upper]+ind]

            # Check again
            needs_update = (xi > x_upper) & (i_upper+1<len(_x))

        # Express the position of xi relative to its neighbours
        xj = (xi-x_lower)/(x_upper - x_lower)

        # Determine where there is a valid interpolation range
        within_bounds = (_x[0, ...] < xi) & (xi < _x[-1, ...])

        if method == 'linear':
            f0, f1 = _y[[i_lower]+ind], _y[[i_upper]+ind]
            a = f1 - f0
            b = f0

            newy[i, ...] = np.where(within_bounds, a*xj+b, np.nan)

        elif method=='cubic':
            f0, f1 = _y[[i_lower]+ind], _y[[i_upper]+ind]
            df0, df1 = ydx[[i_lower]+ind], ydx[[i_upper]+ind]

            a = 2*f0 - 2*f1 + df0 + df1
            b = -3*f0 + 3*f1 - 2*df0 - df1
            c = df0
            d = f0

            newy[i, ...] = np.where(within_bounds, a*xj**3 + b*xj**2 + c*xj + d, np.nan)

        else:
            raise ValueError("invalid interpolation method"
                             "(choose 'linear' or 'cubic')")

    if inverse:
        newy = newy[::-1, ...]

    return np.moveaxis(newy, 0, axis)


def make_colormap(seq):
    """Return a LinearSegmentedColormap
    seq: a sequence of floats and RGB-tuples. The floats should be increasing
    and in the interval (0,1).
    """
    seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
#   print "seq",seq
    cdict = {'red': [], 'green': [], 'blue': []}
    for i, item in enumerate(seq):
        if isinstance(item, float):
            r1, g1, b1 = seq[i - 1]
            r2, g2, b2 = seq[i + 1]
            cdict['red'].append([item, r1, r2])
            cdict['green'].append([item, g1, g2])
            cdict['blue'].append([item, b1, b2])
    return mcolors.LinearSegmentedColormap('CustomMap', cdict)


def diverge_map(high=(1., 0., 0.), low=(0.0, 0., 1.)):
#def diverge_map(high=(0.565, 0.392, 0.173), low=(0.094, 0.310, 0.635)):
    '''
    low and high are colors that will be used for the two
    ends of the spectrum. they can be either color strings
    or rgb color tuples
    '''
    c = mcolors.ColorConverter().to_rgb
    if isinstance(low, basestring): low = c(low)
    if isinstance(high, basestring): high = c(high)
    return make_colormap([low, c('white'), 0.5, c('white'), high])

# frequently used constants
gravity = 9.80665 # std gravity
pi = np.pi  # pi to 
arad = 6.37e6 # radius of earth in meters

        
def onedplot(lat1, data1, lat2, data2):
    """Make 2 lines and the difference plot between the lines."""
    plt.rc('lines', linewidth=4)
    fig, (ax0, ax1)  = plt.subplots(nrows=2)
    # You can get some information about the variables by
    #lat2.info()
    #data2.info()
    # You can also find out what class the object x or y belong to
    print ("class", data2.__class__)
    print ("lat2.size", np.size(lat2))
    print ("aa")
    print ("lat1.size", np.size(lat1))
    print ("bb")
    plt.rc('axes', color_cycle=['r', 'g', 'b', 'y'])
    x2 = np.array(lat2);
    y2 = np.array(data2);
    x1 = np.array(lat1)
    y1 = np.array(data1)
    l1, = ax0.plot(x1, y1,label='CORE.v2')
    l2, = ax0.plot(x2, y2,label='EBAFS-Sfc')
    print ("cc")
    ax0.legend(handles=[l1,l2])
    #ax0.plot(lat2,data2)
    ax0.set_title('Downward Shortwave')
    print ("xxx", np.size(lat2),  np.size(lat1))
    if np.size(lat2) > np.size(lat1):
        xh = x2;
        y1 = np.interp(xh, x1, y1)
    else:
        xh = x1;
        y2 = np.interp(xh, x2, y2)
        print ("regrid data2.size", np.size(y2))
        print ("regrid data1.size", np.size(y1))
    dy = np.array(y2-y1);
    ax1.plot(xh, dy);
        
def findNiceContours(data,nlevs=None,rmClev=None,sym=None,verbose=None):
    """Find Nice Contours
    data = 2d numpy array (or data structure base on numpy) ordered (latitude, pressure)
    nlevs = approximate number of contour levels to return (default 10)
    rmClev = if defined delete the contour level near this value
    sym = if defined make the contour intervals symmetric about zero
    verbose = if defined, print out some info to help debug
    """
    if nlevs is None: nlevs = 9

    dsub = data[np.isfinite(data)] # ignore NaN
    zmax = dsub.max()
    zmin = dsub.min()
    clevs = np.linspace(zmin, zmax, nlevs)
#    if nozero is None: nozero=0
#    if sym is None: sym=0

    if not verbose is None: print ("zmin, zmax", zmin, zmax)
    if zmax == zmin: 
        if (zmax != 0): 
          return np.array([0.9*zmax,zmax,1.1*zmax]);
        else:
          return np.array([-0.1,0.0, 0.1]);
    if not sym is None: 
        zcom = max(abs(zmin),abs(zmax));
        zmin = -zcom;
        zmax = +zcom;
        if not verbose is None: print ("zmin, zmax in sym", zmin, zmax)
    zrange = zmax/(zmin+1.e-20)
    okints = [0.1, 0.2, 0.5,1,2,5,10,20] # acceptable values for contour intervals
    okints = np.array(okints)
    zinc = (zmax-zmin)/nlevs;
    if not verbose is None: print ("zinc", zinc)
    pow = int(np.log10(zinc))
    if not verbose is None: print ("pow", pow)
    cints = okints*(10.**pow);
    if not verbose is None: print ("cints", cints)
    nlevsout = (zmax-zmin)/cints;
    if not verbose is None: print ("nlevsout",nlevsout)
#    np.where(V0 == V0.max())
    f1 = abs(nlevsout-nlevs)
    if not verbose is None: print ("f1", f1)
    f2 = np.where(f1 == f1.min());
    if not verbose is None: print ("f2", f2)
    nlevbest = int(f2[0])
    if not verbose is None: print ("nlevbest, cintsbest",nlevbest, cints[nlevbest])
    if not verbose is None: print ("f3", zmin/cints[nlevbest])
    zminout = int(zmin/cints[nlevbest])*cints[nlevbest];
    zmaxout = int(zmax/cints[nlevbest])*(cints[nlevbest]);
    ninc = int((zmaxout-zminout)/cints[nlevbest] + 1.00000001);
    if not verbose is None: print ("ninc, zminout, zmaxout", ninc, zminout, zmaxout)
    clevs = np.arange(zminout, zmaxout*1.001, cints[nlevbest]);
    if not rmClev is None:    
        f4 = abs((clevs-rmClev)/max(abs(clevs)))
        if not verbose is None: print ("f4", f4)
        alist = np.argmin(f4);
        if not verbose is None: print ("alist", alist)
        if (np.size(alist) > 0): 
            zlist = np.where(clevs != clevs[alist])
            if not verbose is None: print ("zlist", zlist)
            if (np.size(zlist)> 0): 
                if not verbose is None: print ("list of nonzero contour levels", zlist)
                clevs = clevs[zlist];
    return clevs;


def plotZMf(data, x, y, plotOpt=None, modelLevels=None, surfacePressure=None, axesa=None, fig=None):
    """Create a zonal mean contour plot of one variable
    axesa = the axes that we make the plot on 
    data = 2d numpy array (or data structure base on numpy) ordered (latitude, pressure)
    x = 1d numpy array of latitude
    y = 1d numpy array of pressures (or pseudo pressure (like eta))
    plotOpt is a optional dictionary with plotting options:
      'scale_factor': multiply values with this factor before plotting
      'units': a units label for the colorbar
      'clevs': use list of values as contour intervals
      'cmap': the color map to use
      'cabv': the above color
      'cbel': the below color
      'colorbar': location of colorbar ('bot','top','left','right','None')
      'rmClev': contour level to delete; frequently Zero, see findNiceContours
      'title': a title for the plot
      'ybot': if present, the pressure at the plot bottom
      'ytop': if present, the pressure at the top
    modelLevels:  If present a small side panel will be drawn with lines for each model level
    surfacePressure: a list (dimension len(x)) of surface pressure values. If present, these will
        be used to mask out regions below the surface
    """
    # explanation of axes:
    #   ax1: primary coordinate system latitude vs. pressure (left ticks on y axis)
    #   ax2: twinned axes for altitude coordinates on right y axis
    #   axm: small side panel with shared y axis from ax2 for display of model levels
    # right y ticks and y label will be drawn on axr if modelLevels are given, else on ax2
    #   axr: pointer to "right axis", either ax2 or axm

    if plotOpt is None: plotOpt = {}

    labelFontSize = "small"
    # create figure and axes if they are not supplied
    if fig is None:
        fig = plt.gcf() # get current figure
    if axesa is None: axesa = plt.gca()
    ax1 = axesa


    # scale data if requested
    scale_factor = plotOpt.get('scale_factor', 1.0)
    pdata = data * scale_factor
    
    # determine contour levels to be used; nice contours are chosen if not supplied
    # the contour value rmClev will be removed if requested (often the zero contour)
    clevs = plotOpt.get("clevs") 
    rmClev = plotOpt.get("rmClev")

#    print "plotzmf rmClev", rmClev
    if clevs is None:
        clevs = findNiceContours(data,rmClev=rmClev)
#    print "data range",data.min(), data.max()
#    print "clevs", clevs

# choose a colormap if not supplied
    cmap = plotOpt.get("cmap")
    if cmap is None:
        cmap = mpl.cm.get_cmap()
    cmap = copy.copy(cmap)
    norm = mpl.colors.BoundaryNorm(clevs,cmap.N)

    # draw the (filled) contours
    contour = ax1.contourf(x, y, pdata, levels=clevs, norm=norm,
                        cmap=cmap,
                        extend='both')
    cabv = plotOpt.get('cabv','yellow')
    cbel = plotOpt.get('cbel','magenta')
    contour.cmap.set_over(cabv)
    contour.cmap.set_under(cbel)
    # mask out surface pressure if given
    if not surfacePressure is None: 
        ax1.fill_between(x, surfacePressure, surfacePressure.max(), color="white")    
    # add a title
    if hasattr(data,'long_name'):
        deftitle = data.long_name
#        print "deftitle set to longname", deftitle
        #    variance.units = '(%s)^2'%var.units
    else:
#        print "no long_name attribute"
        deftitle = ''
    ltitle = plotOpt.get('ltitle',deftitle)
    ax1.set_title(ltitle,loc='left')
    defrtitle = ''
    rtitle = plotOpt.get('rtitle',defrtitle)
    ax1.set_title(rtitle,loc='right')

# steal some space and allocate axes for a z coordinate at right of plot
    divider = make_axes_locatable(ax1) # create a divider
    ax_z = divider.new_horizontal(size="1%", pad=0.25,frameon=False) # steal z space from end
    ax_z.xaxis.set_major_locator(plt.NullLocator())
#    ax_z.set_axis_off()
    fig.add_axes(ax_z)

    if hasattr(data,'units'):
        defunits = '('+data.units+')'
 #       print "has units", defunits
        #    variance.units = '(%s)^2'%var.units
    else:
#        print "no units attribute"
        defunits = ''
        # if type(data) is np.ndarray:
        #    print "ndarray"
        #else:
        #    defunits = 'a title'
    units = plotOpt.get('units',defunits)
    # add colorbar
    # Note: use of the ticks keyword forces colorbar to draw all labels
    colorbar = plotOpt.get('colorbar', 'bot')
    fmt = mpl.ticker.FormatStrFormatter("%g")
    if colorbar == 'bot':
        cbar = fig.colorbar(contour, ax=ax1, orientation='horizontal', shrink=1.05, pad=0.2,
                            ticks=clevs, format=fmt)
        cbar.set_label(units)
        for t in cbar.ax.get_xticklabels():
            t.set_fontsize(labelFontSize)
    if (colorbar == 'right' or colorbar == 'rightnd'):
        ax_cb2 = divider.new_horizontal(size="10%", pad=0.3) # steal colorbar space from end
        fig.add_axes(ax_cb2)
        if colorbar == 'right':
            cbar = fig.colorbar(contour, cax=ax_cb2, ticks=clevs, orientation="vertical")
            for t in cbar.ax.get_yticklabels():
                t.set_fontsize(labelFontSize)
            #cbar.set_label(units)
            cbar.ax.set_title(units,pad=5)
        else:
            plt.axis('off')


        # create an axes on the right side of ax. The width of cax will be 5%
        # of ax and the padding between cax and ax will be fixed at 0.05 inch.
        #divider = make_axes_locatable(ax1)
        #cax = divider.append_axes("right", size="5%", pad=1.05)
        #cax = fig.add_axes([0.8, 0.1, 0.03, 0.8]) 
        #cbar = fig.colorbar(contour, ax=ax1, pad=0.2, use_gridspec=True)

    # set up y axes: log pressure labels on the left y axis, altitude labels
    # according to model levels on the right y axis
    ax1.set_ylabel("Pressure [hPa]")
    ax1.set_yscale('log')
    ax1.yaxis.set_label_coords(-0.15, 0.5)
#    print "y", y
    xmesh,ymesh = np.meshgrid(x, y)
#    print "ymesh range", ymesh.min(), ymesh.max()
    ybotd = 10.*np.ceil(ymesh.max()/10.)
    ytopd = ymesh.min()
#    print "ybot, ytop defaults", ybotd, ytopd
    ybot = plotOpt.get('ybot', ybotd)
    ytop = plotOpt.get('ytop', ytopd)
#    ybot = 1000.
#    ytop = 1.
    ax1.set_ylim(ybot, ytop) # avoid truncation of 1000 hPa
    subs = [1,2,5]
    if ybot/ytop < 20.:
        subs = [1,2,3,4,5,6,7,8,9]
#    print "subs", subs
    ax1.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())
    y1loc = mpl.ticker.LogLocator(base=10., subs=subs)
    ax1.yaxis.set_major_locator(y1loc)
    fmt = mpl.ticker.FormatStrFormatter("%g")
    ax1.yaxis.set_major_formatter(fmt)
    for t in ax1.get_yticklabels():
        t.set_fontsize(labelFontSize)
   # change values and font size of x labels
    ax1.set_xlabel('Latitude [degrees]')
    ax1.xaxis.labelpad = 0.05
    xloc = mpl.ticker.FixedLocator(np.arange(-90.,91.,30.))
    ax1.xaxis.set_major_locator(xloc)
    for t in ax1.get_xticklabels():
        t.set_fontsize(labelFontSize)
  
    
    # calculate altitudes from pressure values (use fixed scale height)
    z0 = 8.400    # scale height for pressure_to_altitude conversion [km]
    altitude = z0 * np.log(1000./ymesh)
#    altitude = z0 * np.log(1015.23/ymesh)
    # add second y axis for altitude scale 
#   ax2 = ax1.twinx()
   # draw horizontal lines to the right to indicate model levels
    if not modelLevels is None:
        pos = ax1.get_position()
        ax_z.set_xlim(0., 1.)
        ax_z.xaxis.set_visible(False)
        modelLev = ax_z.hlines(altitude, 0., 1., color='0.5')
        axr = ax_z     # specify y axis for right tick marks and labels
        # turn off tick labels of ax2
        for t in ax_z.get_yticklabels():
            t.set_visible(False)
        label_xcoor = 4.7
    else:
        axr = ax_z
        label_xcoor = 7.
    axr.set_ylabel("Altitude [km]")
    axr.yaxis.set_label_coords(label_xcoor, 0.5)
#    axr.yaxis.labelpad = 42.
    axr.yaxis.set_ticks_position('left')
#   axr.yaxis.set_label_coords(10.,0.5)
    altsub = altitude[(ymesh >= ytop) & (ymesh <= ybot)]
#    altsub = altitude
#    axr.set_ylim(altitude.min(), altitude.max())
#    print "pressure and altitude range", ybot, ytop, altsub.min(), altsub.max()
    axr.set_ylim(altsub.min(), altsub.max())
    yrloc = mpl.ticker.MaxNLocator(steps=[1,2,5,10])
    axr.yaxis.set_major_locator(yrloc)
    axr.yaxis.set_tick_params(pad=0)
    axr.yaxis.tick_left()
    for t in axr.yaxis.get_majorticklines():
        t.set_visible(False)
    for t in axr.get_yticklabels():
        t.set_fontsize(labelFontSize)

def hy2plev(T, P, pout,verbose=False):
    """
    hy2plev(T, P, pout)
    interpolate a field from hybrid to pressure coordinates
        assumes the input fields are xarray dataarrays
        tested for EAM/CAM cubed sphere and lat/lon grids
        
           T the array on hybrid surfaces
           P the pressures on hybrid surfaces 
           pout the pressures to interpolate to (array or scalar)
    """
    #print("P",P)
    pold = P.values
    #print('pold',pold.shape)
    #print('aaa',np.isscalar(pout), np.ndim(pout) == 0 )
    ptmp = np.asarray(pout)
    if np.ndim(ptmp) == 0:
        ptmp = np.array([ptmp])
    nznew = len(ptmp)
   #print('ynew',ynew.shape)
    clist = list(T.coords.keys())
    #print('clist',clist,len(clist))
    Tdims = T.dims
    #print ('Tdims',Tdims)
    Tshape = T.values.shape
#    print('Tshape',Tshape)
    dlist = list(T.dims)
#    print('dlist',dlist,len(dlist))
    ndlist = dlist
    d2list = list()
    olist = list()
    d2size = 1
    axs = 0
    for index, item in enumerate(dlist):
        if item == 'lev':
            ndlist[index] = 'plev'
            d2list.append(nznew)
            axs = index
            nzold = Tshape[index]
            olist.append(nznew)
            #d2size = d2size*nznew
        else:
            d2list.append(Tshape[index])
            d2size = d2size*Tshape[index]
            olist.append(1)
        #print('xxx',index,ndlist[index],item,d2size)

    #print('ndlist',ndlist)
    #print('d2list',d2list,d2size)
    #print('olist',olist)
    #print('axs',axs)
    #print('nzold',nzold)

    #nz,ny,nx = P.shape
    #print('nz,ny,nx',nz,ny,nx)

    #pnew = np.repeat(ptmp,nx*ny).reshape((nznew,ny,nx))
    #pnew = np.repeat(ptmp,d2size).reshape(d2list)
    pnew = np.zeros(d2list) + ptmp.reshape(olist)
    #print('pnewxxx',pnew.shape,pnew[0,:,0], pnew[1,:,1])
    #Tvals = T.values
    #for index in range(nzold):
    #    print('pold(0),Tvals(0), pold(1), Tvales(1) ',index, pold[0,index,ind], Tvals[0,index,ind], pold[1,index,ind], Tvals[1,index,ind])    
    #print('pnew',pnew)
    #print('pnew2',pnew2.shape)
    ynew = interp_along_axis(T.values, pold, pnew, axis=axs, inverse=False, method='linear',verbose=verbose)
    #print('ynew',ynew)
    #for index in range(nznew):
    #    print('pnew(0),ynew(0), pnew(1), ynew(1) ',index, pnew[0,index,ind], ynew[0,index,ind], pnew[1,index,ind], ynew[1,index,ind])    
  
    cnew = dict(T.coords)
    del cnew['lev'] # delete lev
    #print('xxx',type(cnew))
    #print('T coords', cnew)
    #print('yyy')
    #print('plev',cnew['plev'])
    #print('Tlll', T.reset_coords('lev', drop=True))
    dnew = T.dims
    #print('T dims', dnew)
    #print('zzz')
    #print('T attrs', T.attrs)
    NV = xr.DataArray(
        ynew,
        dims=ndlist,
    )
    #print('NV0',NV)
    NV = NV.assign_coords(cnew)
    #print('NV1',NV)
    NV = NV.assign_coords({'plev': ptmp}) # add pcol as a new coordinate
    #print('NV2',NV)
    NV.attrs = T.attrs
    #print('NV3',NV)
    return NV

def interp_to_latlon(data2d,lat,lon,lat_i,lon_i):
    """
    # interpolating in lat/lon space has issues. interpolate in
    # stereographic projection:
    #
    # input:
    #    unstructured 1D data:  data(ncol),lat(ncol),lon(ncol)
    #    target lat/lon grid:   lon_i(nlon), lat_i(nlat)
    #
    # output 2D interpolated data::
    #   data(nlon,nlat)
    #
    """
    # mesh grid
    dproj=ccrs.PlateCarree()
    nhalf = int(len(lat_i)/2)
    lat_south = lat_i[ :nhalf]
    lat_north = lat_i[ nhalf:]

    # take source data in the correct hemisphere, include extra halo points for interpolation
    # using the full global data sometimes confuses griddata with points being mapped close to infinity
    halo = 15 # degrees
    data2d_h=data2d[lat<halo]

    lon_h=lon[lat<halo]
    lat_h=lat[lat<halo]
    xv,yv=np.meshgrid(lon_i,lat_south)
    coords_in  = ccrs.SouthPolarStereo().transform_points(dproj,lon_h,lat_h)
    coords_out = ccrs.SouthPolarStereo().transform_points(dproj,xv.flatten(),yv.flatten())
    data_s = griddata(coords_in[:,0:2], data2d_h, coords_out[:,0:2], method='linear')

    data2d_h=data2d[lat>-halo]
    lon_h=lon[lat>-halo]
    lat_h=lat[lat>-halo]
    xv,yv=np.meshgrid(lon_i,lat_north)
    coords_in  = ccrs.NorthPolarStereo().transform_points(dproj,lon_h,lat_h)
    coords_out = ccrs.NorthPolarStereo().transform_points(dproj,xv.flatten(),yv.flatten())
    data_n = griddata(coords_in[:,0:2], data2d_h, coords_out[:,0:2], method='linear')
    
    data_i=np.concatenate((data_s,data_n)).reshape(len(lat_i),len(lon_i))
    return data_i


def interp_ap(xt, yt, data2d,lat,lon,method=None):
    """
    # interp an arbitrary set of points at xt, yt 
    # from data on an unstructured mesh at lat, lon
    #
    # interpolating in lat/lon space has issues with triangulation 
    # at pole and wrapping at greenwich, so interpolate in stereographic projection:
    #
    # input:
    #    data2d(ncol,...),lat(ncol),lon(ncol): data and coords on unstructured mesh
    #    data2d can be multidimensional array, but ncol must be first coordinate
    #    xt, yt: lat and lon coordinates of locations to interpolate to
    #    method: optional, use cubic interpolation if method='cubic'
    #
    # output 
    #    returns an array with same shape as xt with interpolated data
    #
    """
    from scipy.interpolate import LinearNDInterpolator
    from scipy.interpolate import CloughTocher2DInterpolator
    
    intp2D = LinearNDInterpolator
    if method == 'cubic':
        intp2D = CloughTocher2DInterpolator

    ld = data2d.shape[0] # length of first coord of input data
    lx = lon.shape[0]
    ly = lat.shape[0]
    if ((ld != lx) | (ld != ly)):
        print('inconsistent data2d, lon, lat arrays', ld, lx, ly)
        raise TypeError("inconsistent input in interp_ap")
    
    # mesh grid
    dproj=ccrs.PlateCarree()

    # select interpolation points located in the nh and sh
    inds = np.where(yt <= 0)
    indn = np.where(yt > 0)

    xtn = xt[indn]
    ytn = yt[indn]
    xts = xt[inds]
    yts = yt[inds]

    # take source data in the correct hemisphere, include extra halo points for interpolation
    # using the full global data sometimes confuses interpolation with points being mapped close to infinity
    halo = 15 # degrees
    data2d_h=data2d[lat<halo]

    lon_h=lon[lat<halo]
    lat_h=lat[lat<halo]
    coords_in  = ccrs.SouthPolarStereo().transform_points(dproj,lon_h,lat_h)

    #data_i = np.empty_like(xt,dtype=data2d.dtype)
    #data_i[:] = np.nan
    dims = list(xt.shape)+list(data2d[0,...].shape)
    #print('dims',dims)
    #data_i = np.empty_like(xt,dtype=data2d.dtype)
    data_i = np.zeros(dims,dtype=data2d.dtype)
    data_i[:] = np.nan
    #print ('bbb',data_i.shape)


    data_s = []
    if len(yts) > 0:
        cto = ccrs.SouthPolarStereo().transform_points(dproj,xts,yts)
        interp = intp2D(coords_in[:,0:2], data2d_h)
        data_s = interp(cto[:,0],cto[:,1])
        data_i[inds] = data_s

    data2d_h=data2d[lat>-halo]
    lon_h=lon[lat>-halo]
    lat_h=lat[lat>-halo]
    coords_in  = ccrs.NorthPolarStereo().transform_points(dproj,lon_h,lat_h)

    data_n = []
    if len(ytn) > 0:
        cto = ccrs.NorthPolarStereo().transform_points(dproj,xtn,ytn)
        interp = intp2D(coords_in[:,0:2], data2d_h)
        data_n = interp(cto[:,0],cto[:,1])
        data_i[indn] = data_n

    return data_i

print ("pjr3.py complete")
#help(findNiceContours)

