# module of phils useful stuff
"""Phils useful functions"""

import matplotlib.pyplot as plt
import numpy as np
import cdms2
import cdutil
from scipy.interpolate import interp1d

from mpl_toolkits.basemap import Basemap
#from netCDF4 import Dataset, date2index
from datetime import datetime
import matplotlib.colors as mcolors
def make_colormap(seq):
    """Return a LinearSegmentedColormap
    seq: a sequence of floats and RGB-tuples. The floats should be increasing
    and in the interval (0,1).
    """
    seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
    print "seq",seq
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
    print "class", data2.__class__
    print "lat2.size", np.size(lat2)
    print "aa"
    print "lat1.size", np.size(lat1)
    print "bb"
    plt.rc('axes', color_cycle=['r', 'g', 'b', 'y'])
    x2 = np.array(lat2);
    y2 = np.array(data2);
    x1 = np.array(lat1)
    y1 = np.array(data1)
    l1, = ax0.plot(x1, y1,label='CORE.v2')
    l2, = ax0.plot(x2, y2,label='EBAFS-Sfc')
    print "cc"
    ax0.legend(handles=[l1,l2])
    #ax0.plot(lat2,data2)
    ax0.set_title('Downward Shortwave')
    print "xxx", np.size(lat2),  np.size(lat1)
    if np.size(lat2) > np.size(lat1):
        xh = x2;
        y1 = np.interp(xh, x1, y1)
    else:
        xh = x1;
        y2 = np.interp(xh, x2, y2)
        print "regrid data2.size", np.size(y2)
        print "regrid data1.size", np.size(y1)
    dy = np.array(y2-y1);
    ax1.plot(xh, dy);
        
