import cdms2, cdutil
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# f = cdms2.open("/Users/d3x345/Desktop/NetCDF_files/COREV2_ncar_rad.15JUNE2009.nc");
f = cdms2.open("/pic/scratch/d3x345/COREV2_ncar_rad.15JUNE2009.nc");
varlist = f.listvariables();
dims = f.listdimension();
print dims;

print varlist;

U =  f('SWDN');
U = cdutil.averager(U, axis='x');
U = cdutil.averager(U, axis='0',weights='equal');
data1 = U[...]
lat1 = f['LAT'];
print "lat1", lat1.shape;
print "data ", data1.shape;

g = cdms2.open("/Users/d3x345/Desktop/NetCDF_files/CERES_EBAF-Surface_Ed2.8_Subset_200003-201409.nc");
varlist = g.listvariables();
dims = g.listdimension();
print dims;

print varlist;

V =  g('sfc_sw_down_all_mon');
V = cdutil.averager(V, axis='x');
V = cdutil.averager(V, axis='0',weights='equal');
data2 = V[...]
lat2 = g['lat'];
print "lat2", lat2.shape;
print "data2", data2.shape;


plt.rc('lines', linewidth=4)
fig, (ax0)  = plt.subplots(nrows=1)
lat2.info()
data2.info()

plt.rc('axes', color_cycle=['r', 'g', 'b', 'y'])
x2 = np.array(lat2);
y2 = np.array(data2);
x1 = np.array(lat1)
y1 = np.array(data1)
l1, = ax0.plot(x1, y1,label='CORE.v2')
l2, = ax0.plot(x2, y2,label='EBAFS-Sfc')
ax0.legend(handles=[l1,l2])
#ax0.plot(lat2,data2)
ax0.set_title('Downward Shortwave')
if np.size(lat2) > np.size(lat1):
    xh = x2;
    y1 = np.interp(xh, x1, y1)
else:
    xh = x1;
    y2 = np.interp(xh, x2, y2)
print "regrid data2.size", np.size(y2)
print "regrid data1.size", np.size(y1)
#dy = np.array(y2-y1);
#ax1.plot(xh, dy);

plt.savefig('test_matplotlib.pdf', format='pdf')
plt.savefig('test_matplotlib.png', format='png')


