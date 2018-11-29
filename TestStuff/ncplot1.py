import cdms2, cdutil, vcs;

s = vcs.init();
xy = s.createxvsy('test','default')
xy.datawc_y1=75 #your min
xy.datawc_y2=275 # your max
xy.datawc_x1=-90 #xmin
xy.datawc_x2=90. #xmax
xy.linewidth=4
# xy2 = s.getxvsy();
#xy2 = s.createxvsy('test2','default')
##xy2.datawc_y1=75 #your min
#xy2.datawc_y2=275 # your max
#xy2.datawc_x1=-90 #xmin
#xy2.datawc_x2=90. #xmax
##xy2.line = 2;
#xy2.linecolor=22
#xy2.linewidth=4

f = cdms2.open("/Users/d3x345/Desktop/NetCDF_files/COREV2_ncar_rad.15JUNE2009.nc");
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
#s.plot(U);
#raw_input("Press Enter to continue...")

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
#s.plot(U,xy);
#s.plot(V,xy);
lat2;
data2;
s.plot(lat1,data1,xy,long_name='test');
#s.plot(lat2,data2,xy2,long_name='test');
xy.linecolor=22
xy.linewidth=4
s.plot(lat2,data2,xy,long_name='test');
raw_input("Press Enter to continue...")

s.pdf('test_vcs.pdf');
s.png('test_vcs.png');
# s.gif('test_vcs.gif');
