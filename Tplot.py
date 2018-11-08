import matplotlib
from mpl_toolkits.basemap import Basemap, cm
# requires netcdf4-python (netcdf4-python.googlecode.com)
from netCDF4 import Dataset as NetCDFFile
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm



nc = NetCDFFile('/Users/d3x345/NetCDF_files/f1850c5_ne30_ne120tunings_JJA_climo.nc')

data = nc.variables['TREFHT']

print "orig shape", data.shape
# matplotlib likes arrays of ny, nx
data = data[0,:,:]
#data = data[0,:,:]

latcorners = nc.variables['lat'][:]
#lats = latcorners
#print "latcorners", latcorners
#print "diffs", np.diff(latcorners)
loncorners = nc.variables['lon'][:]
print "loncorners shape", loncorners.shape
#lons = loncorners
lon_0 = -100.
lat_0 = 45.
# create figure and axes instances
fig = plt.figure(figsize=(8,8))
ax = fig.add_axes([0.1,0.1,0.8,0.8])

m = Basemap(projection='ortho',lon_0=lon_0,lat_0=lat_0)
#m = Basemap(projection='stere',lon_0=lon_0,lat_0=90.,lat_ts=lat_0,\
#            llcrnrlat=latcorners[0],urcrnrlat=latcorners[2],\
#            llcrnrlon=loncorners[0],urcrnrlon=loncorners[2],\
#            rsphere=6371200.,resolution='l',area_thresh=10000)
# draw coastlines, state and country boundaries, edge of map.
m.drawcoastlines()
m.drawstates()
m.drawcountries()
# draw parallels.
#parallels = np.arange(0.,90,10.)
#m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
# draw meridians
meridians = np.arange(180.,360.,10.)
#m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
ny = data.shape[0]; nx = data.shape[1]
print "nx, ny",  nx, ny
#lons, lats = m.makegrid(ny, nx) #return arrays of shape (ny,nx) containing lon,lat coordinates 
lons = np.ones([ny,nx])*loncorners
print "lons.shape", lons.shape
print "lons", lons[0,:]
print "lons[-1]", lons[-1,:]
lats = np.ones([ny,nx])*latcorners.reshape(ny,1)
print "lats.shape", lats.shape
print "lats", lats[:,-1]
x, y = m(lons, lats) # compute map proj coordinates.
print "shape x", x.shape
#quit()
# draw filled contours.
clevs = [0.0000001,0.1,0.2,0.5,1.,2.,3.,4.,5.,10.,50.]
print "clevs", clevs, type(clevs)
#clevs = np.arange(250.,320., 5.).tolist()
#cs = m.contourf(x,y,data,clevs,cmap=cm.s3pcpn)
#cmap = plt.cm.get_cmap("Paired")
#cmap = plt.cm.get_cmap("Dark2_r")
cmap = plt.cm.get_cmap("nipy_spectral",len(clevs))
#colors = cmap(np.linspace(0,1,len(clevs)))
cs = m.contourf(x,y,data,clevs, cmap=cmap, norm=matplotlib.colors.SymLogNorm( linthresh=0.01, linscale=0.01))
#cs = m.contourf(x,y,data,clevs, cmap=cmap)
# add colorbar.
#cbar = m.colorbar(cs,ticks=clevs,location='bottom',pad="5%")
#cbar = m.colorbar(cs,ticks=clevs,location='bottom',format="%.4e",pad="5%")
cbar = m.colorbar(cs,ticks=clevs,location='bottom',pad="5%")
#cbar.set_ticks(np.linspace(0,0.9,0.1))
#cbar.set_ticklabels([str(i) for i in clevs])
cbar.set_label('T(K)')
# add title
#plt.title(precl.long_name+' for period ending '+precl.dateofdata)
plt.show()
