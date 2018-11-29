from mpl_toolkits.basemap import Basemap, cm
# requires netcdf4-python (netcdf4-python.googlecode.com)
from netCDF4 import Dataset as NetCDFFile
import numpy as np
import matplotlib.pyplot as plt

# plot rainfall from NWS using special precipitation
# colormap used by the NWS, and included in basemap.

# nc = NetCDFFile('/Users/d3x345/ACME_NetCDF/for_phil/f1850c5_conusx4v1/f1850c5_conusx4v1_JJA_climo.nc')
nc = NetCDFFile('/Users/d3x345/ACME_NetCDF/for_phil/f1850c5_conusx4v1/f1850c5_conusx4v1_JJA_climo.nc')

precl = nc.variables['PRECL']
precc = nc.variables['PRECC']
# convert from m/s to mm/day

data = 8.64e7*(precl[:]+precc[:])
print data.shape
#data = data[0,:,:].transpose()
data = data[0,:,:]
print data.shape
latcorners = nc.variables['lat'][:]
lats = latcorners
#print "latcorners", latcorners
#print "diffs", np.diff(latcorners)
loncorners = -nc.variables['lon'][:]
lons = loncorners
lon_0 = 0.
lat_0 = 0.
# create figure and axes instances
fig = plt.figure(figsize=(8,8))
ax = fig.add_axes([0.1,0.1,0.8,0.8])
# create polar stereographic Basemap instance.
m = Basemap(projection='stere',lon_0=lon_0,lat_0=lat_0,llcrnrlat=-30,urcrnrlat=30,llcrnrlon=-80,urcrnrlon=+80)
#m = Basemap(projection='stere',lon_0=lon_0,lat_0=90.,lat_ts=lat_0,\
#            llcrnrlat=latcorners[0],urcrnrlat=latcorners[2],\
#            llcrnrlon=loncorners[0],urcrnrlon=loncorners[2],\
#            rsphere=6371200.,resolution='l',area_thresh=10000)
# draw coastlines, state and country boundaries, edge of map.
m.drawcoastlines()
m.drawstates()
m.drawcountries()
# draw parallels.
parallels = np.arange(0.,90,10.)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
# draw meridians
meridians = np.arange(180.,360.,10.)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
ny = data.shape[0]; nx = data.shape[1]
print "nx, ny",  nx, ny
lons, lats = m.makegrid(nx, ny) # get lat/lons of ny by nx evenly space grid.
print "lons.shape", lons
x, y = m(lons, lats) # compute map proj coordinates.
print "shape x", x.shape
# draw filled contours.
clevs = [0,1,2.5,5,7.5,10,15,20,30,40,50,70,100,150,200,250,300,400,500,600,750]
cs = m.contourf(x,y,data,clevs,cmap=cm.s3pcpn)
# add colorbar.
cbar = m.colorbar(cs,location='bottom',pad="5%")
cbar.set_label('mm')
# add title
#plt.title(precl.long_name+' for period ending '+precl.dateofdata)
plt.show()
