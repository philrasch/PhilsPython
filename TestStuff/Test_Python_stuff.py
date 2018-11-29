import matplotlib.pyplot as plt
import numpy as np
import cdms2
import cdutil
from scipy.interpolate import interp1d

from mpl_toolkits.basemap import Basemap
#from netCDF4 import Dataset, date2index
from datetime import datetime

# coding: utf-8

# In[1]:

# get_ipython().magic(u'matplotlib inline')
from pjr import *


# In[2]:

fname = "/Users/d3x345/Desktop/NetCDF_files/vd05_ANN_climo.nc"
g2 = cdms2.open(fname);
varlist = g2.listvariables();
dims = g2.listdimension();
#print dims;
#print varlist;

conv = 86400.*1000. # convert m/s to mm/d
PRECC =  g2('PRECC',squeeze=1)*conv            # extract fields, remove dimensions length one
PRECL =  g2('PRECL',squeeze=1)*conv            # extract fields, remove dimensions length one
PRECSC =  g2('PRECSC',squeeze=1)*conv          # extract fields, remove dimensions length one
PRECSL =  g2('PRECSL',squeeze=1)*conv          # extract fields, remove dimensions length one
PRECT = PRECC + PRECL
PRECST = PRECSC + PRECSL
LANDFRAC2 = g2('LANDFRAC',squeeze=1)
print "V2 range", PRECT.min(), PRECT.max()
lat2 = g2['lat']
lon2 = g2['lon']


# if V2av and V2a_av differ it suggests that the averager is using a different weighting value
print "PRECCav=", cdutil.averager(PRECC,axis="xy",weights="weighted")
print "PRECLav=", cdutil.averager(PRECL,axis="xy",weights="weighted")
print "PRECTav=", cdutil.averager(PRECT,axis="xy",weights="weighted")
print "PRECSCav=", cdutil.averager(PRECSC,axis="xy",weights="weighted")
print "PRECSLav=", cdutil.averager(PRECSL,axis="xy",weights="weighted")
print "PRECSTav=", cdutil.averager(PRECST,axis="xy",weights="weighted")

print "LANDFRAC2av=", cdutil.averager(LANDFRAC2,axis="xy",weights="weighted")


# In[ ]:




# In[3]:

plt.figure(figsize=(16, 4))
#fig, ax = plt.subplots()
plt.subplot(121)
plt.title("land fraction")
cs = plt.contourf(lon2, lat2, LANDFRAC2)
# add colorbar.
cbar = plt.colorbar(cs)
cbar.set_label('land(Fraction)')
print "LF2a_av=", cdutil.averager(LANDFRAC2,axis="xy",weights="weighted")
plt.show()


# In[4]:

#LANDFRAC.info()
#print "LANDFRAC2 =", cdutil.averager(LANDFRAC2,axis="xy",weights="weighted")
#print "masked LANDFRAC2 =", cdutil.averager(LFM,axis="xy",weights="weighted")
#print "dir(LF)", dir(LANDFRAC)
#print "xx", LANDFRAC2.info()
llist = (LANDFRAC2[:] < 0.5)
print "llist",llist[0,0],llist[-1,-1] # print out a few values which should be near NP and SP
#print "land zero", llist.shape
DVM = LANDFRAC2
DVM.mask = llist # add a mask for land points
print "longname",DVM.long_name
plt.figure(figsize=(16, 4))
#fig, ax = plt.subplots()
plt.subplot(121)
gavg = cdutil.averager(DVM,axis="xy",weights="weighted")
#titl = "{0:s} (model-obs) = {1:f}".format("DVM",gavg)
titl = "%s = %.2f"%("DVM", gavg)
print titl
plt.title(titl)
clevs =  np.arange(-60.,70., 10).tolist()
print clevs
cs = plt.contourf(lon2, lat2, DVM, clevs)
# add colorbar.
cbar = plt.colorbar(cs)
plt.show()
#print "masked DV =", cdutil.averager(DVM,axis="xy",weights="weighted")



# In[5]:

# example from http://matplotlib.org/basemap/users/examples.html

date = datetime(2007,12,15,0) # date to plot.
# open dataset.
#dataset = \
#Dataset('http://www.ncdc.noaa.gov/thredds/dodsC/OISST-V2-AVHRR_agg')
#timevar = dataset.variables['time']
#timeindex = date2index(date,timevar) # find time index for desired date.
# read sst.  Will automatically create a masked array using
# missing_value variable attribute. 'squeeze out' singleton dimensions.
#sst = dataset.variables['sst'][timeindex,:].squeeze()
# read ice.
#ice = dataset.variables['ice'][timeindex,:].squeeze()

lons, lats = np.meshgrid(lon2,lat2)
print lons.shape
# create figure, axes instances.
fig = plt.figure()
ax = fig.add_axes([0.05,0.05,0.9,0.9])
# create Basemap instance.
# coastlines not used, so resolution set to None to skip
# continent processing (this speeds things up a bit)
m = Basemap(projection='kav7',lon_0=0,resolution=None)
# draw line around map projection limb.
# color background of map projection region.
# missing values over land will show up this color.
m.drawmapboundary(fill_color='0.3')
# plot sst, then ice with pcolor
# im1 = m.pcolormesh(lons,lats,DVM,shading='flat',cmap=plt.cm.jet,latlon=True)
# cb = m.colorbar(im1,"bottom", size="5%", pad="2%")
# draw parallels and meridians, but don't bother labelling them.
# contour data over the map.
cs = m.contourf(lons,lats,DVM,10,latlon=True)
plt.title('contour lines over filled continent background')
m.drawparallels(np.arange(-90.,99.,30.))
m.drawmeridians(np.arange(-180.,180.,60.))
# add colorbar
cb = m.colorbar(cs,"bottom", size="5%", pad="2%")
# add a title.
ax.set_title('DVM %s'% DVM.long_name)
plt.show()

