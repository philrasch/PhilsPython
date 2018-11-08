from pjr import *
from matplotlib.backends.backend_pdf import PdfPages

dirname = "/pic/scratch/harr152/climo/20161118.beta0.F20TRCOSP.ne30_ne30.edison/reg_climo/yrs1975_2000/"
fname = "20161118.beta0.F20TRCOSP.ne30_ne30.edison_ANN_climo.nc"

jname = dirname+fname

print "jname"
print  jname

f = cdms2.open(jname)
print f

varlist = f.listvariables();
dims = f.listdimension();

V1 = f("TS",squeeze=1)
V2 = f("TS",squeeze=1)

#V1 = V1
#V1.units = 'mm/day'
V2 = V2-273.16
V2.units = 'Degree C'

LF = f('LANDFRAC',squeeze=1)

AREA = f('area',squeeze=1)
AREA = AREA/AREA.sum()
print "AREA sum", AREA.sum()

LFG = cdutil.averager(LF,axis='xy',weights=AREA)
print "LFG", LFG

V1G = cdutil.averager(V1,axis='xy',weights=AREA)
print "V1G, and range", V1G, V1.min(), V1.max()

V2G = cdutil.averager(V2,axis='xy',weights=AREA)
print "V2G, and range", V2G, V2.min(), V2.max() 

#print "V2.info", V2.info()

# a variable that is all ones
UNITY = V2/V2
# a variable that is proportional to land fraction
UL = UNITY*LF
# global average of one
UG =  cdutil.averager(UNITY,axis='xy',weights=AREA)
print "UG, and range", UG, UNITY.min(), UNITY.max() 
# average of land fraction over whole globe
ULG =  cdutil.averager(UL,axis='xy',weight=AREA)
# the next line give the wrong answer
#ULG =  cdutil.averager(UL,axis='xy',weights=['unweighted','unweighted'])
print "ULG, and range", ULG, UL.min(), UL.max() 
# average of land fraction over land (should be one)
ULGL =  cdutil.averager(UL,axis='xy',weight=AREA*LF)
print "ULGL, and range", ULGL, UL.min(), UL.max() 
quit()

V3 = V2*LF

lmax = 30.
lmin = -10.

# Create a Figure object.
fig = plt.figure(figsize=(4, 4))
fig.subplots_adjust(left=0.3)
# Create an Axes object.
ax = fig.add_subplot(1,2,1) # one row, one column, first plot

# Add a title.
ax.set_title("water fluxes (mm/day)")
# Add some axis labels.
ax.set_xlabel("TS")
ax.set_ylabel("TS")
ax.plot([lmin,lmax],[lmin,lmax],'.r-')

ax.scatter(V1,V2)
#ax.set_ylim([lmin,lmax])
#ax.set_xlim([lmin,lmax])

plt.show()
fig.savefig('scatter.png', format='png')
plt.close() # close figure window
plt.clf() # clear figure
#plt.show()

D = V1 - V2
RD = 2*D/(abs(V1)+abs(V2))
print "Relative Error max and min ", RD.max(), RD.min()
print " Error Max and min ", D.max(), D.min()
maxlist  = zip(*np.where(D==D.max()))
minlist  = zip(*np.where(D==D.min()))
print "points where D is a maximum, minimum = ", len(maxlist), len(minlist)

lat = f['lat']
lon = f['lon']
#print "coords of maxerr ", lon[maxerrind], lat[maxerrind]
lat = f['lat']
cs = plt.contourf(lon, lat, V3)
plt.contour(lon, lat, LF,levels=[0.5])
cbar = plt.colorbar(cs)
plt.show()
plt.savefig('latlon.png', format='png')
