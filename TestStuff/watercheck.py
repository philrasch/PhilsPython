from pjr import *

fname="/pic/scratch/d3x345/model_runs/acme_lowres_compare_fluxes/run/acme_lowres_compare_fluxes.cam.h0.0000-10-01-00000.nc"
fname ="/pic/scratch/d3x345/model_runs/acme_lowres_compare_fluxes/run/out.nc"
f = cdms2.open(fname)

varlist = f.listvariables();
dims = f.listdimension();

V1 = f("QFLXC",time=slice(1,2),squeeze=1)
V2 = f("WPRTP_SFC2",time=slice(1,2),squeeze=1)

V1 = V1*86400.
V1.units = 'mm/day'
V2 = V2*86400.
V2.units = 'mm/day'

LF = f('LANDFRAC',time=slice(1,2),squeeze=1)

AREA = f('area',squeeze=1)
AREA = AREA/AREA.sum()
print "AREA sum", AREA.sum()

LFG = cdutil.averager(LF,weights=AREA)
print "LFG", LFG

V1G = cdutil.averager(V1,weights=AREA)
print "V1G, and range", V1G, V1.min(), V1.max()

V2G = cdutil.averager(V2,weights=AREA)
print "V2G, and range", V2G, V2.min(), V2.max() 

lmax = 30.
lmin = -10.

# Create a Figure object.
fig = plt.figure(figsize=(4, 4))
fig.subplots_adjust(left=0.3)
# Create an Axes object.
ax = fig.add_subplot(1,1,1) # one row, one column, first plot

# Add a title.
ax.set_title("water fluxes (mm/day)")
# Add some axis labels.
ax.set_xlabel("QFLXC")
ax.set_ylabel("WPRTPX")
ax.plot([lmin,lmax],[lmin,lmax],'.r-')

ax.scatter(V1,V2)
ax.set_ylim([lmin,lmax])
ax.set_xlim([lmin,lmax])

#fig.savefig('water_scatter.pdf', format='pdf')
plt.show()

D = V1 - V2
RD = 2*D/(abs(V1)+abs(V2))
print "Relative Error max and min ", RD.max(), RD.min()
print " Error Max and min ", D.max(), D.min()
maxerrind = np.where(D>D.max()*0.99)[0][0]
maxerrind = np.where(D<D.min()*0.99)[0][0]
print "maxerrind = ", maxerrind 
print "D, V1, V2 at maxerrind ", D[maxerrind], V1[maxerrind], V2[maxerrind]

lat = f('lat')
lon = f('lon')
print "coords of maxerr ", lon[maxerrind], lat[maxerrind]
