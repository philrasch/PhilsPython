from pjr import *
import string
from timeit import default_timer as timer

dir1 = "/global/cscratch1/sd/leebrent/climos/output/f.e11.FAMIP.f09_f09_same_setting_LENS/"
f1 = "f.e11.FAMIP.f09_f09_same_setting_LENS_ANN_197912_200511_climo.nc"

jname = dir1+f1

print "jname", jname
g2 = cdms2.open(jname);
varlist = g2.listvariables();
dims = g2.listdimension();
print dims;
print "\n".join(s for s in varlist if 'Z' in s) # find vars that contain 'Z'
#print varlist;

T =  g2('T',squeeze=1)            # extract fields, remove dimensions length one

lat2 = g2['lat']
lon2 = g2['lon']
lev2 = g2['lev']

Txav = cdutil.averager(T,axis="x")
print "Txav=", Txav.info()

ps = g2('PS',squeeze=1)
hyam = g2('hyam',squeeze=1)
hybm = g2('hybm',squeeze=1)
eta = hyam+hybm
print "eta", eta
niceeta = findNiceContours(eta*1010.0,20,rmClev=0.)
print "niceeta", niceeta

start = timer()
P = cdutil.reconstructPressureFromHybrid (ps,hyam,hybm,1.e5)
P = P/100.
P.units = "hPa"
TP2 = cdutil.logLinearInterpolation(T,P,levels=niceeta)
print TP2.info()
# ...
end = timer()
print(end - start) # Time in seconds, e.g. 5.38091952400282print ps.units, ps.max()

TP2xav = cdutil.averager(TP2,axis="x")
levsp = TP2xav.getAxis(0)
#plotZMf(TP2xav, lat2, levsp)
plotZMf(Txav, lat2, levs2)

plt.savefig("test.pdf",type="pdf")
