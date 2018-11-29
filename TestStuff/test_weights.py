import matplotlib.pyplot as plt
import numpy as np
import cdms2
import cdutil

dirname = "/pic/scratch/harr152/climo/20161118.beta0.F20TRCOSP.ne30_ne30.edison/reg_climo/yrs1975_2000/"
fname = "20161118.beta0.F20TRCOSP.ne30_ne30.edison_ANN_climo.nc"

jname = dirname+fname

print "jname"
print  jname

f = cdms2.open(jname)

varlist = f.listvariables();
dims = f.listdimension();

V2 = f("TS",squeeze=1)
LF = f('LANDFRAC',squeeze=1)

AREA = f('area',squeeze=1)
AREA = AREA/AREA.sum()
print "AREA sum", AREA.sum()

LFG = cdutil.averager(LF,axis='xy',weights=AREA)
print "global LFG should be about 0.29 ", LFG

# a variable that is all ones
UNITY = V2/V2
# a variable that is proportional to land fraction
UL = UNITY*LF
# global average of one
UG =  cdutil.averager(UNITY,axis='xy',weights=AREA)
print "UG, (should be one)  and range", UG, UNITY.min(), UNITY.max() 
# average of land fraction over whole globe
ULGL1 =  cdutil.averager(UL,axis='xy',weight=AREA)
print "ULGL1 (should be about 0.29) ", ULGL1
# same using weights arg
ULGL2 =  cdutil.averager(UL,axis='xy',weights=AREA)
print "ULGL2 (should also be about 0.29) ", ULGL2
# average of land fraction over land (should be one)
ULGL3 =  cdutil.averager(UL,axis='xy',weight=AREA*LF)
print "ULGL3 should be about 1", ULGL3
ULGL4 =  cdutil.averager(UL,axis='xy',weights=AREA*LF)
print "ULGL4 should be about 1", ULGL4

