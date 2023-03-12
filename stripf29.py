import argparse
import numpy as np
import xarray as xr

parser = argparse.ArgumentParser()
parser.add_argument("infile", help="inputfile", type=str)
parser.add_argument("outfile",help="outputfile", type=str)
args = parser.parse_args()
infile = args.infile
outfile = args.outfile
print(infile,outfile)

DS = xr.open_mfdataset(infile).chunk({'time': 200})

Var = DS.ts

# identify times that are Feb 29
f29 = ((Var.time.dt.month == 2) & (Var.time.dt.day == 29))
if len(np.where(f29)[0]) > 0:
    myinds = np.where(~f29)[0] # select times that are not feb29
    Var_nof29 = Var.isel(time=myinds)
    Var_nof29.to_netcdf(outfile)
else:
    print('no feb 29 on the file')

