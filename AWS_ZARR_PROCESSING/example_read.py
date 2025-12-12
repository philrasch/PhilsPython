import xarray as xr
import fsspec
import s3fs

fs = s3fs.S3FileSystem(anon=False)

# Loads a single Zarr file
def s3load_zarr(s3path:str):
    return xr.open_dataset(fsspec.get_mapper(s3path),engine = 'zarr')

# Queries a directory with all the zarr stores for a given simulation/domain/frequency
# and looks for the zarr stores matching the list of variables given
def multi_s3load_zarr(parent_path:str,exp:str,domain:str,freq:str,variables:list):
    s3path = f'{parent_path}/{exp}/{domain}/proc/tseries/{freq}'
    to_load = []
    for item in fs.listdir(s3path):
        path = item['Key']
        for var in variables:
            if f'.{var}.' in path:
                to_load.append(xr.open_dataset(fsspec.get_mapper('s3://'+path),engine = 'zarr').drop_duplicates("time"))    
    return xr.merge(to_load,compat = 'override')

data = multi_s3load_zarr('s3://mcb-zarr/CESM2','BSSP245.G6-1.5K-MCB.massbal.001','atm','month_1',['TREFHT','PRECT','QFLX','ncl_a1SF'])
#data['time'] = data.indexes['time'].shift(-15,"D")

data.mean("time")