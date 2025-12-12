#!/opt/ncar/home/haruki/.conda/envs/hhenv/bin python3
# Batch system directives
# SBATCH  --job-name=lt_zarr_archive.20240425.LR.WCYCLSSP245.MCB-SSLT-EM.NI
# SBATCH  --nodes=1
# SBATCH  --output=lt_zarr_archive.20240425.LR.WCYCLSSP245.MCB-SSLT-EM.NI.%j 
# SBATCH  --exclusive

import zarr_reshape
case = 'b.e21.BSSP245smbb_MCBss2.5TgYr_NA.f09_g17.LE2-1011.001'
sourcepath = f'/mcbcesm/ec2-user/archive/{case}/'
s3path = 's3://e3sm-reshaped/zarr-conv/CESM2'
if __name__ == '__main__':
    reshape_obj = zarr_reshape.hist_to_zarr(sourcepath,
                                            s3path,
                                            case,
                                            f_domain_info = 'domain_info.CESM2.json')
    reshape_obj.setup_reshape()
    reshape_obj.run_processing()
    print(reshape_obj.const_vars)
    reshape_obj.append_metadata_var()
