'''
Long term archiving for E3SMv2, CESM2 simulations on AWS

Loads E3SMv2, CESM2 history fields into xarray datasets and
writes them to Zarr store on s3 storage.
'''

import glob,json,zarr,fsspec,dask,cftime
import xarray as xr
import s3fs
import datetime
from dask.distributed import Client
from datetime import timedelta

class hist_to_zarr():
    '''
    #############################################################################
    
    Reshape on-disk netCDF history fields to zarr stores
    
    #############################################################################
    '''
    def __init__(self,
                source_dir:str,
                dest_dir:str,
                experiment:str,
                f_domain_info:str = 'domain_info.json',
                f_var_to_process:str = None,
                verbose:bool = True,
                log_to_file: bool = True,
                s3_region:str = 'us-east-2',
                nworker:int = 8,
                nthread:int = 1,
                combine_1d_var = True,
                filesuffix = '.zarr'):
        '''
        source_dir : path to source model archive on disk
        dest_dir : path to destination reshape location on disk (must start with s3://) (AWS)
        experiment : experiment label
        f_domain_info : path to a json file containing a **list** of information of data chunks and frequencies
                    each list entry should hold: [COMPONENT:STR, MODEL:STR, FREQ:STR, FREQNAME:STR, CHUNKING:DICT, TIME DIMENSION:STR, TIME VARIABLES:LIST]
                    CHUNKING dicionary is used by (XARRAY) 
        f_var_to_process (optional) : path to a json file containing a dictionary
                                    with variables to reshape
        verbose : toggles various print statements
        log_to_file : toggle writing reshaping to do list to file
        s3_region : AWS location of s3 bucket (AWS)
        nworker : number of workers (DASK)
        nthread : number of threads per worker (DASK)
        '''
        
        self.source_dir = source_dir
        self.dest_dir = dest_dir
        self.experiment = experiment
        self.verbose = verbose
        self.log_to_file = log_to_file
        self.domain_info = json.load(open(f_domain_info))
        self.s3_region = s3_region
        self.nworker = nworker
        self.nthread = nthread
        self.f_var_to_process = f_var_to_process
        self.combine_1d_var = combine_1d_var
        self.filesuffix = filesuffix
        
        self.max_jobs_inqueue = 128
        
        self.fs = s3fs.S3FileSystem(anon=False)

        assert self.dest_dir[:5] == 's3://'
        assert type(self.domain_info) is list
        
        if log_to_file:
            now = datetime.datetime.now()
            nowstr = now.strftime("%Y_%m_%d-%H:%M:%S")
            self.log_fname = f'zarr_reshape.{experiment}.{nowstr}.log'
            print(self.log_fname)

    def run(self):
        self.setup_reshape()
        reshape_obj.run_processing()
        reshape_obj.append_metadata_var()
        
    def setup_reshape(self):
        '''
        Retrieves/reads "to-do list" of variables to reshape
        '''
        if self.f_var_to_process is None:
            self.to_do_dict = self.get_var_list()
        else:
            self.to_do_dict = json.load(open(self.f_var_to_process))
            print(self.to_do_dict)
            self.const_vars = {key:[] for key in self.to_do_dict}
            self.oneD_vars= {key:[] for key in self.to_do_dict}
            self.with_errors={key:[] for key in self.to_do_dict}
            self.complete_variables= {key:[] for key in self.to_do_dict}

        if self.log_to_file:
            self.log_todo()

    def get_first_last_timestamp(self,path,model,timedim):
        '''
        looks at archive directory to determine the first and last time
        stamp to be processed
        '''
        # Grabs first file for the given component, frequency
        file_list = sorted(list(glob.glob(f'{path}')))
        first_fname = file_list[0]
        last_fname = file_list[-1]

        # load the first and last raw history field files
        data_first = xr.open_dataset(first_fname)
        data_last = xr.open_dataset(last_fname)
        
        if model == 'mpaso' or model == 'mpassi':
            data_first = e3sm_time_conversion(data_first)
            data_last = e3sm_time_conversion(data_last)

        # Get first and last time steps (to compare vs. zarr store)
        time_i = data_first[timedim][0]
        time_f = data_last[timedim][-1]
        return data_first, time_i, time_f
    
    def get_var_list(self, check_existing = True):
        '''
        Retrieves to do list of reshaping jobs
        '''

        if self.verbose:
            print('#############################################################################')
            print('Generating reshaping to do list')
            print('#############################################################################')

        to_do = {}
        self.const_vars = {}
        self.oneD_vars = {}
        self.with_errors = {}
        self.complete_variables = {}
        for domain, model, freq, freq_name, __, timedim, timevars in self.domain_info:
            if self.verbose: print(domain, model, freq)

            path = f'{self.source_dir}/{domain}/hist'
            __,prefix = get_prefix(self.experiment,path,model,freq)
            if self.verbose: print(path, prefix)

            # Grabs first file for the given component, frequency
            data, time_i, time_f = self.get_first_last_timestamp(f'{path}/{prefix}',model,timedim)
            
            if self.verbose: print(time_i, time_f)

            # list of variables in the history files
            init_varlist = list(data.variables)

            # generating to do list
            # EXCLUDING:
            # - Non-time varying variables (i.e, gw)
            # - 1D variables (i.e., co2_vmr)
            # - variables with complete records
            const_vars = []
            oneD_vars = []
            to_do_vars = []
            with_errors = []
            complete_variables = []

            zarr_path = f'{self.dest_dir}/{self.experiment}/{domain}/proc/tseries/{freq_name}/'

            for var in init_varlist:
                if self.verbose: print('Checking: ',var)

                if not timedim in data[var].dims:
                    const_vars.append(var)
                    if self.verbose: print('Constant variable, skipping')
                elif (timedim in data[var].dims and len(data[var].dims) == 1) or \
                    var in timevars:
                    oneD_vars.append(var)
                    if self.verbose: print('one-dimensional, skipping')
                elif var in ['time_bnds','time_bounds','time_bound']:
                    oneD_vars.append(var)
                    if self.verbose: print('one-dimensional, skipping')
                elif check_existing:
                    ## checks for existing zarr stores
                    ## then checks if the current time steps exist
                    zarr_store = zarr_path + f'{self.experiment}.{model}.{freq}.{var}{self.filesuffix}'
                    if self.fs.isdir(zarr_store):
                        if self.verbose: print('Existing Zarr store found, getting time information')
                        if not self.fs.isdir(zarr_store + f'/{var}'):
                            with_errors.append(zarr_store)
                            print(f'Zarr store error, skipping: {zarr_store}')
                            continue;

                        # load zarr store
                        zarr_data = load_s3_zarr(zarr_store,self.s3_region)

                        # checks final time stamp in zarr and compare to
                        # netcdf dataset
                        t0 = zarr_data[timedim][-1]
                        print(t0)
                        
                        if zarr_data[var].shape[0] != zarr_data[timedim].shape[0]:
                            with_errors.append(zarr_store)
                            print(f'Zarr store error, skipping: {zarr_store}')
                            continue;
                        elif t0 == time_f:
                            print(f'Variable complete, skipping {var}')
                            print(f'{zarr_store}')
                            complete_variables.append((var,zarr_store))
                            continue;
                        else:
                            if self.verbose: print(f'Append to: {zarr_store}')
                            mode = 'a'
                            t0 = t0 + timedelta(days=1)
                            varlist = [var for var in zarr_data.variables
                                              if timedim in zarr_data[var].dims]
                        zarr_data.close()
                    else:
                        if self.verbose: print(f'New file write to: {zarr_store}')
                        if self.verbose: print(f'{data[var].dims}')
                        mode = 'w'
                        t0 = '0000-01-01'
                        varlist = [var]
                    to_do_vars.append((var,zarr_store,mode,t0,varlist))
                else:
                    zarr_store = zarr_path + f'{self.experiment}.{model}.{freq}.{var}{self.filesuffix}'
                    
                    if self.verbose: print(f'New file write to: {zarr_store}')
                    if self.verbose: print(f'{data[var].dims}')
                    
                    mode= 'w'
                    t0 = '0000-01-01'
                    to_do_vars.append((var,zarr_store,mode,t0,[var]))

                  
            to_do[f'{domain}.{model}.{freq}'] = to_do_vars
            self.const_vars[f'{domain}.{model}.{freq}'] = const_vars
            self.oneD_vars[f'{domain}.{model}.{freq}'] = oneD_vars
            self.with_errors[f'{domain}.{model}.{freq}'] = with_errors
            self.complete_variables[f'{domain}.{model}.{freq}'] = complete_variables
            
        return to_do
    
    def dask_compute(self, l_jobs):
        ## Starts the dask processes!
        print('BEGINNING DASK DISTRIBUTED COMPUTE')
        print('This will take a while...')
        with Client(threads_per_worker=self.nthread, n_workers=self.nworker) as client:
            dask.compute(*l_jobs)

    def run_processing(self):
        for domain, model, freq, freq_name, chunk, timedim, timevars in self.domain_info:
            print(f'{domain}.{model}.{freq}')
            
            if len(self.to_do_dict[f'{domain}.{model}.{freq}']) == 0:
                print(f'{domain} {model} {freq}: Completed, skipping.')
                continue;
                
            path = f'{self.source_dir}/{domain}/hist'
            file_list,prefix = get_prefix(self.experiment,path,model,freq)
            
            # Load data
            data = xr.open_mfdataset(f'{path}/{prefix}',combine='nested',concat_dim = timedim)
            
            ## MPAS TIME STAMPING
            if model == 'mpaso' or model == 'mpassi':
                data = e3sm_time_conversion(data)
                
            ## Sets up DASK Delayed objects
            ## Creates empty zarr stores on s3 at this stages
            l_pool_delayed = []
            #print(self.to_do_dict[f'{domain}.{model}.{freq}'])
            for var,zarr_store,mode,t0,varlist in self.to_do_dict[f'{domain}.{model}.{freq}']:
                if self.verbose: print(f'Set up DASK delayed: {var}')
                if mode == 'a':
                    in_data = data[varlist].chunk(chunk).sel({timedim:slice(t0,'9999-12-31')})
                elif mode == 'w':
                    in_data = data[varlist + timevars].chunk(chunk).sel({timedim:slice(t0,'9999-12-31')})

                mapper = fsspec.get_mapper(zarr_store,client_kwargs={'region_name':self.s3_region},check=False)
                
                if mode == 'a':
                    writeobj = in_data.to_zarr(store=mapper,mode=mode,compute=False,append_dim = timedim)
                elif mode == 'w':
                    writeobj = in_data.to_zarr(store=mapper,mode=mode,compute=False)
                l_pool_delayed.append(writeobj)
                
                if len(l_pool_delayed) > self.max_jobs_inqueue:
                    self.dask_compute(l_pool_delayed)
                    del l_pool_delayed
                    l_pool_delayed = []
                    
            if len(l_pool_delayed) > 0 :
                self.dask_compute(l_pool_delayed)

            ### combines 1D variables into a single store to avoid too many small files
            varlist = [var for var in self.oneD_vars[f'{domain}.{model}.{freq}']
                      if not var in timevars]

            ## Reshape time series files to zarr
            if len(varlist) > 0 and self.combine_1d_var:
                print("COMBINE 1D TIME SERIES FILES")
                zarr_path = f'{self.dest_dir}/{self.experiment}/{domain}/proc/tseries/{freq_name}'
                zarr_store = f'{zarr_path}/{self.experiment}.{model}.{freq}.1D_vars{self.filesuffix}'
                if self.fs.isdir(zarr_store):
                    mode_1d = 'a'
                else:
                    mode_1d = 'w'
                mapper = fsspec.get_mapper(zarr_store,
                                               client_kwargs={'region_name':self.s3_region},
                                               check=False)

                in_data = data[varlist]
                if mode_1d == 'a':
                    writeobj = in_data.to_zarr(store=mapper,mode=mode_1d,compute=True,append_dim = timedim)
                else:
                    writeobj = in_data.to_zarr(store=mapper,mode=mode_1d,compute=True)

        return 1

    def clean_with_errors(self):
        for domain, model, freq, freq_name, chunk, timedim, timevars in self.domain_info:
            for zarr_store in self.with_errors[f'{domain}.{model}.{freq}']:
                print(f'Clean {zarr_store}')
                #sp.run(f'aws s3 rm --recursive {zarr_store}', shell = True)
    
    def append_const_to_zarr(self,refdata,var,zarr_store,domain,model,freq,timedim):
        # checks to see if metadata variables already exists
        print(f'Checking: {zarr_store}')
        data = load_s3_zarr(zarr_store,self.s3_region)
        append_const = []

        print(list(data.variables))
        for const in self.const_vars[f'{domain}.{model}.{freq}']:
            if not const in list(data.variables):
                append_const.append(const)

        data.close()

        if len(append_const) > 0:
            if self.verbose: print(f'Append to : {zarr_store}')
            if self.verbose: print(append_const)

            mapper = fsspec.get_mapper(zarr_store,
                                   client_kwargs={'region_name':self.s3_region},
                                   check=False)

            outdata = refdata[append_const]
            writeobj = outdata.to_zarr(store=mapper,mode='a',compute=False)
            
            return writeobj
        else:
            return 0
    
    def append_metadata_var(self,include_complete = False):
        ### Appending constant (metadata) variables
        for domain, model, freq, freq_name, chunk, timedim, timevars in self.domain_info:
            path = f'{self.source_dir}/{domain}/hist'
            i = 0 
            file_list,prefix = get_prefix(self.experiment,path,model,freq)
            
            first_fname = file_list[0]
            refdata = xr.open_dataset(first_fname)

            print('APPEND CONSTANT VARIABLES')
            l_pool_delayed = []
            for var,zarr_store,__,__,__ in self.to_do_dict[f'{domain}.{model}.{freq}']:
                writeobj = self.append_const_to_zarr(refdata,var,zarr_store,domain,model,freq,timedim)
                if not type(writeobj) is int:
                    l_pool_delayed.append(writeobj)
            for var,zarr_store in self.complete_variables[f'{domain}.{model}.{freq}']:
                writeobj = self.append_const_to_zarr(refdata,var,zarr_store,domain,model,freq,timedim)
                if not type(writeobj) is int:
                    l_pool_delayed.append(writeobj)
            if len(l_pool_delayed) == 0:
                continue;
            with Client(threads_per_worker=self.nthread, n_workers=self.nworker) as client:
                dask.compute(*l_pool_delayed)

                
    '''
    #############################################################################
    Logging files
    #############################################################################
    '''
    def log_todo(self):
        with open(self.log_fname,'w') as file:
            file.write('----------------------------------------------------------------\n')
            file.write(f'File conversion logs for {self.experiment} \n')
            file.write(f'Source directory: {self.source_dir} \n')
            file.write(f'Destination directory: {self.dest_dir}/{self.experiment} \n')
            file.write('----------------------------------------------------------------\n')

            file.write('\n')
            file.write('Summary : \n\n')
            for domain, model, freq, freq_name, __, timedim, timevars in self.domain_info:
                file.write(f'{domain}, {model}, {freq} \n')
                l_complete_variables = self.complete_variables[f'{domain}.{model}.{freq}']
                l_todo = self.to_do_dict[f'{domain}.{model}.{freq}']
                if len(self.to_do_dict[f'{domain}.{model}.{freq}']) == 0:
                    file.write(f'Complete : {len(l_complete_variables)}\n')
                    file.write('ALL COMPLETE \n')
                else:
                    file.write(f'Complete : {len(l_complete_variables)}\n')
                    file.write(f'Processing : {len(l_todo)}\n')
                    
            file.write('----------------------------------------------------------------\n')
            file.write('\n')
            file.write('List of files to write : \n\n')
            
            for domain, model, freq, freq_name, __, timedim,timevars in self.domain_info:
                file.write(f'{domain}, {model}, {freq}')
                l_todo = self.to_do_dict[f'{domain}.{model}.{freq}']
                for fname in l_todo:
                    file.write(f'{fname}\n')
                file.write('----------------------------------------------------------------\n')
                file.write('----------------------------------------------------------------\n')
        return self.log_fname
        
'''
Helper functions
'''

def load_s3_zarr(zarr_store,s3_region):
    # opens zarr store
    mapper = fsspec.get_mapper(zarr_store,
                        client_kwargs={'region_name':s3_region},
                        check=False)
    zarr_data = xr.open_dataset(mapper,engine = 'zarr')
    return zarr_data

def get_prefix(experiment,path,model,freq):
    # Have to insert in first number of year into filename wildcard
    # b/c of pop file naming convention
    i = 0 
    file_list = []
    while len(file_list) == 0 and i < 10:
        prefix = f'{experiment}.{model}.{freq}.{i}*'
        file_list = sorted(list(glob.glob(f'{path}/{prefix}')))
        i += 1
    return file_list,prefix

def e3sm_time_conversion(data):
    ### hacky way of converting E3SMv2 mpas time into cftime values
    ### so it can be processed using xarray
    ### assumes xtime is a string with values "b'YYYY-MM-DD_TT:TT:TT'"
    for vv in ['xtime','xtime_startMonthly','xtime_startMonthlyMax','xtime_startMonthlyMin',
              'xtime_startDaily']:
        if vv in data.variables:
            timevar = vv
    l_cftime = [cftime.DatetimeNoLeap(*[int(val) 
                    for val in str(t.values)[2:].split('_')[0].split('-')]) 
                             for t in data[timevar]]
    data = data.assign_coords({'time':('Time',l_cftime)})
    data = data.set_index(Time = "time")
    return data