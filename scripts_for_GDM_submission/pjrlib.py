def make_colormap(seq):
    """Return a LinearSegmentedColormap
    seq: a sequence of floats and RGB-tuples. The floats should be increasing
    and in the interval (0,1).
    """
    seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
#   print "seq",seq
    cdict = {'red': [], 'green': [], 'blue': []}
    for i, item in enumerate(seq):
        if isinstance(item, float):
            r1, g1, b1 = seq[i - 1]
            r2, g2, b2 = seq[i + 1]
            cdict['red'].append([item, r1, r2])
            cdict['green'].append([item, g1, g2])
            cdict['blue'].append([item, b1, b2])
    return mcolors.LinearSegmentedColormap('CustomMap', cdict)

def diverge_map(high=(1., 0., 0.), low=(0.0, 0., 1.)):
#def diverge_map(high=(0.565, 0.392, 0.173), low=(0.094, 0.310, 0.635)):
    '''
    low and high are colors that will be used for the two
    ends of the spectrum. they can be either color strings
    or rgb color tuples
    '''
    c = mcolors.ColorConverter().to_rgb
    if isinstance(low, str): low = c(low)
    if isinstance(high, str): high = c(high)
    return make_colormap([low, c('white'), 0.5, c('white'), high])

def findNiceContours(data,nlevs=None,rmClev=None,sym=None,verbose=None):
    """Find Nice Contours
    data = 2d numpy array (or data structure base on numpy) ordered (latitude, pressure)
    nlevs = approximate number of contour levels to return (default 10)
    rmClev = if defined delete the contour level near this value
    sym = if defined make the contour intervals symmetric about zero
    verbose = if defined, print out some info to help debug
    """
    if nlevs is None: nlevs = 9

    dsub = data[np.isfinite(data)] # ignore NaN
    zmax = dsub.max()
    zmin = dsub.min()
    clevs = np.linspace(zmin, zmax, nlevs)
#    if nozero is None: nozero=0
#    if sym is None: sym=0

    if not verbose is None: print ("zmin, zmax", zmin, zmax)
    if zmax == zmin: 
        if (zmax != 0): 
          return np.array([0.9*zmax,zmax,1.1*zmax]);
        else:
          return np.array([-0.1,0.0, 0.1]);
    if not sym is None: 
        zcom = max(abs(zmin),abs(zmax));
        zmin = -zcom;
        zmax = +zcom;
        if not verbose is None: print ("zmin, zmax in sym", zmin, zmax)
    zrange = zmax/(zmin+1.e-20)
    okints = [0.1, 0.2, 0.5,1,2,5,10,20] # acceptable values for contour intervals
    okints = np.array(okints)
    zinc = (zmax-zmin)/nlevs;
    if not verbose is None: print ("zinc", zinc)
    pow = int(np.log10(zinc))
    if not verbose is None: print ("pow", pow)
    cints = okints*(10.**pow);
    if not verbose is None: print ("cints", cints)
    nlevsout = (zmax-zmin)/cints;
    if not verbose is None: print ("nlevsout",nlevsout)
#    np.where(V0 == V0.max())
    f1 = abs(nlevsout-nlevs)
    if not verbose is None: print ("f1", f1)
    f2 = np.where(f1 == f1.min());
    f2 = np.array(f2).flatten()
    if not verbose is None: print ("f2", f2)
    nlevbest = int(f2[0])
    if not verbose is None: print ("nlevbest, cintsbest",nlevbest, cints[nlevbest])
    if not verbose is None: print ("f3", zmin/cints[nlevbest])
    zminout = int(zmin/cints[nlevbest])*cints[nlevbest];
    zmaxout = int(zmax/cints[nlevbest])*(cints[nlevbest]);
    ninc = int((zmaxout-zminout)/cints[nlevbest] + 1.00000001);
    if not verbose is None: print ("ninc, zminout, zmaxout", ninc, zminout, zmaxout)
    clevs = np.arange(zminout, zmaxout*1.001, cints[nlevbest]);
    if not rmClev is None:    
        f4 = abs((clevs-rmClev)/max(abs(clevs)))
        if not verbose is None: print ("f4", f4)
        alist = np.argmin(f4);
        if not verbose is None: print ("alist", alist)
        if (np.size(alist) > 0): 
            zlist = np.where(clevs != clevs[alist])
            if not verbose is None: print ("zlist", zlist)
            if (np.size(zlist)> 0): 
                if not verbose is None: print ("list of nonzero contour levels", zlist)
                clevs = clevs[zlist];
    return clevs;

def make_fvarea(lon,lat,alt=None):
    """make_fvarea(lon,lat)

    make finite volume areas, given 1-D lon and lat arrays
    """
    pio180 = np.pi/180.
    yw = np.zeros(len(lat)+1)
    # walls in degrees
    yw[1:-1] = (lat[1:] + lat[:-1]) / 2
    if (np.abs(lat[-1]-lat[0])-180.) < 1.e-4:
        #print('polepoint')
        yw[0] = lat[0]
        yw[-1] = lat[-1]
    else:
        #print('nopole point')
        yw[0] = lat[0]-(lat[1]-lat[0])/2.
        yw[-1] = lat[-1]+(lat[1]-lat[0])/2.

    if alt == True:
        #print('alt True')
        dy = np.diff(yw*pio180)*np.cos(lat*pio180)
    else:
        # wall in sin(latitude in radians)
        yw = np.sin(yw*pio180)
        dy = np.diff(yw)
    dx = float(lon[1]-lon[0])*pio180
    dxa = np.full([len(lat),len(lon)],dx)
    area = dxa*dy[:,np.newaxis]

    return area

def xr_getvar(Varname, DS, regtag=None,long_name=None):
    """get variable Varname from xarray dataset DS
       some variables are derived from others
       some are modified to more traditional units
       some are returned directly from DS
       
       regtag is an optional "region tag" that is added to dimension names on EAM and CAM regional history file
    """
    if regtag == None:
        regtag = ""
        
    # see whether the variable is on DS
    if Varname in DS:
        on_DS = True
    else:
        on_DS = False

    # identify a variable that has time, and lev, and another spatial dimension for use in constructing 3D vars
    for vv in DS:
        dlist = list(DS[vv].dims)
        ndlist = len(dlist)
        if 'lev' in dlist and 'time' in dlist and ndlist > 2:
            #print('lev and time found and number of dims > 2',vv,dlist)
            nm3dv = vv
            break
        
    try:    # return derived variables, or modified variables, or variables on DS
        if Varname == "CLDLIQ":
            Var = DS['CLDLIQ'+regtag]
            Var = Var*1.e3
            Var.attrs['units'] = 'g/kg'
            Var.attrs["long_name"] = 'grid-avg liq. Mix. Rat.'
        elif Varname == "CLDICE":
            Var = DS['CLDICE'+regtag]
            Var = Var*1.e6
            Var.attrs['units'] = 'mg/kg'
            Var.attrs["long_name"] = 'grid-avg liq. Mix. Rat.'
        elif ((Varname == "CLDLOW") or (Varname == "CLDMED") or
              (Varname == "CLDHGH") or (Varname == "CLDTOT")):
            Var = DS[Varname+regtag]
            Var = Var*1.e2
            Var.attrs['units'] = '%'
        elif Varname == "CLOUD":
            Var = DS['CLOUD'+regtag]
            Var = Var*1.e2
            Var.attrs['units'] = '%'
        elif Varname == "Q":
            Var = DS['Q'+regtag]
            Var = Var*1.e3
            Var.attrs['units'] = 'g/kg'
        elif Varname == "ICWMR":
            Var = DS['ICWMR'+regtag]
            Var = Var*1.e3
            Var.attrs['units'] = 'g/kg'
            Var.attrs["long_name"] = 'in-cloud liq. Mix. Rat.'
        elif Varname == "ICWNC":
            Var = DS['ICWNC'+regtag]
            Var = Var*1.e-6
            Var.attrs['units'] = '/cm3'
        elif Varname == "P3": 
            # special treatment for constructing a 3D pressure from PS and
            # hybrid coefs
            #Var = (DS.hyam*DS.P0 + DS.hybm*DS['PS'+regtag])/100.
            Var = (DS['PS'+regtag]*DS.hybm + DS.hyam*DS.P0)/100.
            Var.attrs["units"] = 'hPa'
            Var.attrs["long_name"] = 'Pressure'
            # make sure the returned quantities have the same coordinate order as standard
            #ldims = list(DS['T'+regtag].dims)
            #Var = Var.transpose(*ldims)
            #print('newPin.dims', Pin.dims)
            #print('newPin.shape', Pin.shape)
        elif Varname == "P3i": 
            # special treatment for constructing a 3D pressure from PS and
            # hybrid coefs
            #Var = (DS.hyai*DS.P0 + DS.hybi*DS['PS'+regtag])/100.
            Var = (DS['PS'+regtag]*DS.hybi + DS.hyai*DS.P0)/100.
            Var.attrs["units"] = 'hPa'
            Var.attrs["long_name"] = 'Pressure(interfaces)'
            # make sure the returned quantities have the same coordinate order as standard
            #ldims = list(DS['T'+regtag].dims)
            #Var = Var.transpose(*ldims)
            #print('newPin.dims', Pin.dims)
            #print('newPin.shape', Pin.shape)
        elif Varname == "DPOG": 
            # special treatment for constructing a 3D pressure from PS and
            # hybrid coefs
            #VarI = (DS.hyai*DS.P0 + DS.hybi*DS['PS'+regtag])/100.
            VarI = (DS['PS'+regtag]*DS.hybi + DS.hyai*DS.P0)
            #print('VarI',VarI)
            #print('VarI col 0',VarI[0,0,:].values)
            #print('PS',DS['PS'+regtag])
            print('using this variable as a template for DPOG',nm3dv)
            Var = DS[nm3dv+regtag].copy()
            Var = Var.rename(Varname)
            #print('Var',Var)
            Varx = VarI.diff("ilev").values/9.8
            #print('Varx',Varx.shape)
            Var.data = Varx
            Var.attrs = {}
            #print('new Var col 0', Var[0,0,:].values)
            Var.attrs["units"] = 'kg/m2'
            #Var.attrs["basename"] = 'DPOG'
            Var.attrs["long_name"] = 'DeltaPressure(interfaces)_over_gravity'
#            latr = var.attrs
#            if 'standard_name' in latr.keys():
#                x = Var.attrs.pop("standard_name")
            #print('VarO',Var)
            # make sure the returned quantities have the same coordinate order as standard
            #ldims = list(DS['T'+regtag].dims)
            #Var = Var.transpose(*ldims)
            #print('newPin.dims', Pin.dims)
            #print('newPin.shape', Pin.shape)
        elif Varname == "ARELx": 
            V1 = DS['AREL']
            V2 = DS['FREQL']
            Var = V1/(V2+1.e-6)
            Var = Var.rename(Varname)
            Var.attrs["long_name"] = 'Estimated Droplet Effective Radius'
        elif Varname == "AREL":
            Var = DS[Varname]
            Var.attrs["long_name"] = 'Partial Droplet Effective Radius'
        elif Varname == "P3i": 
            # special treatment for constructing a 3D pressure from PS and
            # hybrid coefs
            #Var = (DS.hyai*DS.P0 + DS.hybi*DS['PS'+regtag])/100.
            Var = (DS['PS'+regtag]*DS.hybi + DS.hyai*DS.P0)/100.
            Var.attrs["units"] = 'hPa'
            Var.attrs["long_name"] = 'Pressure(interfaces)'
            # make sure the returned quantities have the same coordinate order as standard
            #ldims = list(DS['T'+regtag].dims)
            #Var = Var.transpose(*ldims)
            #print('newPin.dims', Pin.dims)
            #print('newPin.shape', Pin.shape)
        elif Varname == "DPOG": 
            # special treatment for constructing a 3D pressure from PS and
            # hybrid coefs
            #VarI = (DS.hyai*DS.P0 + DS.hybi*DS['PS'+regtag])/100.
            VarI = (DS['PS'+regtag]*DS.hybi + DS.hyai*DS.P0)
            required_dims = ('lon'+regtag, 'lat'+regtag, 'lev')
            mylist = find_vars_bydims(DS,required_dims)
            #print('VarI',VarI)
            #print('VarI col 0',VarI[0,0,:].values)
            #print('PS',DS['PS'+regtag])
            Var = DS[mylist[0]].copy()
            Var = Var.rename(Varname)
            #print('Var',Var)
            Varx = VarI.diff("ilev").values/9.8
            #print('Varx',Varx.shape)
            Var.data = Varx
            #print('new Var col 0', Var[0,0,:].values)
            Var.attrs["units"] = 'kg/m2'
            #Var.attrs["basename"] = 'DPOG'
            Var.attrs["long_name"] = 'DeltaPressure(interfaces)_over_gravity'
            #x = Var.attrs.pop("standard_name")
            #print('VarO',Var)
            # make sure the returned quantities have the same coordinate order as standard
            #ldims = list(DS['T'+regtag].dims)
            #Var = Var.transpose(*ldims)
            #print('newPin.dims', Pin.dims)
            #print('newPin.shape', Pin.shape)
        elif Varname == "CDNUMC":
            Var = DS[Varname+regtag]
            Var.attrs['long_name'] = 'Vert Int drop number'
        elif Varname == "PRECC":
            Var = DS['PRECC'+regtag]
            Var = Var*8.64e7
            Var.attrs['units'] = 'mm/day'
        elif Varname == "PRECT":
            if on_DS:
                Var = DS['PRECT'+regtag]
            else:
                Var = DS['PRECC'+regtag]+DS['PRECL'+regtag]
            Var = Var*8.64e7
            Var.attrs['long_name'] = 'Total(liq,ice,conv,strat) Precipitation'
            Var.attrs['units'] = 'mm/day'
        elif Varname == "area_unfinished":
            if on_DS:
                Var = DS['area'+regtag]
            else:
                lat = Var1['lat'].values
                lon = Var1['lon'].values
                area = make_fvarea(lon,lat)
                Var = DS['PRECC'+regtag]
            Var = Var*8.64e7
            Var.attrs['long_name'] = 'Total(liq,ice,conv,strat) Precipitation'
            Var.attrs['units'] = 'mm/day'
        elif Varname == "PRECL":
            if on_DS:
                Var = DS['PRECL'+regtag]
            else:
                Var = (DS['PRECT'+regtag]-DS['PRECC'+regtag]).rename(Varname)
            Var = Var*8.64e7
            Var.attrs['long_name'] = 'Stratiform (liq,ice) Precipitation'
            Var.attrs['units'] = 'mm/day'
        elif Varname == "RESTOM":
            Var = (DS['FSNT'+regtag]-DS['FLNT'+regtag]).rename(Varname)
            #Var.attrs['basename'] = Varname
            Var.attrs['long_name'] = 'Residual TOA flux'
            Var.attrs['units'] = 'W/m2'
        elif Varname == 'PS':
            Var = DS['PS'+regtag]
            Var = Var/100.
            Var.attrs['units'] = 'hPa'
        elif Varname == 'PCONVT':
            Var = DS['PCONVT'+regtag]
            Var = Var/100.
            Var.attrs['units'] = 'hPa'
        elif Varname == 'ZCONVT':
            Var = DS['PCONVT'+regtag].copy()
            Var = Var.rename(Varname)
            Var = Var/100.
            #Var.attrs['basename'] = Varname
            Var.attrs['units'] = 'm'
            Var = -8.1e3*np.log(Var/1012.)  # rough conversion to m using 8km scale height
            Var.attrs['long_name'] = 'convection top height'
            Var = Var.rename(Varname)
        elif Varname == 'LWCF':
            Var = DS['LWCF'+regtag]
            Var.attrs['long_name'] = 'LW Cloud Radiative Effect'
        elif Varname == 'SWCF':
            Var = DS['SWCF'+regtag]
            Var.attrs['long_name'] = 'SW Cloud Radiative Effect'
        elif Varname == "TGCLDLWP":
            Var = DS['TGCLDLWP'+regtag]
            Var = Var*1.e3
            Var.attrs['units'] = 'g/m2'
            Var.attrs["long_name"] = 'grid-avg LWP.'
        elif Varname == "AODVIS":
            Var = DS[Varname+regtag]
            Var.attrs['units'] = '1' #print('add units attribute')
        elif Varname == 'XXXXXX':
            print ('example of a complicated calculation')
            if on_DS:
                Var = DS[Varname]
                Var.attrs['long_name'] = 'Vert Int drop number'
                print('from file')
            else:
                VarI = (DS['PS'+regtag]*DS.hybi + DS.hyai*DS.P0) #interface pressures
                dpogdata = VarI.diff("ilev").values/9.8
                print('dpogdata shape', dpogdata.shape)
                print('NUMLIQ shape', DS['NUMLIQ'].shape)
                Varx = DS['NUMLIQ'].transpose(...,'lev')*dpogdata
                Var = DS['PS'].copy()
                Var.data = Varx.sum(dim='lev')
                Var = Var.rename(Varname)
                Var.attrs["units"] = '1/m2'
                Var.attrs["long_name"] = 'Vert. Int drop number'
                print ('derived')
        else:  # look in Xarray dataset DS
            if on_DS:
                Var = DS[Varname+regtag]
            else:
                estr = Varname+regtag+' is not defined or found in DS'
                raise UserWarning(estr)
    except KeyError:  # variable not a derived field or on DS
        estr = Varname+regtag+' is not defined or found in DS'
        raise UserWarning(estr)
    else:
        if not 'units' in Var.attrs:
            Var.attrs['units'] = "?"
        if not 'long_name' in Var.attrs:
            if 'standard_name' in Var.attrs:
                Var.attrs['long_name'] = Var.attrs['standard_name']
            else:
                Var.attrs['long_name'] = Varname
        if not long_name is  None:
            Var.attrs["long_name"] = long_name
        return Var

def xr_getvar_sl(VN, DS1, method='surface', verbose=False):
    """ get a field from a netCDF file.
    If it is a multi-level field, do something useful to provide single level information
    currently all it can do is return the surface value, or the max value of field below 850hPa
    
    This function assumes that the x-coordinate increases monotonically

    """
    Var1 = xr_getvar(VN,DS1)
    dimlist = Var1.dims
    if 'model_level_number' in dimlist:
        level_height = xr_getvar('level_height',DS1)
        sigma = xr_getvar('sigma',DS1)
        surface_altitude = xr_getvar('surface_altitude',DS1)
        altitude = level_height + (sigma * surface_altitude)
        altitude.attrs['long_name'] = 'altitude above mean sea-level'
        if method == 'surface':
            print('method:surface')
            V1 = Var1.isel(model_level_number=0)
            V1.attrs['long_name'] = V1.attrs['long_name'] + ' (surface level)'
        elif method == 'maxb850':
            dz1 = altitude - surface_altitude   # height above surface
            dz2 = (sigma - 1)*surface_altitude + level_height  # height above sea level
            # assume pressure has a 8.4km scale height, find the altitude of 850 hPa
            pmb = 850.
            psmb = 1000.
            scaleheight = 8.4e3
            altmb = -np.log(pmb/psmb)*scaleheight
            # find the max value of field below this altitude
            V1 = Var1.copy()
            V2 = V1.where(altitude <= altmb+50.)
            V3 = V2.max(dim='model_level_number')
            V3.attrs['long_name'] = 'max value below 850hPa of '+V2.name
            V1 = V3
    else:
        V1 = Var1
    
    if V1.attrs['units'] == 'mm/day':
          V1.attrs['units'] = 'mm/d'
          
    return V1


# a few useful code segments for manipulating xarray Datasets and DataArrays
# CODE from ESA CCI project at https://github.com/CCI-Tools/cate

from typing import Optional, Sequence, Union, Tuple

def get_lon_dim_name_impl(ds: Union[xr.Dataset, xr.DataArray]) -> Optional[str]:
    """
    Get the name of the longitude dimension.
    :param ds: An xarray Dataset
    :return: the name or None
    """
    return _get_dim_name(ds, ['lon', 'longitude', 'long'])


def get_lat_dim_name_impl(ds: Union[xr.Dataset, xr.DataArray]) -> Optional[str]:
    """
    Get the name of the latitude dimension.
    :param ds: An xarray Dataset
    :return: the name or None
    """
    return _get_dim_name(ds, ['lat', 'latitude'])


def _get_dim_name(ds: Union[xr.Dataset, xr.DataArray], possible_names: Sequence[str]) -> Optional[str]:
    for name in possible_names:
        if name in ds.dims:
            return name
    return None

def _normalize_lat_lon(ds: xr.Dataset) -> xr.Dataset:
    """
    Rename variables named 'longitude' or 'long' to 'lon', and 'latitude' to 'lon'.
    :param ds: some xarray dataset
    :return: a normalized xarray dataset, or the original one
    """
    lat_name = get_lat_dim_name_impl(ds)
    lon_name = get_lon_dim_name_impl(ds)

    name_dict = dict()
    if lat_name and 'lat' not in ds:
        name_dict[lat_name] = 'lat'

    if lon_name and 'lon' not in ds:
        name_dict[lon_name] = 'lon'

    if name_dict:
        ds = ds.rename(name_dict)

    return ds

def center_time(DS1):
    """center_time(DS1)
    
    correct the time coordinate in DS1 to represent the center of the time bounds
 
    """
    DS = DS1.copy()
    # the time coord is often registered at the end of the time averaging interval
    # reset to the midpoint of the interval
    time = DS['time'].copy()
    #print('time',time)
    #print('xxx',time.values)
    bndname = time.attrs['bounds']
    time_bnds = DS1[bndname]
    #print('time_bnds',time_bnds)
    tbdims = time_bnds.dims
    #print('tbdims',tbdims)
    tbd_name = ''
    # find the bounds name (the dim that isn't time)
    for tbd in tbdims:
        if tbd != 'time':
            tbd_name = tbd
    #print('tbd_name',tbd_name)
    #print('tbdims',tbdims)
    # if no bounds, then do nothing
    if tbd_name == '':
        return DS
    else:
        #tb = time_bnds.values
        #print('time_bnds',time_bnds)
        tbm = time_bnds.mean(dim=tbd_name).values
        #print('yyy',tbm)
        DS.coords["time"] = tbm
        DS['time'].attrs['long_name'] = 'time'
        DS['time'].attrs['bounds'] = 'time_bnds'
        return DS

def reconcile_xr_coords(V1G, V2G, rtol=1.e-8,verbose=False):
    """
    usage: V1, V2 = reconcile_xr_coords(V1, V2, rtol=) 
    check whether the Coordinates in two xarray objects (DataSet or DataArray) are close. 
    if they are identical do nothing 
    if they differ below rtol, then set V2 coordinate to the V1 values 
    if they exceed rtol then issue an error message and stop 
    
    The code is still pretty rough. 
    It could be improved by 
    allowing an exclude list, other than the time coord which is already avoided 
    checking whether coords are the same in V1 and V2 

    """
    #make copies, in case we dont want to overwrite original coords
    V1 = V1G.copy()
    V2 = V2G.copy()
    coords = V1.coords
    coord2 = V2.coords
    for c in coords:
        if verbose:
            print('comparing', c) 
        if c == 'time':
            if verbose:
                print('skipping',c) 
            continue
        if c not in coord2:
            print(c,' absent in second object') 
            continue
        V1cv = V1[c].values
        V2cv = V2[c].values
        relerr = 2.*np.abs(V1cv-V2cv)/(np.abs(V1cv+V2cv)+1.e-16)
        relerrm = np.max(relerr)
        #print('relerrm',relerrm) #xr.testing.assert_allclose(V1[c],V2[c])
        if (relerrm == 0.):
            continue
        if (relerrm < rtol):
            print('over-write coord ', c)
            V2.coords[c] = V1.coords[c]
            continue
        print('large error', relerrm,' in coord ',c)
        print('V1 coord',V1[c].values)
        print('V2 coord',V2[c].values)
        1./0
    return V1, V2;

def tavg_mon_wt(xr_var):
    """
    time avg,, weight by days in each month
    see https://ncar.github.io/esds/posts/2021/yearly-averages-xarray/
    """
    # Determine the month length
    month_length = xr_var['time'].dt.days_in_month
    # Calculate the weights
    wgts = month_length.groupby("time.year") / month_length.groupby("time.year").sum()
    # Make sure the weights in each year add up to 1
    np.testing.assert_allclose(wgts.groupby("time.year").sum(xr.ALL_DIMS), 1.0)
    
    # Subset our dataset for our variable
    obs = xr_var

    # Setup our masking for nan values
    cond = obs.isnull()
    ones = xr.where(cond, 0.0, 1.0)

    # Calculate the numerator
    obs_sum = (obs * wgts).resample(time="YS").sum(dim="time")

    # Calculate the denominator
    ones_out = (ones * wgts).resample(time="YS").sum(dim="time")
    wavg = obs_sum / ones_out
    wavg['time'] = wavg.time+timedelta(days=182.5)
    if len(xr_var['time'].values) > 1:
        if len(wavg['time'].values)*12 != len(xr_var['time'].values):
            raise Exception("time series inconsistency, either variable input partial year provided or input times are not correct")
    # Return the weighted average
    return wavg

def xr_cshplot(xrVar, xrLon, xrLat, plotproj=None, ax=None, cax=None,ylabels=None,clevs=None, cmap=None, title=None,cbar='default'):
    """xr_cshplot xarray cubed sphere horizontal plot
    """

    dinc = 1.  # increment of mesh in degrees
    lon_h=np.arange(np.floor(xrLon.min().values),np.ceil(xrLon.max().values+dinc), dinc)
    lat_h=np.arange(np.floor(xrLat.min().values),np.ceil(xrLat.max().values+dinc), dinc)
    xv,yv=np.meshgrid(lon_h,lat_h)
    data_regridded = interp_ap(xv, yv, xrVar.values,xrLat.values,xrLon.values)
    df = data_regridded.flatten()
    dsub = df[np.isfinite(df)] # ignore NaN
    zmax = dsub.max()
    zmin = dsub.min()
    #print('masked interpolated range',zmin,zmax)
    dataproj=ccrs.PlateCarree()    # data is always assumed to be lat/lon
    #plotproj=ccrs.Orthographic(central_latitude=0,central_longitude=55)   # any projections should work
    #print('plotproj is ',plotproj)
    if plotproj is None: plotproj = ccrs.Mercator
    if ax is None: ax = plt.gca()
    if cax is None: cax = ax
    if ylabels is None: ylabels = True
    if clevs is None:
        clevs = findNiceContours(np.array([zmin,zmax]),nlevs=10)
    #print('clevs',clevs)
    if cmap is None:
        #print('aaa, grabbing cmap default')
        #cmap = mpl.cm.get_cmap()
        cmap = plt.get_cmap()
        #print('bbb',cmap.N)
    #print('cmap',cmap)
    extend = 'both'
    norm = mpl.colors.BoundaryNorm(clevs,cmap.N,extend=extend)
    #print('norm',norm(clevs))

    plotproj=ccrs.Mollweide(central_longitude=200)   # any projections should work 
    #clat = (lat.values.min()+lat.values.max())/2.
    #clon = (lon.values.min()+lon.values.max())/2.
    #plotproj=ccrs.NearsidePerspective(central_longitude=clon, central_latitude=clat)
    plotproj=ccrs.Mercator()
    #ax = plt.axes(projection=plotproj)
    #ax.set_extent([lon.values.min(), 260., lat.values.min(), lat.values.max()])
    #ax.set_global()
    pl = ax.contourf(xv, yv, data_regridded, levels=clevs, # vmin=zmin, vmax=zmax,
                     norm=norm, cmap=cmap,
                     extend=extend, transform=ccrs.PlateCarree())
    if cbar == 'default':
        # Add colorbar to plot
        cb = plt.colorbar(
            pl, orientation='horizontal',ticks=clevs,ax=ax,
            label='%s (%s)'%(xrVar.long_name, xrVar.units), pad=0.1
        )
        cb.ax.tick_params(labelsize=8)
        
    if not title is None:
        ax.set_title(title)
        
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,
                      linewidth=2, color='gray', alpha=0.5)
    gl.left_labels=ylabels
    gl.right_labels=ylabels
    ax.coastlines(linewidth=1,color='blue')
    return pl

def interp_ap(xt, yt, data2d,lat,lon,method=None):
    """
    # interp an arbitrary set of points at xt, yt 
    # from data on an unstructured mesh at lat, lon
    #
    # interpolating in lat/lon space has issues with triangulation 
    # at pole and wrapping at greenwich, so interpolate in stereographic projection:
    #
    # input:
    #    data2d(ncol,...),lat(ncol),lon(ncol): data and coords on unstructured mesh
    #    data2d can be multidimensional array, but ncol must be first coordinate
    #    xt, yt: lat and lon coordinates of locations to interpolate to
    #    method: optional, use cubic interpolation if method='cubic'
    #
    # output 
    #    returns an array with same shape as xt with interpolated data
    #
    """
    from scipy.interpolate import LinearNDInterpolator
    from scipy.interpolate import CloughTocher2DInterpolator
    
    intp2D = LinearNDInterpolator
    if method == 'cubic':
        intp2D = CloughTocher2DInterpolator

    ld = data2d.shape[0] # length of first coord of input data
    lx = lon.shape[0]
    ly = lat.shape[0]
    if ((ld != lx) | (ld != ly)):
        print('inconsistent data2d, lon, lat arrays', ld, lx, ly)
        raise TypeError("inconsistent input in interp_ap")
    
    # mesh grid
    dproj=ccrs.PlateCarree()

    # select interpolation points located in the nh and sh
    inds = np.where(yt <= 0)
    indn = np.where(yt > 0)

    xtn = xt[indn]
    ytn = yt[indn]
    xts = xt[inds]
    yts = yt[inds]

    # take source data in the correct hemisphere, include extra halo points for interpolation
    # using the full global data sometimes confuses interpolation with points being mapped close to infinity
    halo = 15 # degrees
    data2d_h=data2d[lat<halo]

    lon_h=lon[lat<halo]
    lat_h=lat[lat<halo]
    coords_in  = ccrs.SouthPolarStereo().transform_points(dproj,lon_h,lat_h)

    #data_i = np.empty_like(xt,dtype=data2d.dtype)
    #data_i[:] = np.nan
    dims = list(xt.shape)+list(data2d[0,...].shape)
    #print('dims',dims)
    #data_i = np.empty_like(xt,dtype=data2d.dtype)
    data_i = np.zeros(dims,dtype=data2d.dtype)
    data_i[:] = np.nan
    #print ('bbb',data_i.shape)


    data_s = []
    if len(yts) > 0:
        cto = ccrs.SouthPolarStereo().transform_points(dproj,xts,yts)
        interp = intp2D(coords_in[:,0:2], data2d_h)
        data_s = interp(cto[:,0],cto[:,1])
        data_i[inds] = data_s

    data2d_h=data2d[lat>-halo]
    lon_h=lon[lat>-halo]
    lat_h=lat[lat>-halo]
    coords_in  = ccrs.NorthPolarStereo().transform_points(dproj,lon_h,lat_h)

    data_n = []
    if len(ytn) > 0:
        cto = ccrs.NorthPolarStereo().transform_points(dproj,xtn,ytn)
        interp = intp2D(coords_in[:,0:2], data2d_h)
        data_n = interp(cto[:,0],cto[:,1])
        data_i[indn] = data_n

    return data_i
