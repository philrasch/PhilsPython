import cdms2

#f = cdms2.open("/Users/d3x345/Desktop/NetCDF_files/CERES_EBAF-Surface_Ed2.8_Subset_200003-201409.nc");
f = cdms2.open("/Users/d3x345/Desktop/NetCDF_files/COREV2_ncar_rad.15JUNE2009.nc");
atts = f.attributes.keys;
print "global attributes "
for a in f.attributes.keys():
    x = getattr(f,a);
    print "a", a;
    print "x", x;
#   print '\t:'+a+'\t= "'+repr(getattr(f,a))+'";'
    

# Get the names of the dimensions present in files
dims = f.listdimension();

print dims;
for dim in dims:
    print "dim is ", dim;
    D = f[dim];
#    print "dir(dim)", dir(dim), "x";
#    print "dir(D)", dir(D), "x";
    print "D.attributes", D.attributes;
    datts = D.attributes.keys();
    dattv = D.attributes.values();
#    print "dir(D.attributes)", dir(D.attributes);
    print "datts", datts;
    print "dattv", dattv;
#    dimtype = dattv[3].dtype;
#    print "x", dimtype;
#    for da,dv in datts, dattv:
#        print "da,dv ", da, dv;
    data = D[:];
    t = data.dtype;
    print "type", t;
#   print "data ", data;
    print "range: ", data.min(), ":",data.max();


# Get attributes attached to a variable (here clt)
#clt_att = f.listattribute('clt')
#print clt_att
# ['units', 'time_statistic', 'long_name', 'grid_name', 'comments', 'missing_value', 'grid_type']
varlist = f.listvariables();

print varlist;
for var in varlist:
    V = f[var][0];
    print "var", var, V.info();
    
quit();
