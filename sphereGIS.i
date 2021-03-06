%module sphereGIS

%{
  #define SWIG_FILE_WITH_INIT  /* To import_array() below */
  #include "sphereGIS.h"
%}

%include "numpy.i"


%init %{
    import_array();
%}


/* Applying the typemaps */

%apply (double * IN_ARRAY1, int DIM1) {
    (double* x, int len_x),
    (double* y, int len_y),
    (double* z, int len_z),
    (double* gc_x, int len_x),
    (double* gc_y, int len_y),
    (double* gc_z, int len_z),
    (double* lon, int len_lon),
    (double* lat, int len_lat),
    (double* lon_points, int len_lon_points),
    (double* lat_points, int len_lat_points),                  
    (double* lon_nodes, int len_lon_nodes), 
    (double* lat_nodes, int len_lat_nodes)
}



%apply (int * INPLACE_ARRAY1, int DIM1 ) {
    (int* convex_nodes, int len_nodes),
    (int* intersects, int len_intersects) 
}


%pythoncode %{
import numpy

def find_first_indices(lon_nodes, lat_nodes):
    out = numpy.full(2, [-1], dtype=numpy.int32)
    _find_first_indices(lon_nodes, lat_nodes, out)    
    return out[0], out[1]
    
def lonlat2convex(lon_nodes, lat_nodes):
    out = numpy.full(lon_nodes.shape, [-1], dtype=numpy.int32)
    _find_convex_hull(lon_nodes, lat_nodes, out)
    out = out[out!=-1]
    return out
    
def intersects(lon_points, lat_points, lon_nodes, lat_nodes):
    out = numpy.full(lon_points.shape, [0], dtype=numpy.int32)
    _intersects(lon_points, lat_points, lon_nodes, lat_nodes, out)
    out = numpy.array(out, dtype=numpy.bool)
    return out
    
def intersects_convex_gc(lon_points, lat_points, gc_x, gc_y, gc_z):
    out = numpy.full(lon_points.shape, [1], dtype=numpy.int32)
    _intersects_convex(lon_points, lat_points, gc_x, gc_y, gc_z, out)
    out = numpy.array(out, dtype=numpy.bool)
    return out
    


%}

%include "sphereGIS.h"
