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
    (double* z, int len_z)
}


%apply (int * INPLACE_ARRAY1, int DIM1 ) {
    (int* convex_nodes, int len_nodes)
}


%pythoncode %{
import numpy

def xyz2convex(x, y, z):
    out = numpy.full(x.shape, [-1], dtype=numpy.int32)
    _find_convex_hull(x, y, z, out)
    out = out[out!=-1]
    return out
%}

%include "sphereGIS.h"
