void _find_convex_hull( double* lon, int len_lon,
                        double* lat, int len_lat,                          
                        int* convex_nodes, int len_nodes);


void _intersects_convex(double* lon, int len_lon, 
                 double* lat, int len_lat,                  
                 double* x, int len_x, 
                 double* y, int len_y, 
                 double* z, int len_z, 
                 int* intersects, int len_intersects);



void _intersects(double* lon_points, int len_lon_points, 
                 double* lat_points, int len_lat_points,                  
                 double* lon_nodes, int len_lon_nodes, 
                 double* lat_nodes, int len_lat_nodes, 
                 int* intersects, int len_intersects);
