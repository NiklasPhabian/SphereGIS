import sys
sys.path.insert(0, '../') 
import sphereGIS 
import datastructure
import geopandas


# 3:4 is trinidad
polygons = geopandas.read_file(geopandas.datasets.get_path('naturalearth_lowres'))#[4:5]
polygons = polygons[polygons.continent=='North America']


polygons = geopandas.read_file('../data/caribbean.gpkg')[20:21]
geom = polygons.iloc[0].geometry
polygon = datastructure.SphericalPolygon(geom)
polygon.get_convex()
