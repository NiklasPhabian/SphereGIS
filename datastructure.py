import sphereGIS
import random
import numpy
import math
import shapely
import geopandas
import pandas
import pyhdf.SD


def to_ecef(lon, lat):
    lon = numpy.array(lon)
    lat = numpy.array(lat)
    x = numpy.cos(lon/360*math.pi*2) * numpy.cos(lat/360*math.pi*2)
    y = numpy.sin(lon/360*math.pi*2) * numpy.cos(lat/360*math.pi*2)
    z = numpy.sin(lat/360*math.pi*2)
    return numpy.array([x, y, z]).transpose()



class SphereGeoDataFrame(geopandas.GeoDataFrame):
    
    def __init__(self, *args, **kwargs):
        super(SphereGeoDataFrame, self).__init__(*args, **kwargs)
        self.bootstrap()
        
    def bootstrap(self):
        spherical_polygons = []
        convex_hulls = []
        for index, row in self.iterrows():
            spherical_polygon = SphericalPolygon(row.geometry)
            spherical_polygons.append(spherical_polygon)
            convex_hulls.append(spherical_polygon.convex_nodes.as_polygon())
        self['sphere_geom'] = spherical_polygons
        self['convex_hull'] = convex_hulls
    

class SphericalPolygon:

    def __init__(self, geom):
        self.nodes = Nodes()
        self.edges = Edges()
        self.convex_nodes = Nodes()
        self.convex_edges = Edges()
        if geom:
            self.from_geom(geom)
            #self.get_convex()
    
    def from_geom(self, geom):
        geom_type = geom.type
        if  geom_type == 'Polygon':
            self.from_polygon(geom)
        elif geom_type == 'MultiPolygon':
            self.from_multipolygon(geom)
        
    def from_polygon(self, geom):
        self.nodes.from_lonlat(geom.exterior.xy[0], geom.exterior.xy[1])
        self.edges.from_lonlat(geom.exterior.xy[0], geom.exterior.xy[1])
    
    def from_multipolygon(self, geom):  
        for p in list(geom):          
            lon = p.exterior.xy[0]
            lat = p.exterior.xy[1]
            new_nodes = Nodes()
            new_nodes.from_lonlat(lon, lat)
            self.nodes.add(new_nodes)
            new_edges = Edges()
            new_edges.from_lonlat(lon, lat)
            self.edges.add(new_edges)
        
    def get_convex(self):
        x, y, z = self.nodes.as_ecef().transpose()
        convex_node_indices = sphereGIS.xyz2convex(x,y,z)
        convex_lon = self.nodes.lon[convex_node_indices]
        convex_lat = self.nodes.lat[convex_node_indices]
        self.convex_nodes.from_lonlat(convex_lon, convex_lat)
        self.convex_edges.from_lonlat(convex_lon, convex_lat)
        
    def get_convex_indices(self):
        x, y, z = self.nodes.as_ecef().transpose()
        convex_node_indices = sphereGIS.xyz2convex(x,y,z)
        return convex_node_indices 
        
        
class Edges:
    
    def __init__(self):
        self.gcs = numpy.empty([0,3], dtype=numpy.float)
        self.left_gcs = numpy.empty([0,3], dtype=numpy.float)
        self.right_gcs = numpy.empty([0,3], dtype=numpy.float)
        
    def from_lonlat(self, lon, lat):
        from_nodes = to_ecef(lon, lat)
        to_nodes = numpy.roll(from_nodes, -1, axis=0) # Moving the first element to last
        self.gcs = numpy.cross(from_nodes, to_nodes)
        self.left_gcs = numpy.cross(self.gcs, from_nodes)
        self.right_gcs = numpy.cross(to_nodes, self.gcs)
        
    def add(self, other):
        self.gcs = numpy.concatenate((self.gcs, other.gcs))
        self.left_gcs = numpy.concatenate((self.left_gcs, other.left_gcs))
        self.right_gcs = numpy.concatenate((self.right_gcs, other.right_gcs))

        
class Nodes:
    
    def __init__(self):
        self.lon = numpy.empty([0], dtype=numpy.float)
        self.lat = numpy.empty([0], dtype=numpy.float)
        
    def from_lonlat(self, lon, lat):
        self.lon = numpy.array(lon)
        self.lat = numpy.array(lat)
        
    def drop_repetitive_unsorted(self):
        # Dirty trick to get rid of duplicate points
        self.lon, self.lat = zip(*(set(zip(self.lon, self.lat))))
                
    def drop_repetitive(self):
        # Even dirtier; this keeps the nodes in order.
        self.lat, self.lon = zip(*list(dict.fromkeys(zip(self.lat, self.lon))))
        
    def shuffle(self):
        a = list(zip(self.lat, self.lon))
        random.shuffle(a)
        self.lat, self.lon = zip(*a)
        
    def as_ecef(self):
        return to_ecef(self.lon, self.lat)
    
    def to_csv(self, name):
        x, y, z = self.as_ecef().transpose()
        df = pandas.DataFrame({'x':x, 'y': y, 'z': z})
        df.to_csv(name, index=None, header=None)
        
    def as_polygon(self):
        return shapely.geometry.Polygon(zip(self.lon, self.lat))
    
    def as_polygon_df(self):
        geom = self.as_polygon()        
        gdf = geopandas.GeoDataFrame({'geom': [geom]})
        gdf = gdf.set_geometry('geom')
        return gdf
        
    def as_point_df(self):
        points = []
        for lon, lat in zip(self.lon, self.lat):
            points.append(shapely.geometry.Point(lon, lat))
        gdf = geopandas.GeoDataFrame({'geom': points})
        gdf = gdf.set_geometry('geom')
        return gdf
    
    def as_line_df(self):
        # This will give surprising results for MultiPolygon
        edges = shapely.geometry.LineString(zip(self.lon, self.lat))
        gdf = geopandas.GeoDataFrame({'geom': [edges]})
        gdf = gdf.set_geometry('geom')
        return gdf
    
    def add(self, other):
        self.lon = numpy.concatenate((self.lon, other.lon))
        self.lat = numpy.concatenate((self.lat, other.lat))


class Granule:
    
    def __init__(self):
        self.lat = None
        self.lon = None
        self.ecef = None
        
    def intersects_convex(self, convex_edges):
        dots = numpy.einsum('ji,mi->jm', self.ecef, convex_edges.gcs)
        n_constraints = len(convex_edges.gcs)
        inside_convex = numpy.where(numpy.sum(dots>0, axis=1)==n_constraints)[0]
        return inside_convex
            
    def inside_polygon(self, polygon):
        inside_convex = self.intersects_convex(polygon.convex_edges)
        lon_nodes = polygon.nodes.lon
        lat_nodes = polygon.nodes.lat     
        candidate_lon = self.lon[inside_convex]
        candidate_lat = self.lat[inside_convex]
        mask = sphereGIS.intersects(candidate_lon, candidate_lat, lon_nodes, lat_nodes)
        return inside_convex[numpy.where(mask)[0]]
        

class Mod09(Granule):
    
    def __init__(self, fname):
        hdf = pyhdf.SD.SD(fname)
        self.lat = hdf.select('Latitude').get().flatten()
        self.lon = hdf.select('Longitude').get().flatten()
        self.ecef = to_ecef(self.lon, self.lat)
        


