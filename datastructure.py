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
    

class Polygon:

    def __init__(self):
        self.nodes = Nodes()
        self.edges = Edges()
        self.convex_nodes = Nodes()
        self.convex_edges = Edges()
    
    def from_polygon(self, geom):
        lon = numpy.array(geom.exterior.xy[0])
        lat = numpy.array(geom.exterior.xy[1])
        self.nodes.from_lonlat(lon, lat)
        self.edges.from_lonlat(lon, lat)
        
    def get_convex(self):
        x, y, z = self.nodes.as_ecef().transpose()
        convex_node_indices = sphereGIS.xyz2convex(x,y,z)
        convex_lon = self.nodes.lon[convex_node_indices]
        convex_lat = self.nodes.lat[convex_node_indices]
        self.convex_nodes.from_lonlat(convex_lon, convex_lat)
        self.convex_edges.from_lonlat(convex_lon, convex_lat)
        
        
class MultiPolygon:
    
    def from_multipolygon(self, geom):  
        for p in list(geom):
            
            self.nodes.lon += p.exterior.xy[0]
            self.nodes.lat += p.exterior.xy[1]


class Edges:
    
    def __init__(self):
        self.lat = []
        self.lon = []
        self.gcs = []
        self.left_gcs = []
        self.right_gcs = []
        
    def from_lonlat(self, lon, lat):
        self.lon = lon 
        self.lat = lat
        from_nodes = to_ecef(lon, lat)
        to_nodes = numpy.roll(from_nodes, -1, axis=0) # Moving the first element to last
        self.gcs = numpy.cross(from_nodes, to_nodes)
        self.left_gcs = numpy.cross(self.gcs, from_nodes)
        self.right_gcs = numpy.cross(to_nodes, self.gcs)
        
    def as_df(self):
        edges = shapely.geometry.LineString(zip(self.lon, self.lat))
        gdf = geopandas.GeoDataFrame({'geom': [edges]})
        gdf = gdf.set_geometry('geom')
        return gdf


        
class Nodes:
    
    def __init__(self):
        self.lon = []
        self.lat = []
        
    def from_lonlat(self, lon, lat):
        self.lon = lon 
        self.lat = lat 
        
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
        
    def as_df(self):
        points = []
        for lon, lat in zip(self.lon, self.lat):
            points.append(shapely.geometry.Point(lon, lat))
        gdf = geopandas.GeoDataFrame({'geom': points})
        gdf = gdf.set_geometry('geom')
        return gdf
        
    


class Storage:
    def from_geom(self, geom):
        geom_type = geom.type
        if  geom_type == 'Polygon':
            self.from_polygon(geom)
        elif geom_type == 'MultiPolygon':
            self.from_multipolygon(geom)
        




class Mod09:
    def __init__(self, fname):
        #fname = '../data/MOD09.A2020032.1940.006.2020034015024.hdf'
        hdf = pyhdf.SD.SD(fname)
        self.lat = hdf.select('Latitude').get().flatten()
        self.lon = hdf.select('Longitude').get().flatten()
        self.points = to_ecef(self.lon, self.lat)
        
    def intersects_convex(self, convex_edges):
        dots = numpy.einsum('ji,mi->jm', self.points, convex_edges.gcs)
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

