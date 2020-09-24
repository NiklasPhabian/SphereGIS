#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <math.h>       
#include <set>

struct ECEF { 
    double x, y, z;
    ECEF(double x, double y, double z): x(x), y(y), z(z) {}
    ECEF() {}
};

ECEF cross(ECEF & v1, ECEF& v2) {
    double x = v1.y * v2.z - v1.z * v2.y;
    double y = v1.z * v2.x - v1.x * v2.z;
    double z = v1.x * v2.y - v1.y * v2.x;
    return ECEF(x, y, z);
}

double dot(ECEF & v1, ECEF & v2) {
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

bool equals(ECEF & v1, ECEF & v2) {
    return (v1.x==v2.x) && (v1.y==v2.y) && (v1.z==v2.z);
    //return fabs(v1.x-v2.x)<e && fabs(v1.y-v2.y)<e && fabs(v1.z-v2.z)<e;
        
}

typedef std::vector<ECEF> ECEFs;

ECEFs make_ecefs(double* x, double* y, double* z, int len) {
    ECEFs ecef;        
    for (int i = 0; i < len; i++) {
        ecef.push_back(ECEF(x[i], y[i], z[i]));   
    }
    return ecef;
}


ECEFs make_ecefs(double* lon, double* lat, int len) {
    ECEFs ecef;
    for (int i = 0; i < len; i++) {      
        double x = cos(lon[i]/360*M_PI*2) * cos(lat[i]/360*M_PI*2);
        double y = sin(lon[i]/360*M_PI*2) * cos(lat[i]/360*M_PI*2);
        double z = sin(lat[i]/360*M_PI*2);
        ecef.push_back(ECEF(x, y, z));
    }
    return ecef;
}


class Vertex;
typedef std::vector<Vertex> Vertices;

class Vertex: public ECEF {
    public:        
        int index;
        
        Vertex() {}
        
        Vertex(double x_, double y_, double z_) {
            x = x_;
            y = y_; 
            z = z_;
        }
        
        Vertex(double x_, double y_, double z_, int index_) {
            x = x_;
            y = y_; 
            z = z_;
            index = index_;
        }
    
        Vertex(Vertices &vertices, int index_) {
            x = vertices[index_].x;
            y = vertices[index_].y;
            z = vertices[index_].z;
            index = index_;
        }
};



Vertices read_ecef_csv(std::string file_path) {
    Vertices vertices;
    std::ifstream file(file_path);
    std::string row; 
    int index = 0;
    while (std::getline(file, row)) {
        int first_comma = row.find(",");
        int second_comma = row.find(",", row.find(",")+1);
        double x = stof(row.substr(0, first_comma));        
        double y = stof(row.substr(first_comma+1));
        double z = stof(row.substr(second_comma+1));
        vertices.push_back(Vertex(x, y, z, index)); 
        index ++;
    }
    file.close();    
    return vertices;
}


Vertices make_vertices(double* x, double* y, double* z, int len) {
    Vertices vertices;        
    for (int i = 0; i < len; i++) {        
        vertices.push_back(Vertex(x[i], y[i], z[i], i));   
    }
    return vertices;
}


Vertices make_vertices(double* lon, double* lat, int len) {
    Vertices vertices;
    for (int i = 0; i < len; i++) {     
        double x = cos(lon[i]/360*M_PI*2) * cos(lat[i]/360*M_PI*2);
        double y = sin(lon[i]/360*M_PI*2) * cos(lat[i]/360*M_PI*2);
        double z = sin(lat[i]/360*M_PI*2);
        vertices.push_back(Vertex(x, y, z, i));        
    }
    return vertices;
}


class ConvexEdge {
    public:        
        Vertex from_vertex;
        Vertex to_vertex;
        ECEF gc;
        
        ConvexEdge(Vertex from_vertex_, Vertex to_vertex_) {            
            from_vertex = from_vertex_;
            to_vertex = to_vertex_;
            gc = cross(from_vertex_, to_vertex_); 
        }
        
        bool point_outside(Vertex vertex) {            
            if (dot(gc, vertex) <= 1e-18) {                
                return true;
            } else {
                return false;
            }
        }
};

class ConvexEdges {
    public:
        std::vector<ConvexEdge> edges;
                        
        int circular_index(int index) {
            index = index % edges.size();               
            if (index <0) index *=-1 ;
            return index;
        }
    
        void push_back(ConvexEdge edge){
            edges.push_back(edge);
        }
        
        void erase(int index) {            
            index = circular_index(index);
            edges.erase(edges.begin() + index);           
        }
        
        void insert(int index, ConvexEdge edge) {
            if (index==edges.size()) {
                edges.push_back(edge);
            } else {
                index = circular_index(index);
                edges.insert(edges.begin() + index, edge);            
            }
        }
        
        ConvexEdge &get(int index) {            
            index = circular_index(index);              
            return edges[index];
        }
        
        int size() {
            return edges.size();
        }
};



class ConvexHull {
    public:
        ConvexEdges edges;
        
        ConvexHull () {}
        
        ConvexHull (Vertex first_node, Vertex second_node) {
            edges.push_back(ConvexEdge(first_node, second_node));
            edges.push_back(ConvexEdge(second_node, first_node));
        }
        
        void insert_vertex(Vertex vertex) {            
            if (already_in_hull(vertex)) return;
            
            int insertion_edge_index = find_insertion_edge_index(vertex);            
            if (insertion_edge_index == -1) return;           
            
            pop_edges(insertion_edge_index, vertex);
            span_new_edge(insertion_edge_index, vertex);           
        }
        
        int find_insertion_edge_index(Vertex vertex) {
            for (int i=0; i<= edges.size(); i++) {
                if (edges.get(i).point_outside(vertex)) {                    
                    return i;
                }
            }            
            return -1;
        }
        
        bool already_in_hull(Vertex candidate) {         
            for (int i = 0; i<= edges.size(); i++){                                
                if (equals(edges.get(i).from_vertex, candidate)) {
                    std::cout << candidate.index << " " <<  edges.get(i).from_vertex.index << std::endl;
                    return true;                
                }
            }
            return false;
        }
        
        void pop_edges(int insertion_edge_index, Vertex vertex) {               
            edges.erase(insertion_edge_index);            
            while (edges.get(insertion_edge_index).point_outside(vertex)) {
                edges.erase(insertion_edge_index);
            }
        }
                
        void span_new_edge(int insertion_edge_index, Vertex vertex) {
            Vertex previous_vertex = edges.get(insertion_edge_index-1).to_vertex;
            Vertex next_vertex = edges.get(insertion_edge_index).from_vertex;
            
            ConvexEdge new_edge_1 = ConvexEdge(previous_vertex, vertex);
            ConvexEdge new_edge_2 = ConvexEdge(vertex, next_vertex);
            edges.insert(insertion_edge_index, new_edge_2);
            edges.insert(insertion_edge_index, new_edge_1);
        }
};




//////////////////////////
// To find first edge  //
/////////////////////////


bool test_edge(int from_index, int to_index, Vertices &vertices) {
    if (equals(vertices[from_index], vertices[to_index])) return false;        
    ECEF great_circle = cross(vertices[from_index], vertices[to_index]);
    
    for (int node_index = 0; node_index< vertices.size(); node_index++){    
        if (node_index == from_index) continue;
        if (node_index == to_index) continue;
        if (equals(vertices[node_index], vertices[to_index])) continue;
        if (equals(vertices[node_index], vertices[from_index])) continue;
        if (dot(great_circle, vertices[node_index]) < 0) {
            return false;
        }
    }   
    return true;
}

int find_to_index(int from_index, Vertices &vertices) {
    int to_index = from_index;
    int sign = 1;
    sign *= -1;
    int n_vertices=  vertices.size();
    for (int i = 1; i < n_vertices; i++) {      
        // We move from the from_index left and right in index space
        sign *= -1;        
        to_index = (to_index + i * sign) % n_vertices;
        if (to_index < 0) to_index = n_vertices + to_index;
        if (test_edge(from_index, to_index, vertices)) return to_index;
    }
    return -1;
}

ConvexEdge find_first_edge(Vertices &vertices) {
    int to_index;
    int from_index;
    
    // We first try the polygons edges in both directions
    for (from_index = 1; from_index < vertices.size()-1; from_index++) {        
        to_index = from_index + 1;
        if (test_edge(from_index, to_index, vertices)){
            return ConvexEdge(vertices[from_index], vertices[to_index]);
        }
//         
        to_index = from_index - 1;
        if (test_edge(from_index, to_index, vertices)) {
            return ConvexEdge(vertices[from_index], vertices[to_index]);
        }
    }
        
    // If that fails, we try to draw edges
    for (from_index = 0; from_index < vertices.size(); from_index++){
        to_index = find_to_index(from_index, vertices);        
        if ( to_index >= 0) {
            return ConvexEdge(vertices[from_index], vertices[to_index]);
        } 
    }       
    
    return ConvexEdge(vertices[0], vertices[0]);
}


ConvexHull find_convex_hull(Vertices &vertices) {        
    ConvexEdge first_edge = find_first_edge(vertices);       
    ConvexHull convex_hull = ConvexHull(first_edge.from_vertex, first_edge.to_vertex);    
    for (auto const& vertex: vertices) {    
        //std::cout << vertex.index << std::endl;
        convex_hull.insert_vertex(vertex);        
    }        
    
    return convex_hull;
}


////////////////////////////////
// Exposed functions
////////////////////////////////

void _find_first_indices(double* lon, int len_lon,
                        double* lat, int len_lat,                          
                        int* convex_nodes, int len_nodes) {
    Vertices vertices = make_vertices(lon, lat, len_lat);
    ConvexEdge first_edge = find_first_edge(vertices);
    convex_nodes[0] = first_edge.from_vertex.index;
    convex_nodes[1] = first_edge.to_vertex.index;
}


void _find_convex_hull( double* lon, int len_lon,
                        double* lat, int len_lat,                          
                        int* convex_nodes, int len_nodes) {
    std::string file_name = "data/santa_barbara.csv";
    Vertices vertices1 = read_ecef_csv(file_name);
    Vertices vertices2 = make_vertices(lon, lat, len_lat);
    
    //ConvexHull convex_hull1 = find_convex_hull(vertices1);   
    ConvexHull convex_hull2 = find_convex_hull(vertices2);   
    //std::cout << convex_hull1.edges.edges.size() << " "<<  convex_hull2.edges.edges.size()<< std::endl;
    int i =0;
    for(auto const& edge: convex_hull2.edges.edges) {        
        convex_nodes[i] = edge.from_vertex.index;
        i ++;        
    }    
    
        
}


void _intersects_convex(double* lon_points, int len_lon_points, 
                        double* lat_points, int len_lat_points,                  
                        double* gc_x, int len_x, 
                        double* gc_y, int len_y, 
                        double* gc_z, int len_z, 
                        int* intersects, int len_intersects) {
    Vertices points = make_vertices(lon_points, lat_points, len_lon_points);
    ECEFs gcs = make_ecefs(gc_x, gc_y, gc_z, len_x);
    for (int i = 0; i < len_lon_points; i++) {
        for (int j = 0; j < len_x; j++) { 
            if (dot(gcs[j], points[i]) < 0) {
                intersects[i] = 0;
                break;
            }            
        }        
    }
}




void _intersects(double* lon_points, int len_lon_points, 
                 double* lat_points, int len_lat_points,                  
                 double* lon_nodes, int len_lon_nodes, 
                 double* lat_nodes, int len_lat_nodes, 
                 int* intersects, int len_intersects) {
    Vertices points = make_vertices(lon_points, lat_points, len_lon_points);              
    Vertices nodes = make_vertices(lon_nodes, lat_nodes, len_lon_nodes);    
    
    ECEFs edge_gcs;
    ECEFs left_gcs;
    ECEFs right_gcs;
    for (int j = 0; j < len_lat_nodes; j++) {
        ECEF edge_gc = cross(nodes[j], nodes[j+1]);
        ECEF left_gc = cross(edge_gc, nodes[j]);
        ECEF right_gc = cross(nodes[j+1], edge_gc);
        
        edge_gcs.push_back(edge_gc);        
        left_gcs.push_back(left_gc);
        right_gcs.push_back(right_gc);
    }


    ECEF ray_dest = ECEF(0, 0, 1);    
    
    for (int i = 0; i< len_lat_points; i++) {
        int cross_count = 0;        
        ECEF ray = cross(points[i], ray_dest);
        
        for (int j = 0; j < len_lat_nodes; j++) {
            ECEF intersection = cross(edge_gcs[j], ray);
            if ((dot(intersection, left_gcs[j]) * dot(intersection, right_gcs[j])) > 0) {                
                if (dot(points[i], edge_gcs[j]) < 0) {
                    cross_count += 1;
                } else {
                    cross_count -= 1;
                }
            }
        }
        if (cross_count > 0) intersects[i] = 1;
    }
    
}


int main(int argc, char *argv[]) {
    std::string file_name = "contrib/data/santa_barbara.csv";
    Vertices vertices = read_ecef_csv(file_name);
    ConvexHull convex_hull = find_convex_hull(vertices);    
    for(auto const& edge: convex_hull.edges.edges) {
        std::cout << edge.from_vertex.index << std::endl;        
    }    
        
}
