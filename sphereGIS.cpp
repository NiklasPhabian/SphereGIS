#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <math.h>       
#include <set>

struct ECEF { 
    double x, y, z;
    ECEF(double x, double y, double z): x(x), y(y), z(z) {}
};
typedef std::vector<ECEF> ECEFVector;



struct Edge {
    int from_node, to_node;
    Edge(int from_node, int to_node): from_node(from_node), to_node(to_node) {}
};
typedef std::vector<Edge> EdgeVector;



ECEFVector read_ecef_csv(std::string file_path) {
    ECEFVector cornerVector;
    std::ifstream file(file_path);
    std::string row; 
    while (std::getline(file, row)) {
        int first_comma = row.find(",");
        int second_comma = row.find(",", row.find(",")+1);
        double x = stof(row.substr(0, first_comma));        
        double y = stof(row.substr(first_comma+1));
        double z = stof(row.substr(second_comma+1));
        cornerVector.push_back(ECEF(x, y, z));       
    }
    file.close();    
    return cornerVector;
}


ECEF cross(ECEF & v1, ECEF & v2) {
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
}


bool test_edge2(size_t from_index, size_t to_index, ECEFVector &nodes) {
    if (equals(nodes[from_index], nodes[to_index])) return false;        
    ECEF great_circle = cross(nodes[from_index], nodes[to_index]);
    for (size_t node_index = 0; node_index < nodes.size(); node_index++){   
        if (node_index == from_index) continue;
        if (node_index == to_index) continue;
        if (equals(nodes[node_index], nodes[to_index])) continue;
        if (equals(nodes[node_index], nodes[from_index])) continue;
        if (dot(great_circle, nodes[node_index]) < 0) return false;
    }   
    return true;
}

std::set<size_t> offenders;

bool test_edge(size_t from_index, size_t to_index, ECEFVector &nodes) {
    if (equals(nodes[from_index], nodes[to_index])) return false;        
    ECEF great_circle = cross(nodes[from_index], nodes[to_index]);
    for (auto const& offender: offenders) {
        if (offender == from_index) continue;
        if (offender == to_index) continue;
        if (equals(nodes[offender], nodes[to_index])) continue;
        if (equals(nodes[offender], nodes[from_index])) continue;
        if (dot(great_circle, nodes[offender]) < 0) return false;
    }
    
    size_t node_index;
    for (size_t i = 0; i< nodes.size(); i++){    // This is where it is at    
        node_index = i;// + 0) % n_nodes;
        if (node_index == from_index) continue;
        if (node_index == to_index) continue;
        if (equals(nodes[node_index], nodes[to_index])) continue;
        if (equals(nodes[node_index], nodes[from_index])) continue;
        if (dot(great_circle, nodes[node_index]) < 0) {
            offenders.insert(node_index);
            return false;
        }
    }   
    return true;
}



int find_to_index2(size_t from_index, ECEFVector &nodes) {
    for (size_t to_index = 0; to_index < nodes.size(); to_index++) {
        if (from_index == to_index) continue;
        if (test_edge(from_index, to_index, nodes)) { 
            return to_index;
        }
    }
    return -2;
}


int find_to_index(size_t from_index, ECEFVector &nodes) {
    int to_index = from_index;
    int sign = 1;
    sign *= -1;
    int n_nodes =  nodes.size();
    //std::cout << to_index << std::endl;
    for (int i = 1; i < n_nodes; i++) {      
        // We move from the from_index left and right in index space
        sign *= -1;        
        to_index = (to_index + i * sign) % n_nodes;
        if (to_index < 0) {
            to_index = n_nodes + to_index;            
        }
        if (from_index == 9936){ //18468
            //std::cout << to_index << std::endl;
        }


        if (test_edge(from_index, to_index, nodes)) { 
            return to_index;
        }
    }
    return -1;
}


int find_first_from_index(ECEFVector &nodes) {
    for (size_t from_index = 0; from_index < nodes.size(); from_index++){
        if (find_to_index(from_index, nodes) >= 0) {
            return from_index;
        }
    }
    return -1;
}


ECEFVector make_ecefs(double* x, double* y, double* z, int len) {
    ECEFVector ecef_nodes;        
    for (int i = 0; i < len; i++) {      
        ecef_nodes.push_back(ECEF(x[i], y[i], z[i]));   
    }
    return ecef_nodes;
}

ECEFVector make_ecefs(double* lat, double* lon, int len) {
    ECEFVector ecef_nodes;
    for (int i = 0; i < len; i++) {      
        double x = cos(lon[i]/360*M_PI*2) * cos(lat[i]/360*M_PI*2);
        double y = sin(lon[i]/360*M_PI*2) * cos(lat[i]/360*M_PI*2);
        double z = sin(lat[i]/360*M_PI*2);
        ecef_nodes.push_back(ECEF(x, y, z));
    }
    return ecef_nodes;
}



EdgeVector find_convex_hull(ECEFVector &nodes) {
    EdgeVector convex_edges;
    int first_from_index = find_first_from_index(nodes);
    int from_index = first_from_index;
    int to_index = find_to_index(from_index, nodes);
    while (to_index != first_from_index) {
        to_index = find_to_index(from_index, nodes);
        convex_edges.push_back(Edge(from_index, to_index));
        from_index = to_index;
    }
    return convex_edges;
}


void _find_convex_hull_xyz (double* x, int len_x, 
                        double* y, int len_y, 
                        double* z, int len_z, 
                        int* convex_nodes, int len_nodes) {       

        ECEFVector nodes = make_ecefs(x, y, z, len_x);
       
        EdgeVector convex_edges = find_convex_hull(nodes);
        

        int i = 0;
        for (auto const& convex_edge: convex_edges) {
            convex_nodes[i] = convex_edge.from_node;
            i += 1;
        }
}

void _find_convex_hull (double* lon, int len_lon,
                        double* lat, int len_lat,                          
                        int* convex_nodes, int len_nodes) {
    ECEFVector ecef_nodes = make_ecefs(lat, lon, len_lat);
    EdgeVector convex_edges = find_convex_hull(ecef_nodes);
    int i = 0;
    for (auto const& convex_edge: convex_edges) {
        convex_nodes[i] = convex_edge.from_node;
        i++;
    }
}


void _intersects_convex(double* lon, int len_lon, 
                 double* lat, int len_lat,                  
                 double* x, int len_x, 
                 double* y, int len_y, 
                 double* z, int len_z, 
                 int* intersects, int len_intersects) {
    ECEFVector points = make_ecefs(lat, lon, len_lat);
    ECEFVector gcs = make_ecefs(x, y, z, len_x);
    for (int i = 0; i < len_lat; i++) {
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
    ECEFVector points = make_ecefs(lon_points, lat_points, len_lon_points);              
    ECEFVector nodes = make_ecefs(lon_nodes, lat_nodes, len_lon_nodes);    
    
    ECEFVector edge_gcs;
    ECEFVector left_gcs;
    ECEFVector right_gcs;
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
    //std::string file_name = "data/germany.csv";
    std::string file_name = "data/brazil.csv";
    ECEFVector nodes = read_ecef_csv(file_name);
    EdgeVector convex_edges = find_convex_hull(nodes);    
    for(auto const& convex_edge: convex_edges) {
        std::cout << convex_edge.from_node << " " << convex_edge.to_node << std::endl;
    }
    
        
}
