#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <math.h>       


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








bool test_edge(size_t from_index, size_t to_index, ECEFVector &nodes) {
    ECEF great_circle = cross(nodes[from_index], nodes[to_index]);
    for (size_t node_index = 0; node_index < nodes.size(); node_index++){
        if (node_index == from_index) continue;
        if (node_index == to_index) continue;
        if (dot(great_circle, nodes[node_index]) < 0.0) return false;
    }
    return true;
}


int find_to_index(size_t from_index, ECEFVector &nodes) {
    for (size_t to_index = 0; to_index < nodes.size(); to_index++) {
        if (from_index == to_index) continue;
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



EdgeVector _find_convex_hull(ECEFVector &nodes) {
    int first_from_index = find_first_from_index(nodes);
    int from_index = first_from_index;
    int to_index = find_to_index(from_index, nodes);
    EdgeVector convex_edges;
    while (to_index != first_from_index) {
        to_index = find_to_index(from_index, nodes);
        convex_edges.push_back(Edge(from_index, to_index));
        from_index = to_index;
    }
    return convex_edges;
}


void _find_convex_hull (double* x, int len_x, double* y, int len_y, double* z, int len_z, int* convex_nodes, int len_nodes) {       
        ECEFVector ecef_nodes = make_ecefs(x, y, z, len_x);
        EdgeVector convex_edges = _find_convex_hull(ecef_nodes);
        int i = 0;
        for (auto const& convex_edge: convex_edges) {
            convex_nodes[i] = convex_edge.from_node;
            i += 1;
        }
}

void _find_convex_hull (double* lat, int len_lat, double* lon, int len_lon, int* convex_nodes, int len_nodes) {
    ECEFVector ecef_nodes = make_ecefs(lat, lon, len_lat);
    EdgeVector convex_edges = _find_convex_hull(ecef_nodes);
    int i = 0;
    for (auto const& convex_edge: convex_edges) {
        convex_nodes[i] = convex_edge.from_node;
        i += 1;
    }
}




int main(int argc, char *argv[]) {
    ECEFVector nodes;    
    std::string file_name = "data/trinidad_sorted.csv";
    nodes = read_ecef_csv(file_name);
    EdgeVector convex_edges = _find_convex_hull(nodes);    
    for(auto const& convex_edge: convex_edges) {
        std::cout << convex_edge.from_node << " " << convex_edge.to_node << std::endl;
    }
}
